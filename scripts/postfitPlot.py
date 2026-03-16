import CombineHarvester.CombineTools.plotting as plot
import ROOT
import argparse
import math
import ctypes
from scipy.stats import chi2
import os

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)

def getHistogram(fname, histname, dirname='', allowEmpty=False):
    outname = fname.GetName()
    for key in fname.GetListOfKeys():
        histo = fname.Get(key.GetName())
        dircheck = False
        if dirname == '' : dircheck=True
        elif dirname in key.GetName(): dircheck=True
        if isinstance(histo,ROOT.TH1F) and key.GetName()==histname:
            return [histo,outname]
        elif isinstance(histo,ROOT.TDirectory) and dircheck:
            return getHistogram(histo,histname, allowEmpty=allowEmpty)
    print('Failed to find histogram with name %(histname)s in file %(fname)s '%vars())
    if allowEmpty:
        return [ROOT.TH1F('empty', '', 1, 0, 1), outname]
    else:
        return None

def backgroundComp(leg,plots,colour):
    return dict([('leg_text',leg),('plot_list',plots),('colour',colour)])

def createAxisHists(n,src,xmin=0,xmax=499):
    result = []
    for i in range(0,n):
        res = src.Clone()
        res.Reset()
        res.SetTitle("")
        res.SetName("axis%(i)d"%vars())
        res.SetAxisRange(xmin,xmax)
        res.SetStats(0)
        result.append(res)
    return result

def PositionedLegendUnrolled(width, height, pos, offset):
    o = offset
    w = width
    h = height
    l = ROOT.gPad.GetLeftMargin()
    t = ROOT.gPad.GetTopMargin()
    b = ROOT.gPad.GetBottomMargin()
    r = ROOT.gPad.GetRightMargin()
    if pos == 1:
        return ROOT.TLegend(l + o, 1 - t - o - h, l + o + w, 1 - t - o, '', 'NBNDC')
    if pos == 2:
        c = l + 0.5 * (1 - l - r)
        return ROOT.TLegend(c - 0.5 * w, 1 - t - o - h, c + 0.5 * w, 1 - t - o, '', 'NBNDC')
    if pos == 3:
        return ROOT.TLegend(1 - r - o - w, 1 - t - o - h, 1 - r - o, 1 - t - o, '', 'NBNDC')
    if pos == 4:
        return ROOT.TLegend(l + o, b + o, l + o + w, b + o + h, '', 'NBNDC')
    if pos == 5:
        c = l + 0.5 * (1 - l - r)
        return ROOT.TLegend(c - 0.5 * w, b + o, c + 0.5 * w, b + o + h, '', 'NBNDC')
    if pos == 6:
        return ROOT.TLegend(1 - r - o - w, b + o, 1 - r - o, b + o + h, '', 'NBNDC')
    if pos == 7:
        return ROOT.TLegend(1 - o - w, 1 - t - o - h, 1 - o, 1 - t - o, '', 'NBNDC')

def DrawTitleUnrolled(pad, text, align, scale=1):
    pad_backup = ROOT.gPad
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()

    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if pad_ratio < 1.:
        pad_ratio = 1.

    textSize = 0.6
    textOffset = 0.2

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(textSize * t * pad_ratio * scale)

    y_off = 1 - t + textOffset * t + 0.01
    if align == 1:
        latex.SetTextAlign(11)
        latex.DrawLatex(l, y_off, text)
    if align == 2:
        latex.SetTextAlign(21)
        latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text)
    if align == 3:
        latex.SetTextAlign(31)
        latex.DrawLatex(1 - r, y_off, text)
    pad_backup.cd()

def garwood_interval(n, cl=0.683):
    """
    Compute Garwood interval (exact Poisson confidence interval).
    Returns (low, high) for observed count n.
    """
    alpha = 1 - cl
    if n == 0:
        low = 0.0
    else:
        low = 0.5 * chi2.ppf(alpha / 2, 2 * n)
    high = 0.5 * chi2.ppf(1 - alpha / 2, 2 * (n + 1))
    return low, high

def MakeRatioGraph(num,denom):
    # make a ratio graph from a TGraphAsymmErrors and a TH1
    ratio = num.Clone()
    for i in range(0, ratio.GetN()):
        x = ctypes.c_double(0.)
        y = ctypes.c_double(0.)
        ratio.GetPoint(i, x, y)
        x = x.value
        y = y.value
        bin = denom.FindBin(x)
        denom_y = denom.GetBinContent(bin)
        if denom_y != 0:
            ratio.SetPoint(i, x, y/denom_y)
            # keep the relative errors the same
            y_high = ratio.GetErrorYhigh(i)/denom_y
            y_low = ratio.GetErrorYlow(i)/denom_y
            ratio.SetPointError(i, ratio.GetErrorXlow(i), ratio.GetErrorXhigh(i), y_low, y_high)
        else:
            ratio.SetPoint(i, x, 0)
            ratio.SetPointError(i, 0., 0., 0., 0.)
    return ratio

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Post-fit plot script for Htautau CP analysis')
    parser.add_argument('--input_file', '-i',  type=str, default='shapes_output.root',
                        help='Input ROOT file containing histograms')
    parser.add_argument('--input_file_ps', type=str,
                        help='Input ROOT file containing pseudoscalar histograms (optional)')
    parser.add_argument('--bin_name', '-b', type=str, default='htt_tt_3_13p6TeV',
                        help='Name of the bin to plot')
    parser.add_argument('--autoblind', type=float, default=0.05,
                        help='Autoblinding threshold for signal/background ratio')
    parser.add_argument('--norm_bins', action='store_true',
                        help='Normalize bins to bin width')
    parser.add_argument('--lumi', type=str, default='62.4 fb^{-1} (13.6 TeV)',
                        help='Luminosity label for the plot')
    parser.add_argument('--dir', type=str, default="",
                        help='Output directory for the plots')
    parser.add_argument('--prefit', action='store_true',
                        help='prefit_plot')
    parser.add_argument('--x_title', default=None,
                        help='Custom x-axis title (optional)')

    args = parser.parse_args()

    # for unrolled 2D plots this map will determine where the vertical lines are drawn to separate the BDT bins
    # so the numbers should equal the number of phiCP bins in the unrolled plot
    nxbins_map = {
        'htt_et_3_13p6TeV': 10,
        'htt_et_4_13p6TeV': 8,
        'htt_et_5_13p6TeV': 10,
        'htt_et_6_13p6TeV': 8,
        'htt_mt_3_13p6TeV': 10,
        'htt_mt_4_13p6TeV': 8,
        'htt_mt_5_13p6TeV': 10,
        'htt_mt_6_13p6TeV': 8,
        'htt_tt_3_13p6TeV': 10,
        'htt_tt_4_13p6TeV': 10,
        'htt_tt_5_13p6TeV': 10,
        'htt_tt_6_13p6TeV': 4,
        'htt_tt_7_13p6TeV': 10,
        'htt_tt_8_13p6TeV': 4,
        'htt_tt_9_13p6TeV': 4,
        'htt_tt_10_13p6TeV': 4,
        'htt_tt_11_13p6TeV': 4,
    }

    # these can eventually be set by command line arguments:
    lumi = args.lumi
    bin_name = args.bin_name
    channel = bin_name[4:6]


    if 'tt_6' in bin_name:
        ch_label = 'a_{1}^{3pr}a_{1}^{3pr}'
    elif 'tt_11' in bin_name:
        ch_label = 'a_{1}^{1pr}a_{1}^{3pr}'
    elif 'tt_10_' in bin_name:
        ch_label = '#pi a_{1}^{1pr}'
    elif 'tt_7' in bin_name:
        ch_label = '#pi#rho'
    elif 'tt_4' in bin_name:
        ch_label = '#rho a_{1}^{1pr}/a_{1}^{1pr}a_{1}^{1pr}'
    elif 'tt_9' in bin_name:
        ch_label = '#pi a_{1}^{3pr}'
    elif 'tt_8' in bin_name:
        ch_label = '#pi#pi'
    elif 'tt_3' in bin_name:
        ch_label = '#rho#rho'
    elif 'tt_5' in bin_name:
        ch_label = '#rho a_{1}^{3pr}'
    elif 'et_6' in bin_name:
        ch_label = 'ea_{1}^{1pr}'
    elif 'et_3' in bin_name:
        ch_label = 'e#rho'
    elif 'et_5' in bin_name:
        ch_label = 'ea_{1}^{3pr}'
    elif 'et_4' in bin_name:
        ch_label = 'e#pi'
    elif 'mt_6' in bin_name:
        ch_label = '#mua_{1}^{1pr}'
    elif 'mt_3' in bin_name:
        ch_label = '#mu#rho'
    elif 'mt_5' in bin_name:
        ch_label = '#mua_{1}^{3pr}'
    elif 'mt_4' in bin_name:
        ch_label = '#mu#pi'
    elif 'tt_1' in bin_name or 'tt_2' in bin_name or 'tt_100' in bin_name:
        ch_label = '#tau_{h}#tau_{h}'
    elif 'mt_1' in bin_name or 'mt_2' in bin_name or 'mt_100' in bin_name:
        ch_label = '#mu#tau_{h}'
    elif 'et_1' in bin_name or 'et_2' in bin_name or 'et_100' in bin_name:
        ch_label = 'e#tau_{h}'

    if '_tt_8' in bin_name:
        ch_bins = ['(BDT 0.33-0.7)', '(BDT 0.7-0.8)', '(BDT 0.8-1.0)']
    elif '_tt_' in bin_name:
        ch_bins = ['(BDT 0.33-0.7)', '(BDT 0.7-0.8)', '(BDT 0.8-0.9)', '(BDT 0.9-1.0)']
    elif '_mt_' in bin_name:
        ch_bins = ['(BDT 0.33-0.6)', '(BDT 0.6-0.7)', '(BDT 0.7-0.8)', '(BDT 0.8-0.9)', '(BDT 0.9-1.0)']
    elif '_et_' in bin_name:
        ch_bins = ['(BDT 0.33-0.6)', '(BDT 0.6-0.8)', '(BDT 0.8-1.0)']


    input_file = ROOT.TFile(args.input_file)
   
    if args.input_file_ps:
        input_file_ps = ROOT.TFile(args.input_file_ps)
    else: input_file_ps = None
    autoblind = args.autoblind
    norm_bins = args.norm_bins

    is_1d_bin = bin_name in ['htt_mt_1_13p6TeV', 'htt_mt_2_13p6TeV', 'htt_et_1_13p6TeV', 'htt_et_2_13p6TeV', 'htt_tt_1_13p6TeV', 'htt_tt_2_13p6TeV',
    'htt_tt_100_13p6TeV','htt_mt_100_13p6TeV','htt_et_100_13p6TeV']

    is_control_bin = bin_name in ['htt_tt_100_13p6TeV','htt_mt_100_13p6TeV','htt_et_100_13p6TeV']


    signal_scale = 1.0
    if is_control_bin:
        signal_scale = 20.0
    
    if norm_bins:
        y_title = 'Events / bin width'
    else: 
        y_title = 'Events'
    
    if is_1d_bin:
        ratio_range = '0.75,1.25'
    else:
        ratio_range = '0.0,2.0'


    if is_1d_bin:
        plot.ModTDRStyle(width=900, height=800, r=0.3, l=0.15)
        if 't_1' in bin_name:
            x_title = 'BDT score (Tau)'
        elif 't_2' in bin_name:
            x_title = 'BDT score (Fake)'
        logy = False
        extra_pad = 0.15
    else:
        plot.ModTDRStyle(width=1800, height=700, r=0.3, l=0.16, t=0.12,b=0.15)
        x_title = 'Bin number'
        logy = True
        extra_pad = 0.5

    if args.x_title:
        x_title = args.x_title
    
    canv = ROOT.TCanvas()    
    pads=plot.TwoPadSplit(0.35,0.01,0.01)
    
    bkghist = getHistogram(input_file,'TotalProcs',bin_name)[0]
    bkgonlyhist = getHistogram(input_file,'TotalBkg',bin_name)[0]
    
    sighist = getHistogram(input_file,'TotalSig',bin_name)[0]
    
    # if pseudoscalar file is provided then prdoduce a signal histogram from it and average it with the scalar one for the autoblinding
    # if it isn't provided then just use the scalar one
    if input_file_ps:
        sighist_ps = getHistogram(input_file_ps,'TotalSig',bin_name)[0]
        ave_sig = sighist.Clone()
        ave_sig.Add(sighist_ps)
        ave_sig.Scale(0.5)
    else :
        sighist_ps = None
        ave_sig = sighist.Clone()
    
    datahist = getHistogram(input_file,"data_obs",bin_name)[0]
    datahist.SetMarkerStyle(20)
    
    if autoblind>0:
        # get rid of any data bins when the ave_sig/bkghist fraction is > autoblind by settinfg the yield to -0.1 and the error to 0
        for i in range(1, datahist.GetNbinsX()+1):
            if bkghist.GetBinContent(i) == 0: blind = True
            elif ave_sig.GetBinContent(i) / bkghist.GetBinContent(i) > autoblind: blind = True
            else: blind = False
            if blind:
                datahist.SetBinContent(i, -0.1)
                datahist.SetBinError(i, 0)


    #convert datahist into datagraph using TGraphAsymmErrors for proper Poissonian error bars
    datagraph = ROOT.TGraphAsymmErrors(datahist)
    datagraph.SetName("data_obs_graph")

    for i in range(0, datagraph.GetN()):
        x = ctypes.c_double(0.)
        y = ctypes.c_double(0.)
        datagraph.GetPoint(i, x, y)
        x = x.value
        y = y.value
        x_err_low = datagraph.GetErrorXlow(i)
        x_err_high = datagraph.GetErrorXhigh(i)
        if y < 0:
            datagraph.SetPoint(i, x, 0)
            datagraph.SetPointError(i, x_err_low, x_err_high, 0., 0.)
        else:
            low, high = garwood_interval(int(y), cl=0.683)
            datagraph.SetPointError(i, x_err_low, x_err_high, y - low, high - y)
    datahist = datagraph
    
    bkghist.SetMarkerSize(0)
    bkghist.SetFillColor(2001)
    bkghist.SetLineColor(0)
    bkghist.SetLineWidth(1)
    bkghist.SetFillColor(plot.CreateTransparentColor(12,0.4))
    bkgonlyhist.SetMarkerSize(0)
    bkgonlyhist.SetFillColor(2001)
    bkgonlyhist.SetLineColor(0)
    bkgonlyhist.SetLineWidth(1)
    bkgonlyhist.SetFillColor(plot.CreateTransparentColor(12,0.4))
    if norm_bins: 
        bkghist.Scale(1.0,"width")
        bkgonlyhist.Scale(1.0,"width")
        sighist.Scale(1.0,"width")
        datahist.Scale(1.0,"width")
    
    if args.prefit:
        h_label = "H#rightarrow#tau#tau (#alpha^{H#tau#tau}=0^{#circ})"
        if is_control_bin:
            h_label = "H#rightarrow#tau#tau"
        background_schemes = {
            'et':[
                    backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)), ###TODO changed this to TT once the template is named correctly!!!
                    backgroundComp("Jet#rightarrow#tau_{h}",["JetFakes"],ROOT.TColor.GetColor(160,193,114)),
                    backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(255,169,14)),
                    backgroundComp(f"{h_label}",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
                    ],
            'mt':[
                    backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)), ###TODO changed this to TT once the template is named correctly!!!
                    backgroundComp("Jet#rightarrow#tau_{h}",["JetFakes"],ROOT.TColor.GetColor(160,193,114)),
                    backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(255,169,14)),
                    backgroundComp(f"{h_label}",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
                    ],
            'tt':[
                    backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)),
                    backgroundComp("Jet#rightarrow#tau_{h}",['JetFakes','JetFakesSublead'],ROOT.TColor.GetColor(160,193,114)),
                    backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(255,169,14)),
                    backgroundComp(f"{h_label}",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
                    ],
        }
    else:
        background_schemes = {
            'et':[
                    backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)), ###TODO changed this to TT once the template is named correctly!!!
                    backgroundComp("Jet#rightarrow#tau_{h} fakes",["JetFakes"],ROOT.TColor.GetColor(160,193,114)),
                    backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(255,169,14)),
                    backgroundComp("H#rightarrow#tau#tau (Best Fit)",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
                    ],
            'mt':[
                    backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)), ###TODO changed this to TT once the template is named correctly!!!
                    backgroundComp("Jet#rightarrow#tau_{h} fakes",["JetFakes"],ROOT.TColor.GetColor(160,193,114)),
                    backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(255,169,14)),
                    backgroundComp("H#rightarrow#tau#tau (Best Fit)",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
                    ],
            'tt':[
                    backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)),
                    backgroundComp("Jet#rightarrow#tau_{h} fakes",['JetFakes','JetFakesSublead'],ROOT.TColor.GetColor(160,193,114)),
                    backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(255,169,14)),
                    backgroundComp("H#rightarrow#tau#tau (Best Fit)",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
                    ],
        }

    
    stack_histos = []
    
    for t in background_schemes[channel]:
        plots = t['plot_list']
        isHist = False
        h = ROOT.TH1F()
        for k in plots:
            if h.GetEntries()==0 and getHistogram(input_file,k, bin_name,False) is not None:
                isHist = True
                h = getHistogram(input_file,k, bin_name)[0]
                h.SetName(k)
            else:
                if getHistogram(input_file,k, bin_name,False) is not None:
                    isHist = True
                    h.Add(getHistogram(input_file,k, bin_name)[0])
    
        h.SetFillColor(t['colour'])
        h.SetLineColor(ROOT.kBlack)
        h.SetMarkerSize(0)
        
        if norm_bins: h.Scale(1.0,"width")
        if isHist:
            stack_histos.append(h)
    
    stack = ROOT.THStack("hs","")
    for hists in stack_histos:
        stack.Add(hists)
    
    signal_stack = ROOT.THStack("hs_sig","")
    if args.prefit:
        if sighist.GetEntries() > 0:
            sighist.SetFillStyle(0)
            sighist.SetLineColor(ROOT.kRed)
            sighist.SetLineWidth(2)
            if norm_bins: sighist.Scale(1.0,"width")
            sighist.Scale(signal_scale)
            signal_stack.Add(sighist)


    if sighist_ps and sighist_ps.GetEntries() > 0:
        sighist_ps.SetFillStyle(0)
        sighist_ps.SetLineColor(ROOT.kBlue)
        sighist_ps.SetLineWidth(2)
        if norm_bins: sighist_ps.Scale(1.0,"width")
        sighist_ps.Scale(signal_scale)
        signal_stack.Add(sighist_ps)
    
    x_min = bkghist.GetXaxis().GetXmin()
    if is_control_bin:
        # check if any of the first bins are empty and adjuts x_min if so
        for i in range(1, bkghist.GetNbinsX()+1):
            if bkghist.GetBinContent(i) > 10**(-8):
                x_min = bkghist.GetXaxis().GetBinLowEdge(i)
                break
    axish = createAxisHists(2,bkghist,x_min,bkghist.GetXaxis().GetXmax()-0.01)
    if not is_1d_bin:
        axish[0].SetMaximum(3*bkghist.GetMaximum())
    else:
        axish[0].SetMaximum(bkghist.GetMaximum())
    # axish[0].SetMaximum(bkghist.GetMaximum())
    axish[0].GetYaxis().SetTitle(y_title)
    axish[0].GetXaxis().SetLabelSize(0)
    axish[0].GetXaxis().SetTitleSize(0)
    axish[1].GetXaxis().SetTitle(x_title)
    if not is_1d_bin:
        axish[1].GetXaxis().SetTitleOffset(1.0)
    
    
    pads[0].cd()
    pads[0].SetTicks(1)
    if logy:
        pads[0].SetLogy()
        sig_min = sighist.GetMinimum()
        if sig_min > 0:
            axish[0].SetMinimum(0.5*sig_min)
        else:
            axish[0].SetMinimum(0.01)
    
        #axish[0].SetMinimum(0.01)
    axish[0].Draw("hist")
    
    stack.Draw("histsame")
    bkghist.Draw("e2same")
    signal_stack.Draw("histsamenostack")
    datahist.Draw("P Z 0 same")

    pads[0].RedrawAxis()
    
    # now we produce the ratio plot
    # For now the ratio just shows the ratio of the observed data to the signal+background expectation
    # this can be changed in the future to show the ratio to background only and signal+background etc seperatly 
    
    pads[1].cd()
    pads[1].SetTicks(1)
    pads[1].SetGrid(0,1)
    axish[1].Draw("axis")
    axish[1].SetMinimum(float(ratio_range.split(',')[0]))
    axish[1].SetMaximum(float(ratio_range.split(',')[1]))
    axish[1].GetYaxis().SetTitle("Data/Exp.")
    # draw dashed line at 1
    line = ROOT.TLine(x_min, 1, axish[1].GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kBlack)
    line.Draw("same")

    # note bkgonlyhist has issues with uncertainties right now (large as uncertainty on signal is large)
    
    ratio_bkghist = plot.MakeRatioHist(bkghist,bkghist,True,False)
    ratio_bkghist.Draw("e2same")
    if isinstance(datahist,ROOT.TGraphAsymmErrors):
        ratio_datahist = MakeRatioGraph(datahist,bkghist)
    else:
        ratio_datahist = plot.MakeRatioHist(datahist,bkghist,True,False)
    ratio_datahist.Draw("P Z 0 same")

    #sbhist = bkghist.Clone()
    #sbhist.SetLineColor(ROOT.kRed)
    #sbhist.SetLineWidth(2)
    #sbhist.SetFillStyle(0)
    #ratio_sighist = plot.MakeRatioHist(sbhist,bkgonlyhist,True,False)
    #ratio_sighist.Draw("histsame")
    pads[1].RedrawAxis()

    Nxbins = nxbins_map[bin_name] if bin_name in nxbins_map else None
    if not is_1d_bin and Nxbins:
        # add lines to seperate BDT bins
        yline = ROOT.TLine()
        yline.SetLineWidth(2)
        yline.SetLineStyle(3)
        yline.SetLineColor(1)
        x = int(bkghist.GetNbinsX()/Nxbins)
        for l in range(1,x):
            pads[0].cd()
            ymax = axish[0].GetMaximum()
            ymin = axish[0].GetMinimum()
            fraction = extra_pad if extra_pad>0 else 0.15
            ymax = (ymax - fraction * ymin) / (1. - fraction)
            yline.DrawLine(l*Nxbins,ymin,l*Nxbins,ymax)

            pads[1].cd()
            ymax = axish[1].GetMaximum()
            ymin = axish[1].GetMinimum()
            yline.DrawLine(l*Nxbins,ymin,l*Nxbins,ymax) 

        # add extra axis for the phiCP angle
        scale = abs(float(ratio_range.split(',')[0])-float(ratio_range.split(',')[1]))/2
        y_axis = float(ratio_range.split(',')[0]) - 0.7*scale
        extra_axis = ROOT.TGaxis(0,y_axis,Nxbins,y_axis,0,2*math.pi,400,"NS")
        extra_axis.SetLabelSize(0)
        #extra_axis.SetLabelFont(42)
        extra_axis.SetMaxDigits(2)
        extra_axis.SetTitle("#phi_{CP}")
        extra_axis.SetTitleFont(42)
        extra_axis.SetTitleSize(0.035)
        extra_axis.SetTitleOffset(1.05)
        extra_axis.SetTickSize(0.08)
        extra_axis.Draw()  

        positions = [0, math.pi/2, math.pi, 3*math.pi/2, 2*math.pi]
        labels = ["0", "#pi/2", "#pi", "3#pi/2", "2#pi"]
        
        # Use TLatex to manually draw labels
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(21)

        # Map axis positions to pad x-coordinates
        for i, xval in enumerate(positions):
            x_pixel = (xval) / (2 * math.pi) * Nxbins  # from axis range to pixel
            latex.DrawLatex(x_pixel, y_axis - 0.3 * scale, labels[i])  # adjust vertical offset as needed

        latex.SetNDC()
        latex.SetTextSize(0.04)

        for ib, b in enumerate(ch_bins):
            latex.DrawLatex((ib+0.83)/(len(ch_bins)+1), 0.82, "#font[42]{" + str(b) + "}")
    # draw the title, labels, and legend
    
    #CMS and lumi labels
    if not is_1d_bin:
        plot.FixTopRange(pads[0], 3*bkghist.GetMaximum(), extra_pad if extra_pad>0 else 0.15)
    else:
        plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), extra_pad if extra_pad>0 else 0.15)
    extra='Preliminary'
    if is_1d_bin:
        plot.DrawCMSLogo(pads[0], 'CMS', extra, 0, 0.17, -0.0, 2.0, '', 0.85)
        if is_control_bin: # fraf the label inside the frame on the top left for the control region to save space
            ch_tex = ROOT.TLatex()
            ch_tex.SetNDC()
            ch_tex.SetTextFont(42)
            ch_tex.SetTextSize(0.6 * pads[0].GetTopMargin())
            ch_tex.SetTextAlign(11)
            ch_tex.DrawLatex(0.15, 0.9, ch_label)
            plot.DrawTitle(pads[0], lumi, 3, textSize=0.6)
        else: plot.DrawTitle(pads[0], ch_label + "   " + lumi, 3, textSize=0.6)
    else:
        plot.DrawCMSLogo(pads[0], 'CMS', extra, 0, 0.075, -0.0, 2.0, '', 0.6)
        DrawTitleUnrolled(pads[0], ch_label + " "*50 + lumi, 3, scale=0.7)
    if is_1d_bin:
        legend = PositionedLegendUnrolled(0.25,0.5,7,0.)
    else:
        legend = PositionedLegendUnrolled(0.11,0.5,7,0.)
    legend.SetTextSize(0.03) 
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    
    legend.AddEntry(datahist,"Data","PEl")
    # loop backwards over the stack_histos by index to add them in the right order
    for legi in range(len(stack_histos)-1,-1,-1):
        hists = stack_histos[legi]
        legend.AddEntry(hists,background_schemes[channel][legi]['leg_text'],"f")
    
    legend.AddEntry(bkghist,"Uncertainty","f")
    # comment to remove signal
    if args.prefit:
        if is_control_bin:
            legend.AddEntry(sighist,f"H#rightarrow#tau#tau #times {signal_scale:.0f}", "l")
        else: legend.AddEntry(sighist,"H#rightarrow#tau#tau (#alpha^{H#tau#tau}=0^{#circ})","l")
    if sighist_ps:
        legend.AddEntry(sighist_ps,"H#rightarrow#tau#tau (#alpha^{H#tau#tau}=90^{#circ})","l")
    legend.Draw("same")
    
    canv.Print(os.path.join(args.dir, bin_name + '.pdf'))
    canv.Print(os.path.join(args.dir, bin_name + '.png'))
    
