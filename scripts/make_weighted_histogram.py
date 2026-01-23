import CombineHarvester.CombineTools.ch as ch
import argparse
import ROOT
import math
import os

output_dir = 'weighted_phi_CP_plots'
os.makedirs(output_dir, exist_ok=True)

ROOT.gROOT.SetBatch(True)

# keep track of number of phiCP bins for each channel/category
bins_map ={
}

for i in [3,4,5,7]:
    bins_map[f'htt_tt_{i}_13p6TeV'] = 10
for i in [6,8,9,10,11]:
    bins_map[f'htt_tt_{i}_13p6TeV'] = 4
for i in [3,5]:
    bins_map[f'htt_mt_{i}_13p6TeV'] = 10
    bins_map[f'htt_et_{i}_13p6TeV'] = 10
for i in [4,6]:
    bins_map[f'htt_mt_{i}_13p6TeV'] = 8
    bins_map[f'htt_et_{i}_13p6TeV'] = 8

for year in [2016,2017,2018]:
    for i in [3,7]:
        bins_map[f'htt_tt_{year}_{i}_13TeV'] = 10
    for i in [4,5,6,9,10,11]:
        bins_map[f'htt_tt_{year}_{i}_13TeV'] = 4
    bins_map[f'htt_mt_{year}_3_13TeV'] = 10
    bins_map[f'htt_et_{year}_3_13TeV'] = 10
    bins_map[f'htt_mt_{year}_4_13TeV'] = 8
    bins_map[f'htt_et_{year}_4_13TeV'] = 8
    for i in [5,6]:
        bins_map[f'htt_mt_{year}_{i}_13TeV'] = 4
        bins_map[f'htt_et_{year}_{i}_13TeV'] = 4


# negative phase flip categories have the phase shift applies in the opposite direction
negative_phase_flip_categories = []
negative_phase_flip_categories.append('htt_tt_6_13p6TeV')
for year in [2016,2017,2018]:
    negative_phase_flip_categories.append(f'htt_tt_{year}_6_13TeV')

phase_flip_categories = []

for year in [2016,2017,2018]:
    for i in [5,6,9,11]:
        phase_flip_categories.append(f'htt_tt_{year}_{i}_13TeV')
    for i in [3,4,5,6]:
        phase_flip_categories.append(f'htt_mt_{year}_{i}_13TeV')
        phase_flip_categories.append(f'htt_et_{year}_{i}_13TeV')
for i in [5,6,9,11]:
    phase_flip_categories.append(f'htt_tt_{i}_13p6TeV')
for i in [3,4,6]:
    phase_flip_categories.append(f'htt_mt_{i}_13p6TeV')
    phase_flip_categories.append(f'htt_et_{i}_13p6TeV')

parser = argparse.ArgumentParser(description='Post-fit plot script for Htautau CP analysis')
parser.add_argument('--fitresult', '-f', help= 'Path to a RooFitResult, only needed for postfit', default=None)
parser.add_argument('--workspace', '-w', help= 'The input workspace-containing file [REQUIRED]')
parser.add_argument('--output-folder', '-o', help= 'Output folder for datacards', default='cards_weighted_histograms')
parser.add_argument('--unblind', action='store_true', help='Unblind the data, if not set the data will be set the the Asimov dataset')
parser.add_argument('--alt_weights', action='store_true', help='Use alternative weights based on A*ln(1+S/B) instead of A*S/(S+B)')
parser.add_argument('--extra_plots', action='store_true', help='Make extra plots for each category separately, and for Run-2 only')
args = parser.parse_args()

cb = ch.CombineHarvester()
infile = ROOT.TFile(args.workspace)
ws = infile.Get('w')

cb.SetFlag('workspaces-use-clone', True)
ch.ParseCombineWorkspace(cb, ws, "ModelConfig", "data_obs", False)

proto1 = ROOT.TH1F("proto1", "proto1", 10, 0,360)
proto2 = ROOT.TH1F("proto2", "proto2", 4, 0,360)
proto3 = ROOT.TH1F("proto3", "proto3", 8, 0,360)

# if fit result is given, update the parameters to postfit values
f_fit = ROOT.TFile(args.fitresult.split(':')[0])
res = f_fit.Get(args.fitresult.split(':')[1])
cb.UpdateParameters(res)

bin_set = cb.bin_set()

def apply_phase_shift(hist, nxbins):
    half = nxbins // 2

    histout = hist.Clone()

    for b in range(1, hist.GetNbinsX()+1):

        new_bin=b+half - (b+half - nxbins*((b-1) // nxbins) > nxbins)*nxbins

        histout.SetBinContent(new_bin, hist.GetBinContent(b))
        histout.SetBinError(new_bin, hist.GetBinError(b))

    return histout

def process_hist(hist, wt_mapping, nxbins, phase_shift):
    hist_mod = hist.Clone()
    if phase_shift != 0:
        hist_mod = apply_phase_shift(hist_mod, nxbins)
    
    histout = ROOT.TH1F("histout", "histout", nxbins, 0,360)
    histout.Sumw2()

    # now combine bins according to weights
    for b in range(1, hist_mod.GetNbinsX()+1):
        new_bin = (b - 1) % nxbins + 1

        wt = wt_mapping[b]

        histout.SetBinContent(
            new_bin,
            histout.GetBinContent(new_bin) + wt * hist_mod.GetBinContent(b)
        )

        histout.SetBinError(
            new_bin,
            math.sqrt(pow(histout.GetBinError(new_bin),2) + pow(wt*hist_mod.GetBinError(b),2))
        )

    if phase_shift < 0:
        n = histout.GetNbinsX()
        for i in range(1, n // 2 + 1):
            j = n - i + 1

            ci, ei = histout.GetBinContent(i), histout.GetBinError(i)
            cj, ej = histout.GetBinContent(j), histout.GetBinError(j)

            histout.SetBinContent(i, cj); histout.SetBinError(i, ej)
            histout.SetBinContent(j, ci); histout.SetBinError(j, ei)  

    if nxbins == 8: histout.Rebin(2)
    return histout

def ZeroErrors(src):
  for x in range(1,src.GetNbinsX()+1):
    src.SetBinError(x, 0.)
  return src    

samples = 500

# get the initial values of all the parameters - note while we call randomizePars we do not use these values in p_vec at this point
rands = res.randomizePars()
p_vec = [None]*len(rands)
for n in range(0,len(rands)):
    p_vec[n] = cb.cp().GetParameter(rands[n].GetName())

histograms = {}
for b in bin_set:
    histograms[b] = {'data': None, 'sm_sig': None, 'ps_sig': None, 'mm_sig': None, 'bkg': None, 'bkg_variations': []}

wt_mapping = {} # will store weights for rescaling histograms based on expected sensitivity
for samp in range(samples+1): # note samp = 0 is nominal

    if samp % 50 ==0:
        print(f"Processing sample {samp} / {samples}")

    for cat in bin_set:
        # skip if not in our bins_map
        if cat not in bins_map:
            continue
 
        if cat in negative_phase_flip_categories:
            phase_flip = -1
        elif cat in phase_flip_categories:
            phase_flip = 1
        else:
            phase_flip = 0

        nxbins = bins_map[cat]
        sel_bin = cb.cp().bin([cat])

        if samp == 0:
            wt_mapping[cat] = {}
            # we only need to get signal and data for the nominal sample as we don't care about the uncertainties for these
    
            # get SM and CP-odd signals
            # to get PS we scale alpha to 90 degrees
            par = cb.GetParameter('alpha')
            par.set_val(90.)
            ps_sig = sel_bin.cp().signals().process(['ggH_ps_prod_sm_htt','ggH_ps_htt','qqH_ps_htt','WH_ps_htt','ZH_ps_htt']).GetShape()
            # now get SM by setting alpha to 0 degrees
            par.set_val(45.)
            mm_sig = sel_bin.cp().signals().process(['ggH_mm_prod_sm_htt','ggH_mm_htt','qqH_mm_htt','WH_mm_htt','ZH_mm_htt']).GetShape()
            par.set_val(0.)
            sm_sig = sel_bin.cp().signals().process(['ggH_sm_prod_sm_htt','ggH_sm_htt','qqH_sm_htt','WH_sm_htt','ZH_sm_htt']).GetShape()
        
            bkg = sel_bin.cp().backgrounds().GetShape()
    
            if args.unblind: 
                data = sel_bin.cp().GetObservedShape()
            else:
                par.set_val(0.) # take CP-even SM hypothesis for Asimov data 
                data = sel_bin.cp().GetShape()
                for b in range(1, data.GetNbinsX()+1):
                    data.SetBinError(b, math.sqrt(data.GetBinContent(b)))
    
            # for the nominal case we also calculate the weights used to reweight the individual histograms 
            s_sb=0
            A_ave=0
            for b in range(1, bkg.GetNbinsX()+1):
                if (b-1) % nxbins ==0 and b+nxbins-1 <= bkg.GetNbinsX():
                    i_sm = sm_sig.Integral(b,b+nxbins-1)
                    i_ps = ps_sig.Integral(b,b+nxbins-1)
                    i_bkg = bkg.Integral(b,b+nxbins-1)
                    i_sig = (i_sm+i_ps)/2
                    s_vs_b_weight = i_sig/(i_sig+i_bkg)
                    if args.alt_weights:
                        s_vs_b_weight = math.log(1+i_sig/i_bkg) # ATLAS style
                    A_tot=0
                    for i in range(b, b+nxbins):
                        b_sm = sm_sig.GetBinContent(i)
                        b_ps = ps_sig.GetBinContent(i)
                        if b_sm+b_ps>1e-6:
                            A_tot += abs(b_sm-b_ps)/(b_sm+b_ps)

                    A_ave = A_tot/nxbins
                wt_mapping[cat][b] = s_vs_b_weight*A_ave

            mm_sig = process_hist(mm_sig, wt_mapping[cat], nxbins, phase_flip)
            sm_sig = process_hist(sm_sig, wt_mapping[cat], nxbins, phase_flip)
            ps_sig = process_hist(ps_sig, wt_mapping[cat], nxbins, phase_flip)
            sm_sig = ZeroErrors(sm_sig)
            ps_sig = ZeroErrors(ps_sig)
            mm_sig = ZeroErrors(mm_sig)
            data = process_hist(data, wt_mapping[cat], nxbins, phase_flip)
        
            bkg = process_hist(bkg, wt_mapping[cat], nxbins, phase_flip)
            bkg = ZeroErrors(bkg)

            histograms[cat]['data'] = data.Clone()
            histograms[cat]['sm_sig'] = sm_sig.Clone()
            histograms[cat]['ps_sig'] = ps_sig.Clone()
            histograms[cat]['mm_sig'] = mm_sig.Clone()
            histograms[cat]['bkg'] = bkg.Clone()

        else: 
            # now we sample the covariance matrix to get variations on the background only

            res.randomizePars()
            for n in range(0,len(rands)):
                if p_vec[n]: p_vec[n].set_val(rands[n].getVal())

            bkg_var = sel_bin.cp().backgrounds().GetShape()
            bkg_var = process_hist(bkg_var, wt_mapping[cat], nxbins, phase_flip)
            bkg_var = ZeroErrors(bkg_var)
            histograms[cat]['bkg_variations'].append(bkg_var.Clone())

def Subtract(h1,h2):
  for i in range(1,h1.GetNbinsX()+1):
    diff = h1.GetBinContent(i) - h2.GetBinContent(i)
    h1.SetBinContent(i,diff)
  return h1

def CombineCats(cats, histograms):
    for i, cat in enumerate(cats):
        if i == 0:
            bkg = histograms[cat]['bkg'].Clone()
            sig_sm = histograms[cat]['sm_sig'].Clone()
            sig_ps = histograms[cat]['ps_sig'].Clone()
            sig_mm = histograms[cat]['mm_sig'].Clone()
            data = histograms[cat]['data'].Clone()

            bkg_variations = []
            for var in histograms[cat]['bkg_variations']:
                bkg_variations.append(var.Clone())
        else:
           # check if number of bins is the same and if not print a warning
           if bkg.GetNbinsX() != histograms[cat]['bkg'].GetNbinsX():
                  print(f"Warning: Number of bins in category {cat} does not match. {bkg.GetNbinsX()} vs {histograms[cat]['bkg'].GetNbinsX()}")
           bkg.Add(histograms[cat]['bkg'])
           sig_sm.Add(histograms[cat]['sm_sig'])
           sig_ps.Add(histograms[cat]['ps_sig'])
           sig_mm.Add(histograms[cat]['mm_sig'])
           data.Add(histograms[cat]['data'])
           for j, var in enumerate(histograms[cat]['bkg_variations']):
                bkg_variations[j].Add(var) 

    bkg = ZeroErrors(bkg) # should be 0 already but just to make absolutely sure
    for i in range(1, bkg.GetNbinsX()+1):
        err_sq = 0.
        for bkg_var in bkg_variations:
            err_sq += (bkg_var.GetBinContent(i)-bkg.GetBinContent(i))**2
        err = (err_sq/len(bkg_variations))**.5 # divide by number of samples to get RMS then sqrt to get error
        bkg.SetBinError(i, err)

    # do bkg subraction here:
    data = Subtract(data, bkg)
    bkg = Subtract(bkg, bkg)

    return bkg, sig_sm, sig_ps, sig_mm, data

#make a list of all the categories with 10 bins
best_cats = [cat for cat in bin_set if cat in bins_map and bins_map[cat]==10]

#make a list of the categories with or 8 bins i.e not = 10
worst_cats = [cat for cat in bin_set if cat in bins_map and bins_map[cat]!=10]

bkg_best, sig_sm_best, sig_ps_best, sig_mm_best, data_best = CombineCats(best_cats, histograms)
bkg_worst, sig_sm_worst, sig_ps_worst, sig_mm_worst, data_worst = CombineCats(worst_cats, histograms)


fout = ROOT.TFile(f'{output_dir}/weighted_phiCP_histograms.root', 'RECREATE')
# make a directory for best categories
fout.mkdir('best_categories')
fout.cd('best_categories')
bkg_best.Write('bkg')
sig_sm_best.Write('sig_sm')
sig_ps_best.Write('sig_ps')
sig_mm_best.Write('sig_mm')
data_best.Write('data_obs')

# make a directory for worst categories
fout.mkdir('worst_categories')
fout.cd('worst_categories')
bkg_worst.Write('bkg')
sig_sm_worst.Write('sig_sm')
sig_ps_worst.Write('sig_ps')
sig_mm_worst.Write('sig_mm')
data_worst.Write('data_obs')

fout.Close()

import CombineHarvester.CombineTools.plotting as plot

def DrawCMSLogoCustom(pad, cmsText, extraText, iPosX, relPosX, relPosY, relExtraDY, relExtraDX, extraText2='', cmsTextSize=0.8):
    """
    
    Args:
        pad (TYPE): Description
        cmsText (TYPE): Description
        extraText (TYPE): Description
        iPosX (TYPE): Description
        relPosX (TYPE): Description
        relPosY (TYPE): Description
        relExtraDY (TYPE): Description
        extraText2 (str): Description
        cmsTextSize (float): Description
    
    Returns:
        TYPE: Description
    """
    pad.cd()
    cmsTextFont = 62  # default is helvetic-bold

    writeExtraText = len(extraText) > 0
    writeExtraText2 = len(extraText2) > 0
    extraTextFont = 52

    # text sizes and text offsets with respect to the top frame
    # in unit of the top margin size
    lumiTextOffset = 0.2
    # cmsTextSize = 0.8
    # float cmsTextOffset    = 0.1;  // only used in outOfFrame version

    # ratio of 'CMS' and extra text size
    extraOverCmsTextSize = 0.76

    outOfFrame = False
    if iPosX / 10 == 0:
        outOfFrame = True

    alignY_ = 3
    alignX_ = 2
    if (iPosX / 10 == 0):
        alignX_ = 1
    if (iPosX == 0):
        alignX_ = 1
    if (iPosX == 0):
        alignY_ = 1
    if (iPosX / 10 == 1):
        alignX_ = 1
    if (iPosX / 10 == 2):
        alignX_ = 2
    if (iPosX / 10 == 3):
        alignX_ = 3
    # if (iPosX == 0): relPosX = 0.14
    align_ = 10 * alignX_ + alignY_

    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    extraTextSize = extraOverCmsTextSize * cmsTextSize
    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if (pad_ratio < 1.):
        pad_ratio = 1.

    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText)

    posX_ = 0
    if iPosX % 10 <= 1:
        posX_ = l + relPosX * (1 - l - r)
    elif (iPosX % 10 == 2):
        posX_ = l + 0.5 * (1 - l - r)
    elif (iPosX % 10 == 3):
        posX_ = 1 - r - relPosX * (1 - l - r)

    posY_ = 1 - t - relPosY * (1 - t - b)
    if not outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, cmsText)
        if writeExtraText:
            latex.SetTextFont(extraTextFont)
            latex.SetTextAlign(align_)
            latex.SetTextSize(extraTextSize * t * pad_ratio)
            latex.DrawLatex(
                posX_+ relExtraDX * cmsTextSize * t, posY_ - relExtraDY * cmsTextSize * t, extraText)
            if writeExtraText2:
                latex.DrawLatex(
                    posX_, posY_ - 1.8 * relExtraDY * cmsTextSize * t, extraText2)
    elif writeExtraText:
        if iPosX == 0:
            posX_ = l + relPosX * (1 - l - r)
            posY_ = 1 - t + lumiTextOffset * t
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)

def propoganda_plot_phicp(sm,ps,mm, bkg,data,plot_name,extra_label='Preliminary',plot_mm=False):

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(0.06)

    data.GetXaxis().SetTitleOffset(1.0)
    data.GetYaxis().SetTitleOffset(0.8)
    data.GetXaxis().SetTitleSize(0.05)
    if args.alt_weights: 
        data.GetYaxis().SetTitle('A#kern[0.1]{ln(1+S/B)} weighted events / bin')
    else: 
        data.GetYaxis().SetTitle('A#kern[0.1]{S}/(S+B) weighted events / bin')
    data.GetXaxis().SetTitle('#phi_{#it{CP}} (degrees)')
    data.GetYaxis().SetTitleSize(0.05)
    data.GetXaxis().SetNdivisions(506,False)

    c1 = ROOT.TCanvas("c1","c1",700,500)

    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.TH1.AddDirectory(False)
    plot.ModTDRStyle(r=0.04, l=0.14)

    pads=plot.OnePad()
    pads[0].cd()

    hs = ROOT.THStack("hs","")

    data.SetMarkerStyle(20)
    data.SetLineColor(1)
    miny=0.
    maxe=0.
    for i in range(1,bkg.GetNbinsX()+1):
     e = bkg.GetBinError(i)
     if e> maxe: maxe=e
    miny=-maxe*1.4
    for i in range(1,data.GetNbinsX()+1): 
      x=data.GetBinContent(i) - data.GetBinError(i)
      if x < miny: miny=x 
    if miny<data.GetMinimum(): data.SetMinimum(miny)
    data.SetMaximum(data.GetMaximum()*1.6)
    data.Draw("E")

    col_sm = ROOT.kRed
    col_ps = ROOT.kBlue
    col_mm = ROOT.kGreen+2

    sm.SetLineWidth(3)
    sm.SetLineColor(col_sm)
    sm.SetMarkerSize(0)
    sm.SetFillStyle(0)

    ps.SetLineWidth(3)
    ps.SetLineColor(col_ps)
    ps.SetMarkerSize(0)
    ps.SetFillStyle(0)

    mm.SetLineWidth(3)
    mm.SetLineColor(col_mm)
    mm.SetMarkerSize(0)
    mm.SetFillStyle(0)

    hs.Add(ps)
    hs.Add(sm)
    if plot_mm:
        hs.Add(mm)

    hs.Draw("nostack hist same")

    bkg.SetFillColor(plot.CreateTransparentColor(12,0.4))
    bkg.SetLineColor(plot.CreateTransparentColor(12,0.4))
    bkg.SetMarkerSize(0)
    bkg.SetMarkerColor(plot.CreateTransparentColor(12,0.4))

    bkg.Draw("e2same")
    data.Draw("E same")

    relExtraDX = 2.0
    if len(extra_label) > 11:
        relExtraDX+=0.2*(len(extra_label)-11)
    DrawCMSLogoCustom(pads[0], 'CMS', extra_label, 11, 0.05, -0.07, 0.2, relExtraDX, '', 1.0)
    plot.DrawTitle(pads[0], '200 fb^{-1} (13 and 13.6 TeV)', 3)

    #Setup legend
    if plot_mm: legend = plot.PositionedLegend(0.8,0.25,1,0.02,0.00)
    else: legend = plot.PositionedLegend(0.8,0.13,1,0.,0.02)
    # use 2 columns
    legend.SetNColumns(4)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    legend.AddEntry(data,'Obs. #minus Bkg.',"lep")
    legend.AddEntry(bkg,'Bkg. unc.',"f")
    legend.AddEntry(sm,'#alpha^{H#tau#tau} = 0#lower[0.9]{^{#circ}}',"l")
    legend.AddEntry(ps,'#alpha^{H#tau#tau} = 90#lower[0.9]{^{#circ}}',"l")
    if plot_mm:
        legend.AddEntry(mm,'#alpha^{H#tau#tau} = 45#lower[0.9]{^{#circ}}',"l")
    legend.Draw("same")

    #latex.SetTextAlign(32)
    #latex.DrawLatex(0.92, 0.87, title)
 #
    #if title2 is not None:
    #  latex.SetTextAlign(32)
    #  #latex.DrawLatex(0.8, 0.80, title2)
    #  latex.DrawLatex(0.92, 0.80, title2)


    line = ROOT.TLine()
    line.SetLineWidth(1)
    line.SetLineStyle(2)
    line.SetLineColor(1)
    line.DrawLine(0.,0.,360.,0.)

    c1.SaveAs(plot_name+'.pdf')

propoganda_plot_phicp(sig_sm_best, sig_ps_best, sig_mm_best, bkg_best, data_best, f'{output_dir}/weighted_phiCP_10bin_categories', extra_label='Preliminary')
propoganda_plot_phicp(sig_sm_worst, sig_ps_worst, sig_mm_worst, bkg_worst, data_worst, f'{output_dir}/weighted_phiCP_other_categories', extra_label='Supplementary')

if args.extra_plots:
    # make a plot of run-2 only with best categories
    run2_cats = [cat for cat in best_cats if '13TeV' in cat ]
    bkg_run2, sig_sm_run2, sig_ps_run2, sig_mm_run2, data_run2 = CombineCats(run2_cats, histograms)
    propoganda_plot_phicp(sig_sm_run2, sig_ps_run2, sig_mm_run2, bkg_run2, data_run2, f'{output_dir}/weighted_phiCP_run2_only', extra_label='Preliminary')

    # make individual plots for each category including mm as well
    for cat in bin_set:
        if cat not in bins_map:
            continue

        bkg, sig_sm, sig_ps, sig_mm, data = CombineCats([cat], histograms)
        propoganda_plot_phicp(sig_sm, sig_ps, sig_mm, bkg, data, f'{output_dir}/weighted_phiCP_{cat}', extra_label='Preliminary', plot_mm=True)
