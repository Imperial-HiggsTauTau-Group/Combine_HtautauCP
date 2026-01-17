import CombineHarvester.CombineTools.ch as ch
import argparse
import ROOT
import math

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

phase_flip_categories = []

for year in [2016,2017,2018]:
    for i in [5,6,9,11]:
        phase_flip_categories.append(f'htt_tt_{year}_{i}_13TeV')
    for i in [3,4,5,6]:
        phase_flip_categories.append(f'htt_mt_{year}_{i}_13TeV')
        phase_flip_categories.append(f'htt_et_{year}_{i}_13TeV')
for i in [5,6,9,11]:
    phase_flip_categories.append(f'htt_tt_{i}_13p6TeV')
for i in [3,4,5]:
    phase_flip_categories.append(f'htt_mt_{i}_13p6TeV')
    phase_flip_categories.append(f'htt_et_{i}_13p6TeV')

parser = argparse.ArgumentParser(description='Post-fit plot script for Htautau CP analysis')
parser.add_argument('--fitresult', '-f', help= 'Path to a RooFitResult, only needed for postfit', default=None)
parser.add_argument('--workspace', '-w', help= 'The input workspace-containing file [REQUIRED]')
parser.add_argument('--output-folder', '-o', help= 'Output folder for datacards', default='cards_weighted_histograms')
parser.add_argument('--unblind', action='store_true', help='Unblind the data, if not set the data will be set the the Asimov dataset')
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
    if phase_shift:
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
for samp in range(samples): # note samp = 0 is nominal

    for cat in bin_set:
        # skip if not in our bins_map
        if cat not in bins_map:
            continue
 
        phase_flip = cat in phase_flip_categories
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
                    print(f"Setting bin {b} content to {data.GetBinContent(b)} +/- {data.GetBinError(b)}")
    
            # for the nominal case we also calculate the weights used to reweight the individual histograms 
            s_sb=0
            A_ave=0
            for b in range(1, bkg.GetNbinsX()+1):
                if (b-1) % nxbins ==0 and b+nxbins-1 <= bkg.GetNbinsX():
                    i_sm = sm_sig.Integral(b,b+nxbins-1)
                    i_ps = ps_sig.Integral(b,b+nxbins-1)
                    i_bkg = bkg.Integral(b,b+nxbins-1)
                    i_sig = (i_sm+i_ps)/2
                    s_sb = i_sig/(i_sig+i_bkg)
                    A_tot=0
                    for i in range(b, b+nxbins):
                        b_sm = sm_sig.GetBinContent(i)
                        b_ps = ps_sig.GetBinContent(i)
                        if b_sm+b_ps>1e-6:
                            A_tot += abs(b_sm-b_ps)/(b_sm+b_ps)

                    A_ave = A_tot/nxbins
                wt_mapping[cat][b] = s_sb*A_ave

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
        for i in range(1, bkg_variations[0].GetNbinsX()+1):
            err = abs(bkg_variations[0].GetBinContent(i)-bkg.GetBinContent(i))
            bkg.SetBinError(i, err*err + bkg.GetBinError(i))

    return bkg, sig_sm, sig_ps, sig_mm, data

#make a list of all the categories with 10 bins
best_cats = [cat for cat in bin_set if cat in bins_map and bins_map[cat]==10]

#make a list of the categories with or 8 bins i.e not = 10
worst_cats = [cat for cat in bin_set if cat in bins_map and bins_map[cat]!=10]

bkg_best, sig_sm_best, sig_ps_best, sig_mm_best, data_best = CombineCats(best_cats, histograms)
bkg_worst, sig_sm_worst, sig_ps_worst, sig_mm_worst, data_worst = CombineCats(worst_cats, histograms)

def Subtract(h1,h2):
  for i in range(1,h1.GetNbinsX()+1):
    diff = h1.GetBinContent(i) - h2.GetBinContent(i)
    h1.SetBinContent(i,diff)
  return h1

data_best = Subtract(data_best, bkg_best)
data_worst = Subtract(data_worst, bkg_worst)
bkg_best = Subtract(bkg_best, bkg_best)
bkg_worst = Subtract(bkg_worst, bkg_worst)

fout = ROOT.TFile('weighted_phiCP_histograms.root', 'RECREATE')
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

def propoganda_plot_phicp(sm,ps,mm, bkg,data,plot_name,extra_label='Preliminary'):

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(0.06)

    data.GetXaxis().SetTitleOffset(1.0)
    data.GetYaxis().SetTitleOffset(0.8)
    data.GetXaxis().SetTitleSize(0.05)
    data.GetYaxis().SetTitle('A#kern[0.1]{S}/(S+B) weighted events / bin')
    data.GetXaxis().SetTitle('#phi_{#it{CP}} (degrees)')
    data.GetYaxis().SetTitleSize(0.05)
    data.GetXaxis().SetNdivisions(506,False)

    c1 = ROOT.TCanvas()

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
    data.SetMaximum(data.GetMaximum()*1.5)
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

    hs.Draw("nostack hist same")

    bkg.SetFillColor(plot.CreateTransparentColor(12,0.4))
    bkg.SetLineColor(plot.CreateTransparentColor(12,0.4))
    bkg.SetMarkerSize(0)
    bkg.SetMarkerColor(plot.CreateTransparentColor(12,0.4))

    bkg.Draw("e2same")
    data.Draw("E same")

    plot.DrawCMSLogo(pads[0], 'CMS', extra_label, 11, 0.001, -0.07, 0.2, '', 1.0)
    plot.DrawTitle(pads[0], '200 fb^{-1} (13 and 13.6 TeV)', 3)

    #Setup legend
    legend = plot.PositionedLegend(0.25,0.25,1,0.02,0.00)
    legend.SetTextFont(42)
    legend.SetTextSize(0.05)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)

    legend.AddEntry(data,'Obs. #minus Bkg.',"lep")
    legend.AddEntry(bkg,'Bkg. unc.',"f")
    legend.AddEntry(sm,'#alpha^{H#tau#tau} = 0#lower[0.9]{^{#circ}}',"l")
    legend.AddEntry(ps,'#alpha^{H#tau#tau} = 90#lower[0.9]{^{#circ}}',"l")
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

propoganda_plot_phicp(sig_sm_best, sig_ps_best, sig_mm_best, bkg_best, data_best, 'weighted_phiCP_10bin_categories', extra_label='Preliminary')
propoganda_plot_phicp(sig_sm_worst, sig_ps_worst, sig_mm_worst, bkg_worst, data_worst, 'weighted_phiCP_other_categories', extra_label='Supplementary')
