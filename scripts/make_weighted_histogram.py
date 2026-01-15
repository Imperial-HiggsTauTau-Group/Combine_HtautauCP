import CombineHarvester.CombineTools.ch as ch
import argparse
import ROOT
import math

ROOT.gROOT.SetBatch(True)

samples=500

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
if args.fitresult:
    f_fit = ROOT.TFile(args.fitresult.split(':')[0])
    res = f_fit.Get(args.fitresult.split(':')[1])

    cb.UpdateParameters(res)

bin_set = cb.bin_set()

def apply_phase_shift(hist, nxbins, bin):
    half = nxbins // 2

    histout = hist.Clone()

    for b in range(1, hist.GetNbinsX()+1):

        new_bin=b+half - (b+half - nxbins*((b-1) // nxbins) > nxbins)*nxbins

        print(b, new_bin)

        histout.SetBinContent(new_bin, hist.GetBinContent(b))
        histout.SetBinError(new_bin, hist.GetBinError(b))


    return histout

def remap_and_set_shape(x, wt_mapping, nxbins, phase_shift, proto1, proto2, proto3, obs=False):

    old_h = x.shape()

    use_proto1 = nxbins == 10
    use_proto2 = nxbins == 4
    use_proto3 = nxbins == 8

    if use_proto1:
        new_h = proto1.Clone()
    elif use_proto2:
        new_h = proto2.Clone()
    elif use_proto3:
        new_h = proto3.Clone()

    new_h.SetDirectory(0)
    new_h.Reset("ICES")

    nbins_old = old_h.GetNbinsX()
    half = nxbins // 2

    for b in range(1, nbins_old + 1):
        new_bin = (b - 1) % nxbins + 1

        if phase_shift:
            new_bin += half
            if new_bin > nxbins:
                new_bin -= nxbins

        wt = wt_mapping[b]

        new_h.SetBinContent(
            new_bin,
            new_h.GetBinContent(new_bin) + wt * old_h.GetBinContent(b)
        )

        new_h.SetBinError(
            new_bin,
            math.sqrt(pow(new_h.GetBinError(new_bin),2) + pow(wt*old_h.GetBinError(b),2))
        )

    # if has 8 bins, rebin to 4
    if use_proto3:
        new_h.Rebin(2)

    new_h.Scale(x.rate())
    x.set_shape(new_h, True)


for cat in bin_set:
    # skip if not in our bins_map
    if cat not in bins_map:
        continue
    print(f'Processing category: {cat}')
    phase_flip = cat in phase_flip_categories
    nxbins = bins_map[cat]

    sel_bin = cb.cp().bin([cat])

    bkg = sel_bin.cp().backgrounds().GetShape()

    print(f'Number of bins in selected for {cat} category: {bkg.GetNbinsX()}')

    # get SM and CP-odd signals
    # to get PS we scale alpha to 90 degrees
    par = cb.GetParameter('alpha')
    par.set_val(90.)
    ps_sig = sel_bin.cp().signals().process(['ggH_ps_prod_sm_htt','ggH_ps_htt','qqH_ps_htt','WH_ps_htt','ZH_ps_htt']).GetShape()
    # now get SM by setting alpha to 0 degrees
    par.set_val(0.)
    sm_sig = sel_bin.cp().signals().process(['ggH_sm_prod_sm_htt','ggH_sm_htt','qqH_sm_htt','WH_sm_htt','ZH_sm_htt']).GetShape()

    if args.unblind: 
        data = sel_bin.cp().GetObservedShape()
    else: 
        data = sel_bin.cp().GetShape()
        for b in range(1, data.GetNbinsX()+1):
            data.SetBinError(b, math.sqrt(data.GetBinContent(b)))

    # apply phase shift if needed
    if phase_flip:
        sm_sig_orig = sm_sig.Clone()
        print('Applying phase flip for category ', cat)
        sm_sig = apply_phase_shift(sm_sig, nxbins, cat)
        ps_sig = apply_phase_shift(ps_sig, nxbins, cat)
        data = apply_phase_shift(data, nxbins, cat)
        bkg = apply_phase_shift(bkg, nxbins, cat)

        #for testing make a plot of the sm_sig after phase shift
        c1 = ROOT.TCanvas("c1", "c1", 800,600)
        sm_sig.Draw("hist")
        sm_sig_orig.SetLineColor(ROOT.kRed)
        sm_sig_orig.Draw("hist same")
        c1.Print(f"test_plots/sm_sig_phaseshift_{cat}.pdf")

    wt_mapping = {}
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
        wt_mapping[b] = s_sb*A_ave

##make a list of all the categories with 10 bins
#best_cats = [cat for cat in bin_set if cat in bins_map and bins_map[cat]==10]

##make a list of the categories with or 8 bins i.e not = 10
#worst_cats = [cat for cat in bin_set if cat in bins_map and bins_map[cat]!=10]

