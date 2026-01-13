# bins to phase shift:

# Run-2: 
    # tt: 5,6,9,11
    # lt: all (3,4,5,6)

# Run-3
    # tt: 5,6,9,11
    # lt: 3,4,5


# number of bins in each histogram

#10 bins: 
    #tt: 3,4,5,7
    #lt: 3,5

#4 bins:
    #tt 6,8,9,10,11

# 8 bins:
    #lt: 4,6 -> rebin these to 4 bins!

import CombineHarvester.CombineTools.ch as ch
import argparse
import ROOT

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

print(bin_set)

for b in bin_set:
    # skip if not in our bins_map
    if b not in bins_map:
        continue
    print(f'Processing bin: {b}')
    phase_flip = b in phase_flip_categories
    nxbins = bins_map[b]

    sel_bin = cb.cp().bin([b])

    bkg = sel_bin.cp().backgrounds().GetShape()

    # get SM and CP-odd signals
    # to get PS we scale alpha to 90 degrees
    par = cb.GetParameter('alpha')
    par.set_val(90.)
    ps_sig = sel_bin.cp().signals().process(['ggH_ps_prod_sm_htt','ggH_ps_htt','qqH_ps_htt','WH_ps_htt','ZH_ps_htt']).GetShape()
    # now get SM by setting alpha to 0 degrees
    par.set_val(0.)
    sm_sig = sel_bin.cp().signals().process(['ggH_sm_prod_sm_htt','ggH_sm_htt','qqH_sm_htt','WH_sm_htt','ZH_sm_htt']).GetShape()


    print(bkg.Integral(), sm_sig.Integral(), ps_sig.Integral())