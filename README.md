## Setup Instructions:


### Combine Harvester

This package requires HiggsAnalysis/CombinedLimit to be in your local CMSSW area. We follow the release recommendations of the combine developers. The CombineHarvester framework is compatible with the CMSSW 14_1_X and 11_3_X series releases. The default branch, main, is for developments in the 14_1_X releases, and the current recommended tag is v3.0.0. The v2.1.0 tag should be used in CMSSW_11_3_X.

A new full release area can be set up and compiled in the following steps:

```
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
cd CombineHarvester
git checkout v3.0.0
scram b -j8
```

### CombineHarvester Repository for the Htautau CP in decay measurement

```
cd $CMSSW_BASE/src/CombineHarvester
git clone git@github.com:Ksavva1021/Combine_HtautauCP.git Combine_HtautauCP
scram b -j8
```

### The magic command

```
ulimit -s unlimited
```

### Pre-process datacards

The current mt datacard format requires the BDT bin categories to be combined into one histogram, this is done by the following script:

```
python3 scripts/process_mt_cards.py -f cpdatacards/mt_2022_2023.root
```


Then apply flattening/symmetrisation to datacards for both mt and tt channels using: 

```
python3 scripts/convertCards.py -f  cpdatacards/mt_2022_2023_merged.root
python3 scripts/convertCards.py -f  cpdatacards/added_histo_Run2Bins.root

```

### Produce txt datacards

Modify the options configs/harvestDatacards.yml as needed and then run
 
```
python3 scripts/harvestDatacards.py -c configs/harvestDatacards.yml 
```

The systematics can be modified in python/systematics.py


### Creating the workspaces

```
combineTool.py -m 125 -M T2W -P CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays -i outputs/cmb -o ws.root --parallel 8
```

### Run maximum likelihood fits

1D fit for alpha:

```
combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90 --points 21 --redefineSignalPOIs alpha  -d outputs/cmb/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1
```

TODO: add instructions for running points as batch jobs

### make plot of alpha scan

```
python3 scripts/plot1DScan.py --main=outputs/cmb/higgsCombine.alpha.MultiDimFit.mH125.root --POI=alpha --output=alpha_cmb --no-numbers --no-box --x-min=-90 --x-max=90 --y-max=8
```

### making prefit plots

Produce ROOT file containing all prefit histograms:
```
python3 python/PostFitShapesCombEras.py -w outputs/cmb/ws.root -d outputs/cmb/combined.txt.cmb 
```

If we want to show pseudoscalar on the same plot then need to run this again and freeze alpha to 90:
```
python3 python/PostFitShapesCombEras.py -w outputs/cmb/ws.root -d outputs/cmb/combined.txt.cmb --output shapes_output_ps.root --freeze alpha=90
```

Make plots from this scripts setting the -b option to the name of the bin you want to plot:

```
python3 scripts/postfitPlot.py -b htt_tt_3_13p6TeV
```

### running approximate impacts

For some reason the Impacts with --approx option does not run the initial fit and just complains it doesn't exist, so we run this ourselves and rename it using: 

```
combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d outputs/cmb/ws.root -t -1 -n .alpha.Impacts --alignEdges 1 --saveFitResult

mv multidimfit.alpha.Impacts.root multidimfit_approxFit_.alpha.Impacts.root
```

Then we make the json with the approximate impacts using:

```
combineTool.py -m 125 -M Impacts --setParameters muV=1,alpha=0,muggH=1,mutautau=1 -d outputs/cmb/ws.root -t -1 -n .alpha.Impacts  --approx=hesse --output impacts.json
```

And then plot using:

```
plotImpacts.py -i impacts.json -o impacts
```

### Running full impacts

Run initial fit:

```
combineTool.py -M Impacts -m 125 --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d outputs/cmb/ws.root -t -1 -n .alpha.Impacts --doInitialFit --robustFit=1 --cminDefaultMinimizerStrategy=0
```

Then run fits for each parameter as batch jobs, in the example below running 10 parameters per job:

```
combineTool.py -M Impacts -m 125 --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d outputs/cmb/ws.root -t -1 -n .alpha.Impacts --doFits --robustFit=1 --cminDefaultMinimizerStrategy=0 --job-mode condor --task-name condor-impacts --sub-opts='+MaxRuntime=10800' --merge 10
```

Collect the outputs and make a json file:

```
combineTool.py -M Impacts -m 125 --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d outputs/cmb/ws.root -t -1 -n .alpha.Impacts --exclude lumi_13p6TeV -o impacts.json
```

Note you may find that some of the fits did not produce a propper output and they get skipped, e.g in my case the lumi parameter. 
You may want to investigate what actually went wrong in the fit and re-run it for this specific parameter.

You can run the fits for a specific parameter using e.g:

```
combine -M MultiDimFit -n _paramFit_.alpha.Impacts_lumi_13p6TeV --algo impact --redefineSignalPOIs alpha -P lumi_13p6TeV --floatOtherPOIs 1 --saveInactivePOI 1 --setParameterRanges alpha=-90,90 -t -1 --robustFit=1 -m 125 -d outputs/cmb/ws.root --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --cminDefaultMinimizerStrategy=0
```
Then debug and fix as needed

Finally, make the plot using:

```
plotImpacts.py -i impacts.json -o impacts
```

It can also be informative to look at impacts for the signal strength, in this example we do t
his for a common signal strength (scale both ggH and VBF with same rate param):

```
combineTool.py -M Impacts -m 125 --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90:mutautau=0,10  --redefineSignalPOIs mutautau --freezeParameters=muggH,muV  -d outputs/cmb/ws.root -t -1 -n .mutautau.Impacts --doInitialFit --robustFit=1 --cminDefaultMinimizerStrategy=0

combineTool.py -M Impacts -m 125 --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90:mutautau=0,10  --redefineSignalPOIs mutautau --freezeParameters=muggH,muV  -d outputs/cmb/ws.root -t -1 -n .mutautau.Impacts --doFits --robustFit=1 --cminDefaultMinimizerStrategy=0 --exclude muV,muggH --job-mode condor --task-name condor-impacts --sub-opts='+MaxRuntime=10800' --merge 10
```
Wait for jobs to finish then collect and plot:
```
combineTool.py -M Impacts -m 125 --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90:mutautau=0,10  --redefineSignalPOIs mutautau --freezeParameters=muggH,muV  -d outputs/cmb/ws.root -t -1 -n .mutautau.Impacts --exclude muV,muggH -o impacts_mu.json

plotImpacts.py -i impacts_mu.json -o impacts_mu
```

### Running combination with Run-2

checkout the Run-2 datacards

```
git clone ssh://git@gitlab.cern.ch:7999/cms-analysis/hig/HIG-20-006/datacards.git Run2_datacards
```

make a directory to collect all the datacards for Run-2 and Run-3

```
mkdir run2run3_comb
```

copy Run-2 and Run-3 datacards to this folder

```
cp -r Run2_datacards/for_hig-25-012_combination/* run2run3_comb/. 
cp outputs/cmb/htt_*.txt run2run3_comb/.
cp outputs/cmb/common/htt_*.root run2run3_comb/common/.
```

don't forget this!
```
ulimit -s unlimited
```

make workspace

```
combineTool.py -m 125 -M T2W -P CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays -i run2run3_comb/ -o ws.root --parallel 8 
```

run the fits. Note these take a while so they are submitted as batch jobs
The fallback options are also needed as especially far from the best fit values the fits can struggle to converge

```
combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90 --points 21 --redefineSignalPOIs alpha  -d run2run3_comb/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10  --job-mode condor --task-name condor-run2run3-scan --sub-opts='+MaxRuntime=10800' --split-points 1
```

then just hadd them, and run the 1D plotting code as normal

### Make weighted phi-CP plot

Note current instructions are for expected plot. Replacing data with azimov. Need to update options after unblinding

First run a ML fit:

```
combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90  --redefineSignalPOIs alpha  -d latest_combined_datacards_Jan07/ws.root -t -1 -n .alpha.BestFit --algo none --saveFitResult --there --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1
```

This will produce 'multidimfit.alpha.BestFit.root'

The weighted plot is then produced using:

```
python3 -u scripts/make_weighted_histogram.py -w latest_combined_datacards_Jan07/ws.root -f latest_combined_datacards_Jan07/multidimfit.alpha.BestFit.root:fit_mdf
```

This replaces data with Asimov by default. To show the real data add the option --unblind
Note this may take a while to run as sampling is used to get the systematic uncertainty band

### Scan of mu vs alpha

Use a directory with all the required datacards inside it, and run T2W with the float mu physics model option

```
combineTool.py -m 125 -M T2W -P CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays -i run2run_muVsalpha/ --PO float_mutautau -o ws_muVsalpha.root --parallel 8
```

Run as batch jobs (navigate to the directory if you want the logs etc to end up there)
```
combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --redefineSignalPOIs alpha,mutautau --points 2000 -d ws_muVsalpha.root --algo grid -t -1 --there -n .muVsalpha --alignEdges 1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10 --job-mode condor --task-name condor-run2run3-muVsalpha --sub-opts='+MaxRuntime=10799' --split-points 5
```

then combine:

```
hadd -v 1 -f higgsCombine.muVsalpha.MultiDimFit.mH125 higgsCombine.muVsalpha.POINTS*.MultiDimFit.mH125.root
```

and make the plot:

```
python3 scripts/plot_2D_scans.py --file run2run_muVsalpha/higgsCombine.muVsalpha.MultiDimFit.mH125  --mutautau
```

### Scan of kappa tilde vs tilde


Use a directory with all the required datacards inside it, and run T2W with the kappas physics model option

```
combineTool.py -m 125 -M T2W -P CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays -i run2run_muVsalpha/ --PO float_mutautau -o ws_muVsalpha.root --parallel 8
```

navigate to dir
```
combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --redefineSignalPOIs alpha,mutautau --points 2000 -d ws_muVsalpha.root --algo grid -t -1 --there -n .muVsalpha --alignEdges 1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10 --job-mode condor --task-name condor-run2run3-muVsalpha --sub-opts='+MaxRuntime=10799' --split-points 5
```


```
hadd -v 1 -f higgsCombine.muVsalpha.MultiDimFit.mH125 higgsCombine.muVsalpha.POINTS*.MultiDimFit.mH125.root
```

```
python3 scripts/plot_2D_scans.py --file run2run3_2Dkappas/higgsCombine.kappas.MultiDimFit.mH125.root --kappas
```
