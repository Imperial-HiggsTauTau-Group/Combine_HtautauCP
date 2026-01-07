import ROOT
import math
import argparse
ROOT.TH1.AddDirectory(False)
from ctypes import c_double
import os
import numpy as np
ROOT.gROOT.SetBatch(True)
 
# note that the merging of bins requires an even numnber of phi_CP bins, and this number must be set to the specific value used in the dictionary below otherwise the method will give incorrect results!

cp_bins = {
        "tt_mva_fake": 1,
        "tt_mva_tau": 1,
        "tt_mva_higgs": 1,
        "tt_higgs_pia1" : 4,
        "tt_higgs_a11pra1" : 4,
        "tt_higgs_rhoa1" : 10,
        "tt_higgs_pia11pr" : 4,
        "tt_higgs_pirho" : 10,
        "tt_higgs_a1a1" : 4,
        "tt_higgs_pipi" : 4,
        "tt_higgs_rhoa11pr" : 10,
        "tt_higgs_rhorho" : 10,
        "mt_mva_tau_mTLt65": 1,
        "mt_mva_fake_mTLt65": 1,
        "mt_higgs_mua1_mTLt65": 10,
        "mt_higgs_mua11pr_mTLt65": 8,
        "mt_higgs_mupi_mTLt65": 8,
        "mt_higgs_murho_mTLt65": 10,
        "et_mva_tau_mTLt65": 1,
        "et_mva_fake_mTLt65": 1,
        "et_higgs_ea1_mTLt65": 10,
        "et_higgs_ea11pr_mTLt65": 8,
        "et_higgs_epi_mTLt65": 8,
        "et_higgs_erho_mTLt65": 10,

        "tt_tau_pia1" : 4,
        "tt_tau_a11pra1" : 4,
        "tt_tau_rhoa1" : 10,
        "tt_tau_pia11pr" : 4,
        "tt_tau_pirho" : 10,
        "tt_tau_a1a1" : 4,
        "tt_tau_pipi" : 4,
        "tt_tau_rhoa11pr" : 10,
        "tt_tau_rhorho" : 10,

        "tt_fake_pia1" : 4,
        "tt_fake_a11pra1" : 4,
        "tt_fake_rhoa1" : 10,
        "tt_fake_pia11pr" : 4,
        "tt_fake_pirho" : 10,
        "tt_fake_a1a1" : 4,
        "tt_fake_pipi" : 4,
        "tt_fake_rhoa11pr" : 10,
        "tt_fake_rhorho" : 10,

        "mt_tau_mua1_mTLt65": 10,
        "mt_tau_mua11pr_mTLt65": 8,
        "mt_tau_mupi_mTLt65": 8,
        "mt_tau_murho_mTLt65": 10,
        "et_tau_ea1_mTLt65": 10,
        "et_tau_ea11pr_mTLt65": 8,
        "et_tau_epi_mTLt65": 8,
        "et_tau_erho_mTLt65": 10,

        "mt_fake_mua1_mTLt65": 10,
        "mt_fake_mua11pr_mTLt65": 8,
        "mt_fake_mupi_mTLt65": 8,
        "mt_fake_murho_mTLt65": 10,
        "et_fake_ea1_mTLt65": 10,
        "et_fake_ea11pr_mTLt65": 8,
        "et_fake_epi_mTLt65": 8,
        "et_fake_erho_mTLt65": 10, 
}


# background cats (to be added for chi2 tests)
# "tt_tau_pia1" : 4,
# "tt_tau_a11pra1" : 4,
# "tt_tau_rhoa1" : 10,
# "tt_tau_pia11pr" : 4,
# "tt_tau_pirho" : 10,
# "tt_tau_a1a1" : 4,
# "tt_tau_pipi" : 4,
# "tt_tau_rhoa11pr" : 10,
# "tt_tau_rhorho" : 10,
# "tt_fake_pia1" : 4,
# "tt_fake_a11pra1" : 4,
# "tt_fake_rhoa1" : 10,
# "tt_fake_pia11pr" : 4,
# "tt_fake_pirho" : 10,
# "tt_fake_a1a1" : 4,
# "tt_fake_pipi" : 4,
# "tt_fake_rhoa11pr" : 10,
# "tt_fake_rhorho" : 10,
# "mt_tau_mua1_mTLt65": 10,
# "mt_tau_mua11pr_mTLt65": 8,
# "mt_tau_mupi_mTLt65": 8,
# "mt_tau_murho_mTLt65": 10,
# "et_tau_ea1_mTLt65": 10,
# "et_tau_ea11pr_mTLt65": 8,
# "et_tau_epi_mTLt65": 8,
# "et_tau_erho_mTLt65": 10,
# "mt_fake_mua1_mTLt65": 10,
# "mt_fake_mua11pr_mTLt65": 8,
# "mt_fake_mupi_mTLt65": 8,
# "mt_fake_murho_mTLt65": 10,
# "et_fake_ea1_mTLt65": 10,
# "et_fake_ea11pr_mTLt65": 8,
# "et_fake_epi_mTLt65": 8,
# "et_fake_erho_mTLt65": 10,

test_results_sym = {}
test_results_flat = {}
test_results_asym = {}

ff_aiso_yields = {}
extra_hists = {}

def FlattenSystematic(nom_hist, syst_hist, nxbins):

  # flatten a systematic histogram

  ratio_hist = syst_hist.Clone()
  ratio_hist.Divide(nom_hist) # divide systematic by nominal to get ratio
  for i in range(1, ratio_hist.GetNbinsX() + 1):
    ratio_hist.SetBinError(i, 0)
  # flatten the ratio histogram
  ratio_flat, _, _ = MergeXBins(ratio_hist, nxbins)

  # get new systematic histogram by multiplying flattened ratio by nominal
  new_syst_hist = nom_hist.Clone()
  new_syst_hist.Multiply(ratio_flat)

  return new_syst_hist


def MergeXBins(hist, nxbins):
  histnew = hist.Clone()
  nbins = hist.GetNbinsX()
  chi2_total = 0.
  ndf_perbin = nxbins-1
  chi2_perbin=[]
  for i in range(1,nbins+1,nxbins):
    tot_err = c_double(0)
    tot = hist.IntegralAndError(i,i+nxbins-1,tot_err)
    tot_err = tot_err.value
    ave = tot / nxbins
    ave_err = tot_err / nxbins
    chi2 = 0.
    for j in range(i,i+nxbins):
      if hist.GetBinError(j) > 0: chi2 += (hist.GetBinContent(j)-ave)**2 / (hist.GetBinError(j)**2)
      histnew.SetBinContent(j,ave)
      histnew.SetBinError(j,ave_err)
    chi2_perbin.append(chi2)
  p_val_perbin = [ROOT.TMath.Prob(x, ndf_perbin) for x in chi2_perbin]
  chi2_total = sum(chi2_perbin)
  ndf_total = ndf_perbin * len(chi2_perbin)
  p_val_total = ROOT.TMath.Prob(chi2_total, ndf_total)

  return histnew, p_val_total, p_val_perbin, chi2_total, ndf_total

def MergeYBins(hist, nxbins):
  histnew = hist.Clone()
  nbins = hist.GetNbinsX()

  nybins = int(nbins/nxbins)

  for i in range(1,nxbins+1):
    tot = 0.
    tot_err = 0.
    for j in range(1,nybins+1):
      bin_num = i+(j-1)*nxbins
      tot += hist.GetBinContent(bin_num)
      tot_err += hist.GetBinError(bin_num)**2
    tot_err = math.sqrt(tot_err)
    # now set all j bins to this value
    for j in range(1,nybins+1):
      bin_num = i+(j-1)*nxbins
      histnew.SetBinContent(bin_num, tot)
      histnew.SetBinError(bin_num, tot_err)

  return histnew


def Symmetrise(hist,nxbins, verbose=False):
  if verbose:
      print('\n Processing', hist.GetName())
  histnew=hist.Clone()
  nbins = hist.GetNbinsX()
  if nbins % 2:
    print('N X bins in 2D histogram is not even so cannot symmetrise!')
    return
  nybins = int(nbins/nxbins)
  chi2_perbin = [0 for i in range(nybins)]
  ndf_perbin = int(nxbins / 2)
  ndf_total_alt=0
  for i in range(1,int(nxbins/2)+1):
    lo_bin = i
    hi_bin = nxbins-i+1
    chi2 = 0.
    for j in range(1,nybins+1):
      lo_bin_ = lo_bin+(j-1)*nxbins
      hi_bin_ = hi_bin+(j-1)*nxbins
      c1 = hist.GetBinContent(lo_bin_)
      c2 = hist.GetBinContent(hi_bin_)
      e1 = hist.GetBinError(lo_bin_)
      e2 = hist.GetBinError(hi_bin_)
      cnew = (c1+c2)/2
      enew = math.sqrt(e1**2 + e2**2)/2
      histnew.SetBinContent(lo_bin_,cnew)
      histnew.SetBinContent(hi_bin_,cnew)
      histnew.SetBinError(lo_bin_,enew)
      histnew.SetBinError(hi_bin_,enew)

      chi2 = 0.
      if verbose: 
          if e1<0: print('e1<0!', e1)
          if e2<0: print('e2<0!', e2)
      if (e1**2 + e2**2) > 0: 
          chi2 = (c1-c2)**2 / (e1**2 + e2**2)
          ndf_total_alt+=1
      chi2_perbin[j-1] += chi2
    p_val_perbin = [ ROOT.TMath.Prob(x, ndf_perbin) for x in chi2_perbin]
    chi2_total = sum(chi2_perbin)
    ndf_total = ndf_perbin * nybins
    p_val_total = ROOT.TMath.Prob(chi2_total, ndf_total)

  if verbose:
    print('New hist name:', histnew.GetName())

  return histnew, p_val_total, p_val_perbin, chi2_total, ndf_total

def ASymmetrise(hist,hsm,hps,nxbins):
  histnew=hist.Clone()
  hsub=hsm.Clone()
  hsub.Add(hps)
  hsub.Scale(0.5)
  for i in range(1,hsub.GetNbinsX()+1): 
    histnew.SetBinContent(i,histnew.GetBinContent(i)-hsub.GetBinContent(i))

  nbins = hist.GetNbinsX()
  if nbins % 2:
    print('N X bins in 2D histogram is not even so cannot symmetrise!')
    return
  nybins = int(nbins/nxbins)
  for i in range(1,int(nxbins/2)+1):
    lo_bin = i
    hi_bin = nxbins-i+1
    for j in range(1,nybins+1):
      lo_bin_ = lo_bin+(j-1)*nxbins
      hi_bin_ = hi_bin+(j-1)*nxbins

      mmi = hist.GetBinContent(lo_bin_)       
      mmj = hist.GetBinContent(hi_bin_)       
      smi = hsm.GetBinContent(lo_bin_) 
      smj = hsm.GetBinContent(hi_bin_)
      psi = hps.GetBinContent(lo_bin_)
      psj = hps.GetBinContent(hi_bin_)

      e_mmi = hist.GetBinError(lo_bin_)
      e_mmj = hist.GetBinError(hi_bin_)
      e_smi = hsm.GetBinError(lo_bin_)
      e_smj = hsm.GetBinError(hi_bin_)
      e_psi = hps.GetBinError(lo_bin_)
      e_psj = hps.GetBinError(hi_bin_) 

      c1_new = ( smj+psj-mmj + mmi)/2
      c2_new = ( smi+psi-mmi + mmj)/2
      e1_new = math.sqrt((e_smj+e_psj-e_mmj)**2 + e_mmi**2)/2 
      e2_new = math.sqrt((e_smi+e_psi-e_mmi)**2 + e_mmj**2)/2       

      histnew.SetBinContent(lo_bin_,c1_new)
      histnew.SetBinContent(hi_bin_,c2_new)
      histnew.SetBinError(lo_bin_,e1_new)
      histnew.SetBinError(hi_bin_,e2_new)

  return histnew

def GetDataAndEstFakes(dirname, infile):
    directory = infile.Get(dirname)
    # check if directory exists
    if not directory:
        print('Directory', dirname, 'not found in file', infile.GetName())
        return None, None
    data = directory.Get('data_obs')
    backgrounds = ['ZTT','ZL','TTT','VVT','JetFakesSublead']
    bkg = directory.Get(backgrounds[0])
    for b in backgrounds[1:]:
        bkg.Add(directory.Get(b))
    data.Add(bkg,-1)
    fake_est = directory.Get('JetFakes')
    if not data or not fake_est:
        print('Data or fake estimate not found in directory', dirname)
        return
    print(data, fake_est)
    return data, fake_est

def GetFFUncerts(dirname, infile):
    print('Getting FF uncertainties for:', dirname)
    dirname_aiso = dirname+'_aiso'

    extra_hists[dirname] = []

    data, fake_est = GetDataAndEstFakes(dirname_aiso, infile)

    ff_aiso_yields[dirname] = (data.Integral(), fake_est.Integral())

    if not data or not fake_est:
        return

    # now we estimate shape uncertainties for the BDT score and the aco-angle
    nxbins = cp_bins[dirname] if dirname in cp_bins else 1

    directory = infile.Get(dirname)
    fake_est_iso = directory.Get('JetFakes')
    rate = fake_est_iso.Integral()

    name_extra=""
    if 'higgs' in dirname:
        name_extra = '_higgs'
    elif 'tau' in dirname:
        name_extra = '_tau'
    elif 'fake' in dirname:
        name_extra = '_fake'

    for x in ['BDT', 'aco']:

        if x == 'BDT':
            # for the BDT shape uncertainties we rebin by nxbins so that the uncertainty does only depends on the BDT score and not the aco angle
            data_merged, _, _, _, _ = MergeXBins(data, nxbins)
            fake_est_merged, _, _,  _, _ = MergeXBins(fake_est, nxbins)
            name_extra=""
            if 'higgs' in dirname:
                name_extra = '_higgs'
            elif 'tau' in dirname:
                name_extra = '_tau'
            elif 'fake' in dirname:
                name_extra = '_fake'
        elif x == 'aco':
            # for the aco angle shape uncertainties we rebin by nybins so that the uncertainty does only depends on the aco angle and not the BDT score
            data_merged = MergeYBins(data, nxbins)
            fake_est_merged = MergeYBins(fake_est, nxbins)
            name_extra = '_'+dirname
            
        uncert_up = fake_est_iso.Clone()
        if uncert_up.GetNbinsX() != data_merged.GetNbinsX():
            print('Number of bins in fake estimate and fake estimate iso do not match, cannot calculate uncertainties!')
            return

        uncert_up.Multiply(data_merged)
        uncert_up.Divide(fake_est_merged)
        uncert_up.Scale(rate/uncert_up.Integral())
        uncert_down = fake_est_iso.Clone()
        uncert_down.Multiply(fake_est_merged)
        uncert_down.Divide(data_merged)
        uncert_down.Scale(rate/uncert_down.Integral())

        uncert_up.SetName(f'JetFakes_CMS_HIG25012_fake_t_syst_{x}shape{name_extra}Up')
        uncert_down.SetName(f'JetFakes_CMS_HIG25012_fake_t_syst_{x}shape{name_extra}Down')
    
        extra_hists[dirname] += [uncert_up, uncert_down]

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname, incData=False):
    directory = infile.Get(dirname)

    histos = []
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo, ROOT.TH1D) or isinstance(histo, ROOT.TH1F):
            histos.append(histo)
    # add extra_hists to histos
    if dirname in extra_hists:
        for histo in extra_hists[dirname]:
            if isinstance(histo, ROOT.TH1D) or isinstance(histo, ROOT.TH1F):
                histos.append(histo)


    print('Processing:', dirname) 
    for histo in histos:
        #print('Processing:', dirname, histo.GetName()) 
        if dirname in cp_bins: nxbins = cp_bins[dirname]
        else: nxbins = 1  
        # we write all the old histograms to the output file
        outfile.cd()
        if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
        ROOT.gDirectory.cd(dirname)
        histo.Write()    
        if nxbins== 1: continue    
        # if data we skip unless incData is True
        if 'data_obs' in histo.GetName() and not incData: continue    
        # if mm signal then we only anti-symmetrise
        elif 'H_mm' in histo.GetName() and 'htt125' in histo.GetName() and nxbins>1:
            hsm = directory.Get(histo.GetName().replace('H_mm_','H_sm_'))
            hps = directory.Get(histo.GetName().replace('H_mm_','H_ps_'))
            histo_asym = ASymmetrise(histo,hsm,hps,nxbins)
            histo_asym.SetName(histo.GetName()+'_sym')
            histo_asym.Write()
            continue    
        # if not mm signal then we always symmetrise
        else:
            histo_sym, p_val_total, p_val_perbin, chi2_total, ndf_total = Symmetrise(histo,nxbins, verbose=histo.GetName()=='data_obs')
            if dirname not in test_results_sym:
                test_results_sym[dirname] = {}
            test_results_sym[dirname][histo.GetName()] = (p_val_total, p_val_perbin[-1], chi2_total, ndf_total)
            histo_sym.SetName(histo.GetName()+'_sym')
            histo_sym.Write() 
            #print(histo.GetName(), histo_sym.GetName())
            # if not signal then we also flatten
            is_fake_cat = ('fake' in dirname)
            if ('htt125' not in histo.GetName() or 'Higgs_flat' in histo.GetName() or 'H_flat' in histo.GetName()) and not (is_fake_cat and 'data_obs' in histo.GetName()):
                print(dirname, histo.GetName(), ' - flattening')
                print(is_fake_cat, 'data_obs' in histo.GetName(), not (is_fake_cat and 'data_obs' in histo.GetName()))
                histo_flat, p_val_total, p_val_perbin, chi2_total, ndf_total = MergeXBins(histo,nxbins)
                if dirname not in test_results_flat:
                    test_results_flat[dirname] = {}
                test_results_flat[dirname][histo.GetName()] = (p_val_total, p_val_perbin[-1], chi2_total, ndf_total)
                histo_flat.SetName(histo.GetName()+'_flat')
                histo_flat.Write()

        ROOT.gDirectory.cd('/')


def getFlattenedSysts(infile,outfile,dirname,write_dirname, incData=False):
    directory = infile.Get(dirname)
    print("Flattening certain systematics:")
    histos = []
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo, ROOT.TH1D) or isinstance(histo, ROOT.TH1F):
            histos.append(histo)

    for histo in histos:
        is_syst = False # track whether this is a systematic variation
        if "_CMS" in histo.GetName():
            if "res_IP" in histo.GetName() or "res_SV" in histo.GetName():
                is_syst = False # expected to affect phiCP
            else:
                if "Up" in histo.GetName():
                    syst_name = '_CMS' + histo.GetName().split('_CMS')[1].split('Up')[0] + "Up"
                    is_syst=True
                elif "Down" in histo.GetName():
                    syst_name = '_CMS' + histo.GetName().split('_CMS')[1].split('Down')[0] + "Down"
                    is_syst=True


        if dirname in cp_bins: nxbins = cp_bins[dirname]
        else: nxbins = 1
        # we write all the old histograms to the output file
        outfile.cd()
        if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
        ROOT.gDirectory.cd(dirname)

        if not is_syst: # write original histo
          histo.Write()
        else:
          if nxbins== 1:
            histo.Write()
            continue
          # Get nominal histogram
          # print('Processing:', histo.GetName())
          nom_histo_name = histo.GetName().replace(syst_name,'')
          # print(f">> Getting nominal histo: {nom_histo_name}")
          nom_histo = directory.Get(nom_histo_name)
          if not nom_histo:
              raise ValueError(f"Could not find nominal histogram for systematic {histo.GetName()}")
          flattened_syst_histo = FlattenSystematic(nom_histo, histo, nxbins)
          flattened_syst_histo.SetName(histo.GetName())
          # print(f">> Writing flattened systematic histogram: {flattened_syst_histo.GetName()}"  )
          flattened_syst_histo.Write()



parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which we want to merge bins')
parser.add_argument('--test', '-t', action='store_true', help= 'Run statistical tests on the histograms to check compatability of smoothed and unsmoothed histograms')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','-mergeXbins.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename,"RECREATE") 
output_flatsyst_file = ROOT.TFile(filename.replace('.root','-mergeXbins_flatSyst.root'),"RECREATE")

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        incData = dirname.startswith("tt_tau_") or dirname.startswith("mt_tau_") or dirname.startswith("et_tau_") or dirname.startswith("tt_fake_") or dirname.startswith("mt_fake_") or dirname.startswith("et_fake_")
###        if dirname in cp_bins and dirname.startswith('tt_') and not dirname.startswith("tt_mva_higgs") and not dirname.startswith("tt_tau_") and not dirname.startswith("tt_fake_"): GetFFUncerts(dirname, original_file)
        getHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname, incData)
        # getFlattenedSysts(output_file,output_flatsyst_file,key.GetName(),dirname, incData)  # enable this to flatten systmatic variations

if 'tt_mva_higgs' in ff_aiso_yields and 'tt_mva_tau' in ff_aiso_yields and 'tt_mva_fake' in ff_aiso_yields:



    # now we estimate the normalization uncertainties for the FFs
    # get the correlated uncertainty from summing mva_higgs, mva_tau, and mva_fake
    data_total = 0.
    fake_total = 0.
    for x in ['tt_mva_higgs', 'tt_mva_tau', 'tt_mva_fake']:
        if x in ff_aiso_yields:
            data_total += ff_aiso_yields[x][0]
            fake_total += ff_aiso_yields[x][1]
    uncert_corr = data_total/fake_total

    # create a yml file to store the lnN uncertainties for the FFs
    with open('configs/ff_lnN_uncertainties.yml', 'w') as f:

        f.write('tt:')
        f.write(f'\n  correlated: {abs(1-uncert_corr)+1:.3f}\n')
    

        # now loop over 'tt' categories and get uncertainty per bin
        # in this case we scale the by the uncert_corr so that we don't double count the correlated uncertainty
        for dirname in ff_aiso_yields:
            if dirname.startswith('tt_'):
                data, fake_est = ff_aiso_yields[dirname]
                uncert = data / fake_est / (uncert_corr)
                uncert = abs(uncert - 1)+1
                f.write(f'  {dirname}: {uncert:.3f}\n')

def MakePlot(hists=[], filename='plot.pdf'):

    c = ROOT.TCanvas()
    hists[0].SetMinimum(0)
    hists[0].SetStats(0)
    hists[0].SetTitle('')
    hists[0].GetXaxis().SetTitle('#phi_{CP} (#circ)')
    hists[0].Draw()
    hists[1].SetLineColor(ROOT.kRed)
    hists[1].Draw('same')
    hists[2].SetLineColor(ROOT.kGreen-2)
    hists[2].Draw('same')

    sym_clone = hists[1].Clone()
    # zero errors for sym_clone
    for i in range(1,sym_clone.GetNbinsX()+1):
        sym_clone.SetBinError(i,0)
    flat_clone = hists[2].Clone()
    for i in range(1,flat_clone.GetNbinsX()+1):
        flat_clone.SetBinError(i,0)
    sym_KS_test = hists[0].KolmogorovTest(sym_clone,"WW")
    flat_KS_test = hists[0].KolmogorovTest(flat_clone,"WW")
    sym_chi_pvalue = test_results_sym[dirname][hists[0].GetName()][0]
    flat_chi_pvalue = test_results_flat[dirname][hists[0].GetName()][0]

    legend = ROOT.TLegend(0.1,0.1,0.9,0.4)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hists[0],'Original','l')
    legend.AddEntry(hists[1],f'Symmetrised, chi2 p-val = {sym_chi_pvalue:.2f}, KS p-val = {sym_KS_test:.2f}','l')
    legend.AddEntry(hists[2],f'Flattened, chi2 p-val = {flat_chi_pvalue:.2f}, KS p-val = {flat_KS_test:.2f}','l')
    legend.Draw()

    c.Print(filename)
    
def MakePlot(hists=[], filename='plot.pdf'):

    # Create canvas and pads
    c = ROOT.TCanvas("c","c",800,800)
    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1.0)  # top pad
    pad2 = ROOT.TPad("pad2","pad2",0,0.0,1,0.3)  # bottom pad
    pad1.SetBottomMargin(0)  # no x-labels on top plot
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad1.Draw()
    pad2.Draw()

    ## ---- TOP PAD ----
    pad1.cd()
    hists[0].SetMinimum(0)
    hists[0].SetStats(0)
    hists[0].SetTitle('')
    hists[0].GetXaxis().SetTitle('')
    hists[0].Draw("e")
    hists[1].SetLineColor(ROOT.kRed)
    hists[1].Draw("e SAME")
    hists[2].SetLineColor(ROOT.kGreen-2)
    hists[2].Draw("e SAME")

    # KS tests and chi2 p-values
    sym_clone = hists[1].Clone()
    for i in range(1,sym_clone.GetNbinsX()+1):
        sym_clone.SetBinError(i,0)
    flat_clone = hists[2].Clone()
    for i in range(1,flat_clone.GetNbinsX()+1):
        flat_clone.SetBinError(i,0)
    sym_KS_test = hists[0].KolmogorovTest(sym_clone,"WW")
    flat_KS_test = hists[0].KolmogorovTest(flat_clone,"WW")
    sym_chi_pvalue = test_results_sym[dirname][hists[0].GetName()][0]
    flat_chi_pvalue = test_results_flat[dirname][hists[0].GetName()][0]

    legend = ROOT.TLegend(0.1,0.1,0.9,0.4)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hists[0],'Original','el')
    legend.AddEntry(hists[1],f'Symmetrised, chi2 p-val = {sym_chi_pvalue:.2f}, KS p-val = {sym_KS_test:.2f}','el')
    legend.AddEntry(hists[2],f'Flattened, chi2 p-val = {flat_chi_pvalue:.2f}, KS p-val = {flat_KS_test:.2f}','el')
    legend.Draw()

    ## ---- BOTTOM PAD ----
    pad2.cd()

    ratio1 = hists[1].Clone()
    ratio1.Divide(hists[0])
    ratio1.SetLineColor(ROOT.kRed)
    ratio1.SetStats(0)
    ratio1.SetTitle("")
    ratio1.GetYaxis().SetTitle("Ratio")
    ratio1.GetYaxis().SetNdivisions(505)
    ratio1.GetYaxis().SetTitleSize(0.1)
    ratio1.GetYaxis().SetTitleOffset(0.5)
    ratio1.GetYaxis().SetLabelSize(0.08)

    ratio1.GetXaxis().SetTitle("#phi_{CP} (#circ)")
    ratio1.GetXaxis().SetTitleSize(0.1)
    ratio1.GetXaxis().SetLabelSize(0.08)

    ratio1.SetMinimum(0.85)
    ratio1.SetMaximum(1.15)
    ratio1.Draw("E")

    ratio2 = hists[2].Clone()
    ratio2.Divide(hists[0])
    ratio2.SetLineColor(ROOT.kGreen-2)
    ratio2.Draw("E SAME")

    #draw a dashed line at y=1
    line = ROOT.TLine(ratio1.GetXaxis().GetXmin(),1,ratio1.GetXaxis().GetXmax(),1)
    line.SetLineStyle(2)
    line.Draw("SAME")

    # Save
    c.Print(filename)

if args.test:

    for_bkg_cats = True # can change to true to make plots for the background only categories
    #if for_bkg_cats:
    #  n_flat = 3
    #  n_sym = 3
    #else:
    #  n_flat = 6
    #  n_sym = 6

    n_flat=6
    n_sym=7 

    if not os.path.exists('test_results'):
        os.makedirs('test_results')

    # make ROOT plots of the p-values
    # 1D plot just showing the p-values to show they are uniform
    # 2D plot showing processes on the x-axis, categories on the y-axis and p-values as the z-axis
    gr_sym_id = ROOT.TGraph()
    gr_sym_last_id = ROOT.TGraph()
    gr_flat_id = ROOT.TGraph()
    gr_flat_last_id = ROOT.TGraph()

    hist_sym_1d = ROOT.TH1D('hist_sym_1d','Symmetrisation tests',20,0,1)
    hist_flat_1d = ROOT.TH1D('hist_flat_1d','Flattening tests',20,0,1)
    hist_sym_last_1d = ROOT.TH1D('hist_sym_last_1d','Symmetrisation tests: last BDT bin',20,0,1)
    hist_flat_last_1d = ROOT.TH1D('hist_flat_last_1d','Flattening tests: last BDT bin',20,0,1)
 

    N_categories = len(test_results_sym)
    hist_sym_2d = ROOT.TH2D('hist_sym_2d','Symmetrisation tests',N_categories,0,N_categories,n_sym,0,n_sym)
    hist_flat_2d = ROOT.TH2D('hist_flat_2d','Flattening tests',N_categories,0,N_categories,n_flat,0,n_flat)
    hist_sym_last_2d = ROOT.TH2D('hist_sym_last_2d','Symmetrisation tests: last BDT bin',N_categories,0,N_categories,n_sym,0,n_sym)
    hist_flat_last_2d = ROOT.TH2D('hist_flat_last_2d','Flattening tests: last BDT bin',N_categories,0,N_categories,n_flat,0,n_flat)

    hist_sym_2d.GetZaxis().SetTitle('P-value')
    hist_flat_2d.GetZaxis().SetTitle('P-value')
    hist_sym_last_2d.GetZaxis().SetTitle('P-value (last BDT bin)')
    hist_flat_last_2d.GetZaxis().SetTitle('P-value (last BDT bin)')

    #if for_bkg_cats:
    #    hist_sym_2d.GetYaxis().SetBinLabel(1,'ZTT')
    #    hist_sym_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    #    hist_sym_2d.GetYaxis().SetBinLabel(3,'data_obs')
#
    #
    #    hist_flat_2d.GetYaxis().SetBinLabel(1,'ZTT')
    #    hist_flat_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    #    hist_flat_2d.GetYaxis().SetBinLabel(3,'data_obs')
#
    #
    #    hist_sym_last_2d.GetYaxis().SetBinLabel(1,'ZTT')
    #    hist_sym_last_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    #    hist_sym_last_2d.GetYaxis().SetBinLabel(3,'data_obs')
#
    #
    #    hist_flat_last_2d.GetYaxis().SetBinLabel(1,'ZTT')
    #    hist_flat_last_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    #    hist_flat_last_2d.GetYaxis().SetBinLabel(3,'data_obs')

    if True:

        hist_sym_2d.GetYaxis().SetBinLabel(1,'ZTT')
        hist_sym_2d.GetYaxis().SetBinLabel(2,'JetFakes')
        hist_sym_2d.GetYaxis().SetBinLabel(3,'qqH_sm_htt125')
        hist_sym_2d.GetYaxis().SetBinLabel(4,'ggH_sm_prod_sm_htt125')
        hist_sym_2d.GetYaxis().SetBinLabel(5,'ggH_ps_prod_ps_htt125')
        hist_sym_2d.GetYaxis().SetBinLabel(6,'ggH_flat_prod_mm_htt125')
        hist_sym_2d.GetYaxis().SetBinLabel(7,'data_obs')
    
        hist_flat_2d.GetYaxis().SetBinLabel(1,'ZTT')
        #hist_flat_2d.GetYaxis().SetBinLabel(2,'JetFakes')
        hist_flat_2d.GetYaxis().SetBinLabel(2,'qqH_flat_htt125')
        hist_flat_2d.GetYaxis().SetBinLabel(3,'ggH_flat_prod_sm_htt125')
        hist_flat_2d.GetYaxis().SetBinLabel(4,'ggH_flat_prod_ps_htt125')
        hist_flat_2d.GetYaxis().SetBinLabel(5,'ggH_flat_prod_mm_htt125')
        hist_flat_2d.GetYaxis().SetBinLabel(6,'data_obs')
    
        hist_sym_last_2d.GetYaxis().SetBinLabel(1,'ZTT')
        hist_sym_last_2d.GetYaxis().SetBinLabel(2,'JetFakes')
        hist_sym_last_2d.GetYaxis().SetBinLabel(3,'qqH_sm_htt125')
        hist_sym_last_2d.GetYaxis().SetBinLabel(4,'ggH_sm_prod_sm_htt125')
        hist_sym_last_2d.GetYaxis().SetBinLabel(5,'ggH_ps_prod_ps_htt125')
        hist_sym_last_2d.GetYaxis().SetBinLabel(6,'ggH_flat_prod_mm_htt125')
        hist_sym_last_2d.GetYaxis().SetBinLabel(7,'data_obs')
    
        hist_flat_last_2d.GetYaxis().SetBinLabel(1,'ZTT')
        #hist_flat_last_2d.GetYaxis().SetBinLabel(2,'JetFakes')
        hist_flat_last_2d.GetYaxis().SetBinLabel(2,'qqH_flat_htt125')
        hist_flat_last_2d.GetYaxis().SetBinLabel(3,'ggH_flat_prod_sm_htt125')
        hist_flat_last_2d.GetYaxis().SetBinLabel(4,'ggH_flat_prod_ps_htt125')
        hist_flat_last_2d.GetYaxis().SetBinLabel(5,'ggH_flat_prod_mm_htt125')
        hist_flat_last_2d.GetYaxis().SetBinLabel(6,'data_obs')

    hist_sym_2d.SetStats(0)
    hist_flat_2d.SetStats(0)
    hist_sym_last_2d.SetStats(0)
    hist_flat_last_2d.SetStats(0)

    hist_sym_1d.GetXaxis().SetTitle('P-value')
    hist_flat_1d.GetXaxis().SetTitle('P-value')
    hist_sym_last_1d.GetXaxis().SetTitle('P-value (last BDT bin)')
    hist_flat_last_1d.GetXaxis().SetTitle('P-value (last BDT bin)')

    hist_sym_1d.SetStats(0)
    hist_flat_1d.SetStats(0)
    hist_sym_last_1d.SetStats(0)
    hist_flat_last_1d.SetStats(0)
    # make sure 1D keeps track of errors
    hist_sym_1d.Sumw2()
    hist_flat_1d.Sumw2()
    hist_sym_last_1d.Sumw2()
    hist_flat_last_1d.Sumw2()

    #print('\n\nSymmetrisation test results:')
    i = 0
    for dirname, results in test_results_sym.items():

        if for_bkg_cats: 
            data = output_file.Get(dirname+'/data_obs')
            ztt = output_file.Get(dirname+'/ZTT')
            fakes = output_file.Get(dirname+'/JetFakes')
            data_sym = output_file.Get(dirname+'/data_obs_sym')
            ztt_sym = output_file.Get(dirname+'/ZTT_sym')
            fakes_sym = output_file.Get(dirname+'/JetFakes_sym')
            data_flat = output_file.Get(dirname+'/data_obs_flat')
            ztt_flat = output_file.Get(dirname+'/ZTT_flat')
            fakes_flat = output_file.Get(dirname+'/JetFakes_flat')

            if 'higgs' not in dirname and '_fake_' not in dirname: MakePlot([data,data_sym,data_flat], 'test_results/'+dirname+'_data.pdf')
            MakePlot([ztt,ztt_sym,ztt_flat], 'test_results/'+dirname+'_ZTT.pdf')
            MakePlot([fakes,fakes_sym,fakes_flat], 'test_results/'+dirname+'_fakes.pdf')

        i+=1

        # set x-labels
        hist_sym_2d.GetXaxis().SetBinLabel(i, dirname)
        hist_sym_last_2d.GetXaxis().SetBinLabel(i, dirname)
        for hist_name, (p_val_total, p_val_last, _, _) in results.items():
            #if not(hist_name in ['JetFakes', 'ZTT','Higgs_flat_htt125'] or hist_name.startswith('ggH') or hist_name.startswith('qqH_sm')) or 'unfiltered' in hist_name: continue
            #print(f'  Histogram: {hist_name}, P-value total: {p_val_total:.4f}, P-value last bin: {p_val_last:.4f}')
            j = None

            if hist_name == 'ZTT': j = 1
            if hist_name == 'JetFakes': j = 2
            if hist_name == 'qqH_sm_htt125': j = 3
            if hist_name == 'ggH_sm_prod_sm_htt125': j = 4
            if hist_name == 'ggH_ps_prod_ps_htt125': j = 5
            if hist_name == 'ggH_flat_prod_mm_htt125': j = 6
            if hist_name == 'data_obs': j = 7

            if j is None: continue

            hist_sym_2d.SetBinContent(i,j,p_val_total)
            hist_sym_last_2d.SetBinContent(i,j,p_val_last)

            if j: # only use statistically independent processes for 1D plot
              gr_sym_id.SetPoint(gr_sym_id.GetN(), p_val_total, 1)
              gr_sym_last_id.SetPoint(gr_sym_last_id.GetN(), p_val_last,1)
              hist_sym_1d.Fill(p_val_total)
              hist_sym_last_1d.Fill(p_val_last)

    #print('\n\nFlattening test results:')
    i=0
    count_removed_bins = 0
    for dirname, results in test_results_flat.items():

        # if dirname ends in pipi, a1a1, pia1, mupi_mTLt65, mua1_mTLt65, epi_mTLt65, ea1_mTLt65 then continue
        if True in [dirname.endswith(x) for x in ['pipi', 'a1a1', 'pia1', 'mupi_mTLt65', 'mua1_mTLt65', 'epi_mTLt65', 'ea1_mTLt65']]:
            count_removed_bins += 1
            continue

        i+=1
        #print(f'Directory: {dirname}')
        # set x-labels
        hist_flat_2d.GetXaxis().SetBinLabel(i, dirname)
        hist_flat_last_2d.GetXaxis().SetBinLabel(i, dirname)
        for hist_name, (p_val_total, p_val_last, _, _) in results.items():
            #if hist_name not in ['JetFakes', 'ZTT', 'qqH_flat_htt125', 'ggH_flat_prod_sm_htt125', 'ggH_flat_prod_ps_htt125', 'ggH_flat_prod_mm_htt125']: continue
            #print(f'  Histogram: {hist_name}, P-value total: {p_val_total:.4f}, P-value last bin: {p_val_last:.4f}')

            j = None

            if hist_name == 'ZTT': j = 1
            #if hist_name == 'JetFakes': j = 2
            if hist_name == 'qqH_flat_htt125': j = 2
            if hist_name == 'ggH_flat_prod_sm_htt125': j = 3
            if hist_name == 'ggH_flat_prod_ps_htt125': j = 4
            if hist_name == 'ggH_flat_prod_mm_htt125': j = 5
            if hist_name == 'data_obs': j = 6

            if j is None: continue

            hist_flat_2d.SetBinContent(i,j,p_val_total)
            hist_flat_last_2d.SetBinContent(i,j,p_val_last)

            if j: # only use statistically independent processes for 1D plot
              hist_flat_1d.Fill(p_val_total)
              hist_flat_last_1d.Fill(p_val_last)
              gr_flat_id.SetPoint(gr_flat_id.GetN(), p_val_total, 1)
              gr_flat_last_id.SetPoint(gr_flat_last_id.GetN(), p_val_last, 1)

    # change range based on number of removed bins
    hist_flat_2d.GetXaxis().SetRangeUser(0, hist_flat_2d.GetNbinsX()-count_removed_bins)
    hist_flat_last_2d.GetXaxis().SetRangeUser(0, hist_flat_last_2d.GetNbinsX()-count_removed_bins)

    # compute the p-value per row and per column for the 2D histograms - use the chi2 and ndf stored in the test_results dictionaries
    def computePValuePerRowAndColumnAndTotal(hist_2d, test_results, add_to_hist=True):
        nrows = hist_2d.GetNbinsY()
        ncols = hist_2d.GetNbinsX()

        chi2_overall = 0
        ndf_overall = 0

        # overall
        for i in range(1, nrows+1):
            for j in range(1, ncols+1):
                dirname = hist_2d.GetXaxis().GetBinLabel(j)
                hist_name = hist_2d.GetYaxis().GetBinLabel(i)
                if dirname in test_results and hist_name in test_results[dirname]:
                    chi2, ndf = test_results[dirname][hist_name][2], test_results[dirname][hist_name][3]
                    chi2_overall += chi2
                    ndf_overall += ndf
        if ndf_overall > 0:
            p_val_overall = ROOT.TMath.Prob(chi2_overall, ndf_overall)
            print(f'Overall: chi2 = {chi2_overall}, ndf = {ndf_overall}, p-value = {p_val_overall:.4f}')

        # if add_to_hist is True then set the by row and by collumn p-values in the hist_2d
        if add_to_hist:
            print('!test!')
            #hist_2d.GetXaxis().SetRangeUser(0, ncols+1)
            #hist_2d.GetYaxis().SetRangeUser(0, nrows+1)
            hist_2d_new = ROOT.TH2D(hist_2d.GetName(), hist_2d.GetTitle(), ncols+1, 0, ncols+1, nrows+1, 0, nrows+1)
            for i in range(1, nrows+1):
                for j in range(1, ncols+1):
                    hist_2d_new.SetBinContent(j, i, hist_2d.GetBinContent(j, i))
                    hist_2d_new.GetXaxis().SetBinLabel(j, hist_2d.GetXaxis().GetBinLabel(j))
                    hist_2d_new.GetYaxis().SetBinLabel(i, hist_2d.GetYaxis().GetBinLabel(i))
            hist_2d = hist_2d_new.Clone()
            hist_2d.SetStats(0)
            hist_2d.GetXaxis().SetBinLabel(ncols+1, 'p-value per row')
            hist_2d.GetYaxis().SetBinLabel(nrows+1, 'p-value per column')
            

        # per row
        for i in range(1, nrows+1):
            chi2_row = 0
            ndf_row = 0
            for j in range(1, ncols+1):
                dirname = hist_2d.GetXaxis().GetBinLabel(j)
                hist_name = hist_2d.GetYaxis().GetBinLabel(i)
                if dirname in test_results and hist_name in test_results[dirname]:
                    chi2, ndf = test_results[dirname][hist_name][2], test_results[dirname][hist_name][3]
                    chi2_row += chi2
                    ndf_row += ndf
            if ndf_row > 0:
                p_val_row = ROOT.TMath.Prob(chi2_row, ndf_row)
                if add_to_hist:
                    hist_2d.SetBinContent(ncols+1, i, p_val_row)
                print(f'Row {hist_2d.GetYaxis().GetBinLabel(i)}: chi2 = {chi2_row}, ndf = {ndf_row}, p-value = {p_val_row:.4f}')

        # per column
        for j in range(1, ncols+1):
            chi2_col = 0
            ndf_col = 0
            for i in range(1, nrows+1):
                dirname = hist_2d.GetXaxis().GetBinLabel(j)
                hist_name = hist_2d.GetYaxis().GetBinLabel(i)
                if dirname in test_results and hist_name in test_results[dirname]:
                    chi2, ndf = test_results[dirname][hist_name][2], test_results[dirname][hist_name][3]
                    chi2_col += chi2
                    ndf_col += ndf
            if ndf_col > 0:
                p_val_col = ROOT.TMath.Prob(chi2_col, ndf_col)
                if add_to_hist:
                    hist_2d.SetBinContent(j, nrows+1, p_val_col)
                print(f'Column {hist_2d.GetXaxis().GetBinLabel(j)}: chi2 = {chi2_col}, ndf = {ndf_col}, p-value = {p_val_col:.4f}')
    
        return p_val_overall, hist_2d

    print('\n\nSymmetrisation 2D test results:')
    pval_overall_sym, hist_sym_2d = computePValuePerRowAndColumnAndTotal(hist_sym_2d, test_results_sym)
    print('\n\nFlattening 2D test results:')
    p_val_overall_flat, hist_flat_2d = computePValuePerRowAndColumnAndTotal(hist_flat_2d, test_results_flat)

    # add this p-value in the title of the 2D histograms
    hist_sym_2d.SetTitle(f'Symmetrisation tests (overall p-value = {pval_overall_sym:.2f})')
    hist_flat_2d.SetTitle(f'Flattening tests (overall p-value = {p_val_overall_flat:.2f})')

    fout = ROOT.TFile('test_results/test_results.root', 'RECREATE')
    fout.cd()
    hist_sym_2d.Write()
    hist_flat_2d.Write()
    hist_sym_last_2d.Write()
    hist_flat_last_2d.Write()
    gr_sym_id.Write('gr_sym_id')
    gr_sym_last_id.Write('gr_sym_last_id')
    gr_flat_id.Write('gr_flat_id')
    gr_flat_last_id.Write('gr_flat_last_id')
    hist_sym_1d.Write()
    hist_flat_1d.Write()
    hist_sym_last_1d.Write()
    hist_flat_last_1d.Write()

    c_2d = ROOT.TCanvas('', '', 1200, 600)
    # decrease number of sf shown in plot
    ROOT.gStyle.SetPaintTextFormat(".2g")
    ROOT.gStyle.SetPalette(ROOT.kBlackBody)

    #increase L and R margins
    c_2d.SetLeftMargin(0.15)
    c_2d.SetRightMargin(0.13)
    hist_sym_2d.GetZaxis().SetRangeUser(0,1)
    hist_sym_2d.Draw('COLZ TEXT')
    # draw lines to seperate the last row and column
    line = ROOT.TLine(hist_sym_2d.GetXaxis().GetXmin(), hist_sym_2d.GetYaxis().GetBinLowEdge(hist_sym_2d.GetNbinsY()), hist_sym_2d.GetXaxis().GetXmax(), hist_sym_2d.GetYaxis().GetBinLowEdge(hist_sym_2d.GetNbinsY()))
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(1)
    line.Draw('SAME')
    line2 = ROOT.TLine(hist_sym_2d.GetXaxis().GetBinLowEdge(hist_sym_2d.GetNbinsX()), hist_sym_2d.GetYaxis().GetXmin(), hist_sym_2d.GetXaxis().GetBinLowEdge(hist_sym_2d.GetNbinsX()), hist_sym_2d.GetYaxis().GetXmax())
    line2.SetLineColor(ROOT.kBlack)
    line2.SetLineWidth(1)
    line2.Draw('SAME')
    c_2d.Print('test_results/chi2tests_symmetrisation_2d.pdf')

    hist_sym_last_2d.GetZaxis().SetRangeUser(0,1)
    hist_sym_last_2d.Draw('COLZ TEXT')
    c_2d.Print('test_results/chi2tests_symmetrisation_last_2d.pdf')

    hist_flat_2d.GetZaxis().SetRangeUser(0,1)
    hist_flat_2d.Draw('COLZ TEXT')
    c_2d.Print('test_results/chi2tests_flattening_2d.pdf')

    hist_flat_last_2d.GetZaxis().SetRangeUser(0,1)
    hist_flat_last_2d.Draw('COLZ TEXT')
    c_2d.Print('test_results/chi2tests_flattening_last_2d.pdf')

    c_1d = ROOT.TCanvas('', '', 800, 600)
    hist_sym_1d.Draw()
    c_1d.Print('test_results/chi2tests_symmetrisation_1d.pdf')

    hist_flat_1d.Draw()
    c_1d.Print('test_results/chi2tests_flattening_1d.pdf')

    hist_sym_last_1d.Draw()
    c_1d.Print('test_results/chi2tests_symmetrisation_last_1d.pdf')

    hist_flat_last_1d.Draw()
    c_1d.Print('test_results/chi2tests_flattening_last_1d.pdf')
