import CombineHarvester.CombineTools.ch as ch
from argparse import ArgumentParser
import yaml
from CombineHarvester.Combine_HtautauCP.helpers import *
from CombineHarvester.Combine_HtautauCP.systematics import AddSMRun3Systematics

# HI
description = '''This script makes datacards with CombineHarvester for performing tau ID SF measurments.'''
parser = ArgumentParser(prog="harvestDatacardsControlPlots",description=description,epilog="Success!")
parser.add_argument('-c', '--config', dest='config', type=str, default='configs/harvestDatacardsControlPlots.yml', action='store', help="set config file")
args = parser.parse_args()

with open(args.config, 'r') as file:
   setup = yaml.safe_load(file)

chans = setup['channels']
if chans == 'all': chans = ['tt','mt','et']
else: chans = chans.split(',')


output_folder = setup['output_folder']
input_folder = setup['input_folder']
var_name = setup['var_name']

# define background processes
bkg_procs_tt = ['ZTT','ZL','TTT','VVT','JetFakes','JetFakesSublead']
bkg_procs_lt = ['ZTT','ZL','TTT','VVT','JetFakes']
fake_procs = ['JetFakes','JetFakesSublead']

# define signal processes, which are the same for every channel
sig_procs = {}
sig_procs['ggH'] = ['ggH_sm_prod_sm_htt','ggH_ps_prod_sm_htt','ggH_mm_prod_sm_htt']
sig_procs['qqH'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt','WH_sm_htt','WH_ps_htt','WH_mm_htt','ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']

# define all MC procs
mc_procs = ['ZTT','ZL','TTT','VVT']
for p in sig_procs.values(): mc_procs+=p
# define categories which can depend on the channel
cats = {}
cats['tt'] = [
        (100, 'tt_cp_inclusive'),
]
        
cats['mt'] = [
        (100, 'mt_cp_inclusive_mTLt65'),
]
        
cats['et'] = [
        (100, 'et_cp_inclusive_mTLt65'),
        ]

# Create an empty CombineHarvester instance
cb = ch.CombineHarvester()

# Add processes and observations
for chn in chans:
    # Adding Data,Signal Processes and Background processes to the harvester instance
    cb.AddObservations(['*'], ['htt'], ['13p6TeV'], [chn], cats[chn])
    if chn == 'tt':
        cb.AddProcesses(['*'], ['htt'], ['13p6TeV'], [chn], bkg_procs_tt, cats[chn], False)
    else:
        cb.AddProcesses(['*'], ['htt'], ['13p6TeV'], [chn], bkg_procs_lt, cats[chn], False)
    cb.AddProcesses(['125'], ['htt'], ['13p6TeV'], [chn], sig_procs['ggH'], cats[chn], True)
    cb.AddProcesses(['125'], ['htt'], ['13p6TeV'], [chn], sig_procs['qqH'], cats[chn], True)

# Systematics are added here
cb = AddSMRun3Systematics(cb)

## Populating Observation, Process and Systematic entries in the harvester instance

for chn in chans:
    if chn == 'tt': filename = '%s/combined_earlyRun3/%s/datacard_%s_cp_inclusive_%s_full2223.root' % (input_folder, chn, var_name, chn)
    else: filename = '%s/combined_earlyRun3/%s/datacard_%s_mTLt65_cp_inclusive_%s_full2223.root' % (input_folder, chn, var_name, chn)
    print (">>>   file %s" % (filename))
    cb.cp().channel([chn]).backgrounds().process([]).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC") # add data shapes
    cb.cp().channel([chn]).backgrounds().era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
    for sig_proc in sig_procs.values(): 
        cb.cp().channel([chn]).process(sig_proc).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC")
 
# for QCD scale uncertainties we need to scale the yields to factor out any differences in XS
#TODO: will need updating onces datacards templates are renamed
for proc in ['ggH','qqH']:
    cb.cp().process(sig_procs[proc]).RenameSystematic(cb,"QCDscale_ren_signal",f"QCDscale_ren_{proc}_ACCEPT")
    cb.cp().process(sig_procs[proc]).RenameSystematic(cb,"QCDscale_fac_signal",f"QCDscale_fac_{proc}_ACCEPT")
    cb.cp().process(sig_procs[proc]).RenameSystematic(cb,"ps_isr_signal",f"ps_isr_{proc}")
    cb.cp().process(sig_procs[proc]).RenameSystematic(cb,"ps_fsr_signal",f"ps_fsr_{proc}")

cb.cp().syst_name(["QCDscale_ren_ggH_ACCEPT"]).ForEachSyst(lambda syst: (
      syst.set_value_u(syst.value_u() * 1/0.7605580771666764),
      syst.set_value_d(syst.value_d() * 1/1.2696408372342587)
))

cb.cp().syst_name(["QCDscale_fac_ggH_ACCEPT"]).ForEachSyst(lambda syst: (
      syst.set_value_u(syst.value_u() * 1/1.0605734162962437),
      syst.set_value_d(syst.value_d() * 1/0.9197774810421466)
))

cb.cp().syst_name(["QCDscale_ren_qqH_ACCEPT"]).ForEachSyst(lambda syst: (
      syst.set_value_u(syst.value_u() * 1/1.0025941737902164),
      syst.set_value_d(syst.value_d() * 1/0.9967738173425197)
))

cb.cp().syst_name(["QCDscale_fac_qqH_ACCEPT"]).ForEachSyst(lambda syst: (
      syst.set_value_u(syst.value_u() * 1/1.0057565776872635),
      syst.set_value_d(syst.value_d() * 1/0.9991435604512692)
))

# rename IP sig uncertainties to decorrelate electrons and muons
for tautype in ['prompt', 'tauDecay']:
    for eta in ['Lt1p0', '1p0to1p6', 'Gt1p6']:
        cb.cp().process(mc_procs).channel(['mt']).RenameSystematic(cb,f'CMS_HIG25012_eff_IPSigCut_{tautype}_eta_{eta}', f'CMS_HIG25012_eff_mu_IPSigCut_{tautype}_eta_{eta}')
        cb.cp().process(mc_procs).channel(['et']).RenameSystematic(cb,f'CMS_HIG25012_eff_IPSigCut_{tautype}_eta_{eta}', f'CMS_HIG25012_eff_e_IPSigCut_{tautype}_eta_{eta}')

# uncorrelate some FF uncertainties
cb.cp().process(['JetFakes']).channel(['tt']).RenameSystematic(cb,f'CMS_HIG25012_fake_t_sub_syst', f'CMS_HIG25012_fake_t_tt_sub_syst')
cb.cp().process(['JetFakes']).channel(['mt']).RenameSystematic(cb,f'CMS_HIG25012_fake_t_sub_syst', f'CMS_HIG25012_fake_t_mt_sub_syst')
cb.cp().process(['JetFakes']).channel(['et']).RenameSystematic(cb,f'CMS_HIG25012_fake_t_sub_syst', f'CMS_HIG25012_fake_t_et_sub_syst')

for chn in ['mt','et']:
    # uncorrelate all statistical components
    for ff_type in ['wj','qcd','mc_top']:
        for njets in [0,1,2]:
            cb.cp().process(['JetFakes']).channel([chn]).bin_id([1,2,4]).RenameSystematic(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm0_{njets}j", f"CMS_HIG25012_fake_t_{chn}_{ff_type}_stat_dm0_{njets}j")
            cb.cp().process(['JetFakes']).channel([chn]).bin_id([1,2,3]).RenameSystematic(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm1_{njets}j", f"CMS_HIG25012_fake_t_{chn}_{ff_type}_stat_dm1_{njets}j")
            cb.cp().process(['JetFakes']).channel([chn]).bin_id([1,2,6]).RenameSystematic(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm2_{njets}j", f"CMS_HIG25012_fake_t_{chn}_{ff_type}_stat_dm2_{njets}j")
            cb.cp().process(['JetFakes']).channel([chn]).bin_id([1,2,5]).RenameSystematic(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm10_{njets}j", f"CMS_HIG25012_fake_t_{chn}_{ff_type}_stat_dm10_{njets}j")
    # uncorrelate systematic component of et and mt for qcd only
    cb.cp().process(['JetFakes']).channel([chn]).RenameSystematic(cb, f"CMS_HIG25012_fake_t_qcd_syst", f"CMS_HIG25012_fake_t_{chn}_qcd_syst")

ch.SetStandardBinNames(cb)

cb.SetAutoMCStats(cb, 0., 1, 1)

# Implement fixes for negative bins and yields

# Zero negative bins
print(green(">>> Zeroing negative bins"))
cb.ForEachProc(NegativeBins)

print(green(">>> Zeroing negative yields"))
cb.ForEachProc(NegativeYields)

# Get nominal histograms for all processes (needed when setting systematics)
cb.ForEachProc(GetNominalHisto)
# raise RuntimeError("Stopping here for debugging purposes")
print(green(">>> Zeroing negative systematics"))
cb.ForEachSyst(DetectNegativeSyst)

# Write datacards
print(green(">>> Writing datacards..."))
datacardtxt  = "%s/$TAG/$BIN.txt" % (output_folder)
datacardroot = "%s/$TAG/common/$BIN_input.root" % (output_folder)
writer = ch.CardWriter(datacardtxt,datacardroot)
writer.SetVerbosity(1)
writer.SetWildcardMasses([ ])

writer.WriteCards("cmb", cb)
# Cards per channel
writer.WriteCards("tt", cb.cp().channel({"tt"}))
writer.WriteCards("mt", cb.cp().channel({"mt"}))
writer.WriteCards("et", cb.cp().channel({"et"}))



