import CombineHarvester.CombineTools.ch as ch
from argparse import ArgumentParser
import yaml
from CombineHarvester.Combine_HtautauCP.helpers import *

description = '''This script makes datacards with CombineHarvester for performing DM migration fits'''
parser = ArgumentParser(prog="harvestDatacardsDMMigrations",description=description,epilog="Success!")
parser.add_argument('-c', '--config', dest='config', type=str, default='configs/harvestDatacardsDMMigrations.yml', action='store', help="set config file")
args = parser.parse_args()

with open(args.config, 'r') as file:
   setup = yaml.safe_load(file)

output_folder = setup['output_folder']
input_folder = setup['input_folder']

cats = {
   (1, 'mt_DM0_tau_cp_mt_1_LT_50'),
   (2, 'mt_DM1_tau_cp_mt_1_LT_50'),
   (3, 'mt_DM2_tau_cp_mt_1_LT_50'),
   (4, 'mt_DM10_tau_cp_mt_1_LT_50')
}

sig_procs = ['ZTT_pi','ZTT_rho','ZTT_a11pr','ZTT_a13pr','ZTT_other']
bkg_procs = ['ZJ','ZL','TTT','TTJ','VVT','VVJ','W','QCD']

cb = ch.CombineHarvester()

cb.AddObservations(['*'], ['htt'], ['13p6TeV'], ['mt'], cats)
cb.AddProcesses(['*'], ['htt'], ['13p6TeV'], ['mt'], bkg_procs, cats, False)
cb.AddProcesses(['*'], ['htt'], ['13p6TeV'], ['mt'], sig_procs, cats, True)

filename = f'{input_folder}/dm_migrations_datacard.root'


# add systematics here
cb.cp().process(sig_procs).bin_id([1]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM0PNet', 'shape', ch.SystMap()(1.0))
cb.cp().process(sig_procs).bin_id([2]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM1PNet', 'shape', ch.SystMap()(1.0))
cb.cp().process(sig_procs).bin_id([3]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM2PNet', 'shape', ch.SystMap()(1.0))
cb.cp().process(sig_procs).bin_id([4]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM10PNet', 'shape', ch.SystMap()(1.0))

#Add lnN uncerts for the background split by DM
for i, dm in enumerate([0,1,2,10]):
   cb.cp().process(['ZL']).bin_id([i+1]).AddSyst(cb, f"LFake_uncert_DM{dm}PNet", "lnN", ch.SystMap()(1.2))
   cb.cp().process(['ZJ','VVJ','TTJ','WJ']).bin_id([i+1]).AddSyst(cb, f"JFake_uncert_DM{dm}PNet", "lnN", ch.SystMap()(1.2))
   cb.cp().process(['QCD']).bin_id([i+1]).AddSyst(cb, f"QCD_uncert_DM{dm}PNet", "lnN", ch.SystMap()(1.2))
   cb.cp().process(['TTT','VVT']).bin_id([i+1]).AddSyst(cb, f"T_uncert_DM{dm}PNet", "lnN", ch.SystMap()(1.2))

# Add rate params here

for i, dm in enumerate([0,1,2,10]):
   if dm == 0:     
      cb.cp().process(['ZTT_pi']).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_pi", "rateParam", ch.SystMap()(1.0))
      cb.cp().process(['ZTT_pi'],False).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_others", "rateParam", ch.SystMap()(1.0))
   elif dm == 1:
      cb.cp().process(['ZTT_rho']).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_rho", "rateParam", ch.SystMap()(1.0))
      cb.cp().process(['ZTT_rho'],False).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_others", "rateParam", ch.SystMap()(1.0))
   elif dm == 2:
      cb.cp().process(['ZTT_a11pr']).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_a11pr", "rateParam", ch.SystMap()(1.0))
      #cb.cp().process(['ZTT_a11pr'],False).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_others", "rateParam", ch.SystMap()(1.0)) ###
      cb.cp().process(['ZTT_rho']).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_rho", "rateParam", ch.SystMap()(1.0)) ###
      cb.cp().process(['ZTT_a11pr','ZTT_rho'],False).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_others", "rateParam", ch.SystMap()(1.0)) ###
   else:
      cb.cp().process(['ZTT_a13pr']).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_a13pr", "rateParam", ch.SystMap()(1.0))
      cb.cp().process(['ZTT_a13pr'],False).bin_id([i+1]).AddSyst(cb, f"rate_param_DM{dm}PNet_gen_others", "rateParam", ch.SystMap()(1.0))

# set ranges for rate params
for dm in [0,1,2,10]:
   if dm == 0:
      cb.GetParameter(f"rate_param_DM{dm}PNet_gen_pi").set_range(0.1,2.0)
   elif dm == 1:
      cb.GetParameter(f"rate_param_DM{dm}PNet_gen_rho").set_range(0.1,2.0)
   elif dm == 2:
      cb.GetParameter(f"rate_param_DM{dm}PNet_gen_a11pr").set_range(0.1,2.0)

      cb.GetParameter(f"rate_param_DM{dm}PNet_gen_rho").set_range(0.0,2.0) ###
   else:      
      cb.GetParameter(f"rate_param_DM{dm}PNet_gen_a13pr").set_range(0.1,2.0)
   cb.GetParameter(f"rate_param_DM{dm}PNet_gen_others").set_range(0.1,2.0)

print (">>>   file %s" % (filename))
cb.cp().era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")

rebin = ch.AutoRebin()
rebin.SetBinUncertFraction(0.2)
rebin.SetRebinMode(1)
rebin.SetPerformRebin(True)
rebin.SetVerbosity(1) 
rebin.Rebin(cb,cb)

# Zero negative bins
print(green(">>> Zeroing negative bins"))
cb.ForEachProc(NegativeBins)

print(green(">>> Zeroing negative yields"))
cb.ForEachProc(NegativeYields)

ch.SetStandardBinNames(cb)

# add bbb systs here
cb.SetAutoMCStats(cb, 0., 1, 1)

# Write datacards
print(green(">>> Writing datacards..."))
datacardtxt  = "%s/$TAG/$BIN.txt" % (output_folder)
datacardroot = "%s/$TAG/common/$BIN_input.root" % (output_folder)
writer = ch.CardWriter(datacardtxt,datacardroot)
writer.SetVerbosity(1)
writer.SetWildcardMasses([ ])
writer.WriteCards("cmb", cb)
# Cards per category
writer.WriteCards("DM0",   cb.cp().bin_id({1}))
writer.WriteCards("DM1",   cb.cp().bin_id({2}))
writer.WriteCards("DM2",   cb.cp().bin_id({3}))
writer.WriteCards("DM10",   cb.cp().bin_id({4}))
writer.WriteCards("oneprong",   cb.cp().bin_id({1,2,3}))