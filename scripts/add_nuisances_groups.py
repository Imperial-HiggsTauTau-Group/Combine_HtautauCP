import CombineHarvester.CombineTools.ch as ch
import argparse
import os

def get_args():
    description = '''This script can be used to break down the impact of different systematic uncertainties on the final measurement.'''
    parser = argparse.ArgumentParser(prog="breakdownSyst",description=description)
    parser.add_argument('--dir', '--directory', type=str, help="directory containing datacards")
    args = parser.parse_args()
    return args



def split_systematics(s):
    name = s.name()

    # protected parameters:

    special_systs = ['lumi','CMS_eff_m', 'CMS_eff_e', 'eff_t_syst_floor', 'CMS_HIG25012_eff_m_trigger', 'CMS_HIG25012_eff_e_trigger','CMS_scale_j_', 'CMS_res_j_', 'CMS_HIG25012_scale_met', 'CMS_HIG25012_res_met']

    if s.name() not in systematic_groups['all']:

        if True in [s.name().startswith(special) for special in special_systs]:
            return

        # print(f"Adding systematic: {name}")
        systematic_groups['all'].append(name)
        if "bbb" in name:
            systematic_groups['bbb'].append(name)
        elif ("CMS_" in name or 'lumi' in name or name.startswith("ff_") or name.startswith('SV_eff') or name.startswith('dy_pt')):
            systematic_groups['systematic'].append(name)
        elif "QCDscale" in name or "pdf_" in name or "BR_htt" in name or 'top_pt' in name or 'cross_section' in name or 'ps_isr' in name or 'ps_fsr' in name:
            systematic_groups['theory'].append(name)
        else:
            systematic_groups['other'].append(name)




def add_groups():

    print(f"Processing {len(all_cards)} datacards")

    for card_name in all_cards:

        global systematic_groups
        systematic_groups = {
            'theory': [],
            'systematic': [],
            'bbb': [],
            'other': [],
            'all': []
        }

        print(f'Processing datacard: {card_name}')
        cb = ch.CombineHarvester()
        cb.ParseDatacard(card_name, "", mass="125")

        cb.ForEachSyst(split_systematics)

        print(f'>> Identified total of {len(systematic_groups["all"])} systematics: {len(systematic_groups["systematic"])} experimental, {len(systematic_groups["bbb"])} bbb and {len(systematic_groups["theory"])} theory and {len(systematic_groups["other"])} other.')
        # theory_group = []
        if len(systematic_groups['other']) > 0:
            print("WARNING: There are systematics that were not categorized into any group")
            print(systematic_groups['other'])
            raise RuntimeError("Please categorize the above systematics into a group before proceeding.")


        if len(systematic_groups['theory']) > 0: cb.AddDatacardLineAtEnd(f"theory group = {' '.join(systematic_groups['theory'])}")
        if len(systematic_groups['systematic']) > 0: cb.AddDatacardLineAtEnd(f"experimental group = {' '.join(systematic_groups['systematic'])}")
        if len(systematic_groups['bbb']) > 0: cb.AddDatacardLineAtEnd(f"bbb group = {' '.join(systematic_groups['bbb'])}")


        # Write datacards
        print(">>> Writing datacards...")
        datacardtxt  = "%s/$TAG/$BIN.txt" % ('yr18_groups')
        datacardroot = "%s/$TAG/common/$BIN_input.root" % ('yr18_groups')
        writer = ch.CardWriter(datacardtxt,datacardroot)
        writer.SetVerbosity(1)
        writer.SetWildcardMasses([ ])
        writer.WriteCards("cmb", cb)

        print('-'*100)
        del cb




args = get_args()
directory = args.dir

print("Adding nuisance groups to datacards in directory:", directory)
run2_cards = []
run3_cards = []
for f in os.listdir(directory):
    if f.endswith('.txt'):
        if '2016' in f or '2017' in f or '2018' in f:
            run2_cards.append(os.path.join(directory, f))
        else:
            run3_cards.append(os.path.join(directory, f))
all_cards =  run3_cards + run2_cards
all_cards = sorted(all_cards)

add_groups()
