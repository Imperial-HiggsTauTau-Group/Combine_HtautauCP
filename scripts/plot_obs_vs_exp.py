import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep as hep
import argparse
import os

hep.style.use("CMS")
plt.rcParams.update({"font.size": 20})

parser = argparse.ArgumentParser()
parser.add_argument('--combination', action='store_true')

parser.add_argument('--directory', type=str, help='Directory')
parser.add_argument('--expected', type=str)
parser.add_argument('--observed', type=str)
parser.add_argument('--channel', type=str, help='Channel', default ='cmb')
parser.add_argument('--extra_name', type=str, help='Extra name to append to output name', default ='')


args = parser.parse_args()


def find_intersections(x, y, ylevel):
    intersections_x = []
    for i in range(len(x)-1):
        if (y[i] - ylevel)*(y[i+1] - ylevel) < 0: # find point where curve crosses ylevel
            # linear interpolation to estimate exact crossing
            dx = x[i+1] - x[i]
            dy = y[i+1] - y[i]
            slope = dy/dx
            x_int = x[i] + (ylevel - y[i])/slope
            intersections_x.append(x_int)
    return np.array(intersections_x)


confidence = [0.683, 0.9545, 0.9973]
levels = [ROOT.Math.chisquared_quantile_c(1 - cl, 1) for cl in confidence]


channel = args.channel

if 'cmb' in channel:
    ch_label = r'$\tau_{h}\tau_{h}/\mu\tau_{h}/e\tau_{h}$'
elif 'tt' in channel:
    ch_label = r'$\tau_{h}\tau_{h}$'
elif 'mt' in channel:
    ch_label = r'$\mu\tau_{h}$'
elif 'et' in channel:
    ch_label = r'$e\tau_{h}$'
elif 'a1a1' in channel:
    ch_label = r'$a_{1}^{3pr}a_{1}^{3pr}$'
elif 'a11pra1' in channel:
    ch_label = r'$a_{1}^{1pr}a_{1}^{3pr}$'
elif 'pia11pr' in channel:
    ch_label = r'$\pi a_{1}^{1pr}$'
elif 'pirho' in channel:
    ch_label = r'$\pi \rho$'
elif 'rhoa11pr' in channel:
    ch_label = r'$\rho a_{1}^{1pr}/a_{1}^{1pr}a_{1}^{1pr}$'
elif 'pia1' in channel:
    ch_label = r'$\pi a_{1}^{3pr}$'
elif 'pipi' in channel:
    ch_label = r'$\pi \pi$'
elif 'rhorho' in channel:
    ch_label = r'$\rho \rho$'
elif 'rhoa1' in channel:
    ch_label = r'$\rho a_{1}^{3pr}$'
elif 'ea11pr' in channel:
    ch_label = r'$e a_{1}^{1pr}$'
elif 'erho' in channel:
    ch_label = r'$e \rho$'
elif 'ea1' in channel:
    ch_label = r'$e a_{1}^{3pr}$'
elif 'epi' in channel:
    ch_label = r'$e \pi$'
elif 'mua11pr' in channel:
    ch_label = r'$\mu a_{1}^{1pr}$'
elif 'murho' in channel:
    ch_label = r'$\mu \rho$'
elif 'mua1' in channel:
    ch_label = r'$\mu a_{1}^{3pr}$'
elif 'mupi' in channel:
    ch_label = r'$\mu \pi$'
else:
    ch_label = r''


if args.expected is not None:
    file_exp = ROOT.TFile.Open(os.path.join(args.directory, args.expected))
    print(">> Opened expected file:", os.path.join(args.directory, args.expected))
    # Get expected points
    key_exp = file_exp.GetListOfKeys().At(0)
    g_exp = key_exp.ReadObj()
    n_exp = g_exp.GetN()
    alpha_exp = np.array([g_exp.GetX()[i] for i in range(n_exp)])
    nll_exp = np.array([g_exp.GetY()[i] for i in range(n_exp)])

    int_1sig_exp_exists = True if max(nll_exp) >= levels[0] else False
    int_2sig_exp_exists = True if max(nll_exp) >= levels[1] else False

    if int_1sig_exp_exists:
        int_1sig_exp = find_intersections(alpha_exp, nll_exp, levels[0])
        print(f"1 sigma expected intersections at alpha = {int_1sig_exp}")
    if int_2sig_exp_exists:
        int_2sig_exp = find_intersections(alpha_exp, nll_exp, levels[1])
        print(f"2 sigma expected intersections at alpha = {int_2sig_exp}")


if args.observed is not None:
    file_obs = ROOT.TFile.Open(os.path.join(args.directory, args.observed))
    print(">> Opened observed file:", os.path.join(args.directory, args.observed))
    # Get observed points
    key_obs = file_obs.GetListOfKeys().At(0)
    g_obs = key_obs.ReadObj()
    n_obs = g_obs.GetN()
    alpha_obs = np.array([g_obs.GetX()[i] for i in range(n_obs)])
    nll_obs = np.array([g_obs.GetY()[i] for i in range(n_obs)])

    int_1sig_obs_exists = True if max(nll_obs) >= levels[0] else False
    int_2sig_obs_exists = True if max(nll_obs) >= levels[1] else False

    if int_1sig_obs_exists:
        int_1sig_obs = find_intersections(alpha_obs, nll_obs, levels[0])
        print(f"1 sigma observed intersections at alpha = {int_1sig_obs}")
    if int_2sig_obs_exists:
        int_2sig_obs = find_intersections(alpha_obs, nll_obs, levels[1])
        print(f"2 sigma observed intersections at alpha = {int_2sig_obs}")




fig, ax = plt.subplots(figsize=(7.5,6))
# add CLs
ax.axhline(levels[0], color='grey', linestyle='-', linewidth=0.75)
ax.axhline(levels[1], color='grey', linestyle='-', linewidth=0.75)
if args.expected is not None:
    bestfit_exp = alpha_exp[np.argmin(nll_exp)]
    print(f"Best fit expected alpha: {bestfit_exp}")
    # intersections and 1 and 2 sigma
    if int_1sig_exp_exists:
        ax.plot([int_1sig_exp[0], int_1sig_exp[0]], [0, levels[0]], color='grey', linestyle='--', linewidth=0.75)
        ax.plot([int_1sig_exp[1], int_1sig_exp[1]], [0, levels[0]], color='grey', linestyle='--', linewidth=0.75)
        label_string_exp = rf"Expected: $\alpha^{{H\tau\tau}} = {round(np.abs(bestfit_exp),0):.0f}^{{+{round(np.abs(bestfit_exp-int_1sig_exp[1]),0):.0f}}}_{{-{round(np.abs(bestfit_exp-int_1sig_exp[0]),0):.0f}}}$"
    else:
        label_string_exp = rf"Expected: $\alpha^{{H\tau\tau}} = {round(np.abs(bestfit_exp),0):.0f}$"
    if int_2sig_exp_exists:
        ax.plot([int_2sig_exp[0], int_2sig_exp[0]], [0, levels[1]], color='grey', linestyle='--', linewidth=0.75)
        ax.plot([int_2sig_exp[1], int_2sig_exp[1]], [0, levels[1]], color='grey', linestyle='--', linewidth=0.75)



if args.observed is not None:
    # add NLL curve
    bestfit_obs = alpha_obs[np.argmin(nll_obs)]
    print(f"Best fit observed alpha: {bestfit_obs}")

    # intersections and 1 and 2 sigma
    if int_1sig_obs_exists:
        ax.plot([int_1sig_obs[0], int_1sig_obs[0]], [0, levels[0]], color='grey', linestyle='-', linewidth=0.75)
        ax.plot([int_1sig_obs[1], int_1sig_obs[1]], [0, levels[0]], color='grey', linestyle='-', linewidth=0.75)
        label_string_obs = rf"Observed: $\alpha^{{H\tau\tau}} = {round(bestfit_obs,0):.0f}^{{+{round(np.abs(bestfit_obs-int_1sig_obs[1]),0):.0f}}}_{{-{round(np.abs(bestfit_obs-int_1sig_obs[0]),0):.0f}}}$"
    else:
        label_string_obs = rf"Observed: $\alpha^{{H\tau\tau}} = {round(bestfit_obs,0):.0f}$"
    if int_2sig_obs_exists:
        ax.plot([int_2sig_obs[0], int_2sig_obs[0]], [0, levels[1]], color='grey', linestyle='-', linewidth=0.75)
        ax.plot([int_2sig_obs[1], int_2sig_obs[1]], [0, levels[1]], color='grey', linestyle='-', linewidth=0.75)


if args.observed is not None:
    ax.plot(alpha_obs, nll_obs, linestyle='-', color='darkblue', label=label_string_obs)
if args.expected is not None:
    ax.plot(alpha_exp, nll_exp, linestyle='--', color='gray', label=label_string_exp)

if args.expected is None:
    ax.set_ylim(0, max(nll_obs)*1.4)
elif args.observed is None:
    ax.set_ylim(0, max(nll_exp)*1.4)
else:
    ax.set_ylim(0, max([max(nll_obs), max(nll_exp)])*1.4)

ax.set_xlabel(r'$\alpha^{H\tau\tau}$')
ax.set_ylabel(r'-2$\Delta$lnL')

ax.set_xlim(-90, 90)
fig.tight_layout(pad=1.2)
plt.legend(frameon=True, loc='upper right')
plt.title(ch_label, fontsize=18)
if args.combination:
    hep.cms.label(ax=ax, label="Preliminary", data=True, lumi='200', com='13 and 13.6', fontsize=18)
else:
    hep.cms.label(ax=ax, label="Preliminary", data=True, lumi='62.4', com='13.6', fontsize=18)
if args.expected is None:
    plt.savefig(os.path.join(args.directory, args.observed.replace('.root','_formatted.pdf')))
elif args.observed is None:
    plt.savefig(os.path.join(args.directory, args.expected.replace('.root','_formatted.pdf')))
else:
    plt.savefig(os.path.join(args.directory, f'alpha_OBS_vs_EXP_{args.channel}_{args.extra_name}.pdf'))
