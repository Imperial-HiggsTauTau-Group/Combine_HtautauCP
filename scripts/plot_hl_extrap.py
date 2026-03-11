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
parser.add_argument('--run3_systs', type=str)
parser.add_argument('--yr18_systs', type=str)
parser.add_argument('--stat_only', type=str)

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


confidence = [0.683, 0.9545, 0.9973]#, 0.9999367, 0.9999994267]
levels = [ROOT.Math.chisquared_quantile_c(1 - cl, 1) for cl in confidence]


if args.run3_systs is not None:
    file_exp = ROOT.TFile.Open(os.path.join(args.directory, args.run3_systs))
    print(">> Opened run3_systs file:", os.path.join(args.directory, args.run3_systs))
    # Get run3_systs points
    key_exp = file_exp.GetListOfKeys().At(0)
    g_exp = key_exp.ReadObj()
    n_exp = g_exp.GetN()
    alpha_run3_systs = np.array([g_exp.GetX()[i] for i in range(n_exp)])
    nll_exp = np.array([g_exp.GetY()[i] for i in range(n_exp)])

    int_1sig_exp_exists = True if max(nll_exp) >= levels[0] else False
    int_2sig_exp_exists = True if max(nll_exp) >= levels[1] else False

    if int_1sig_exp_exists:
        int_1sig_exp = find_intersections(alpha_run3_systs, nll_exp, levels[0])
        print(f"1 sigma run3_systs intersections at alpha = {int_1sig_exp}")
    if int_2sig_exp_exists:
        int_2sig_exp = find_intersections(alpha_run3_systs, nll_exp, levels[1])
        print(f"2 sigma run3_systs intersections at alpha = {int_2sig_exp}")


if args.stat_only is not None:
    file_obs = ROOT.TFile.Open(os.path.join(args.directory, args.stat_only))
    print(">> Opened stat_only file:", os.path.join(args.directory, args.stat_only))
    # Get stat_only points
    key_obs = file_obs.GetListOfKeys().At(0)
    g_obs = key_obs.ReadObj()
    n_obs = g_obs.GetN()
    alpha_stat_only = np.array([g_obs.GetX()[i] for i in range(n_obs)])
    nll_obs = np.array([g_obs.GetY()[i] for i in range(n_obs)])

    int_1sig_obs_exists = True if max(nll_obs) >= levels[0] else False
    int_2sig_obs_exists = True if max(nll_obs) >= levels[1] else False

    if int_1sig_obs_exists:
        int_1sig_obs = find_intersections(alpha_stat_only, nll_obs, levels[0])
        print(f"1 sigma stat_only intersections at alpha = {int_1sig_obs}")
    if int_2sig_obs_exists:
        int_2sig_obs = find_intersections(alpha_stat_only, nll_obs, levels[1])
        print(f"2 sigma stat_only intersections at alpha = {int_2sig_obs}")

if args.yr18_systs is not None:
    file_yr18_systs = ROOT.TFile.Open(os.path.join(args.directory, args.yr18_systs))
    print(">> Opened yr18_systs file:", os.path.join(args.directory, args.yr18_systs))
    # Get yr18_systs points
    key_yr18_systs = file_yr18_systs.GetListOfKeys().At(0)
    g_yr18_systs = key_yr18_systs.ReadObj()
    n_yr18_systs = g_yr18_systs.GetN()
    alpha_yr18_systs = np.array([g_yr18_systs.GetX()[i] for i in range(n_yr18_systs)])
    nll_yr18_systs = np.array([g_yr18_systs.GetY()[i] for i in range(n_yr18_systs)])

    int_1sig_yr18_systs_exists = True if max(nll_yr18_systs) >= levels[0] else False
    int_2sig_yr18_systs_exists = True if max(nll_yr18_systs) >= levels[1] else False

    if int_1sig_yr18_systs_exists:
        int_1sig_yr18_systs = find_intersections(alpha_yr18_systs, nll_yr18_systs, levels[0])
        print(f"1 sigma yr18_systs intersections at alpha = {int_1sig_yr18_systs}")




fig, ax = plt.subplots(figsize=(7.5,6))
# add CLs
ax.axhline(levels[0], color='grey', linestyle='-', linewidth=0.75)
ax.text(-15, levels[0]+0.1, '68.3%', color='grey', fontsize=14, verticalalignment='bottom', horizontalalignment='right')
ax.axhline(levels[1], color='grey', linestyle='-', linewidth=0.75)
ax.text(-15, levels[1]+0.1, '95.5%', color='grey', fontsize=14,  verticalalignment='bottom', horizontalalignment='right')
ax.axhline(levels[2], color='grey', linestyle='-', linewidth=0.75)
ax.text(-15, levels[2]+0.1, '99.7%', color='grey', fontsize=14, verticalalignment='bottom', horizontalalignment='right')


if args.run3_systs is not None:
    bestfit_exp = alpha_run3_systs[np.argmin(nll_exp)]
    print(f"Best fit run3_systs alpha: {bestfit_exp}")
    # intersections and 1 and 2 sigma
    if int_1sig_exp_exists:
        label_string_run3_systs = rf"With Run-3 syst. uncert.: $\alpha^{{H\tau\tau}} = {round(np.abs(bestfit_exp),0):.0f} \pm {round((np.abs(bestfit_exp-int_1sig_exp[1])+np.abs(bestfit_exp-int_1sig_exp[0]))/2,0):.0f} ^\circ$"
    else:
        label_string_run3_systs = rf"With Run-3 syst. uncert.: $\alpha^{{H\tau\tau}} = {round(np.abs(bestfit_exp),0):.0f} ^\circ$"



if args.stat_only is not None:
    # add NLL curve
    bestfit_obs = alpha_stat_only[np.argmin(nll_obs)]
    print(f"Best fit stat_only alpha: {bestfit_obs}")

    # intersections and 1 and 2 sigma
    if int_1sig_obs_exists:
        label_string_stat_only = rf"With Stat. uncert. only: $\alpha^{{H\tau\tau}} = {round(bestfit_obs,0):.0f} \pm {round((np.abs(bestfit_obs-int_1sig_obs[1])+np.abs(bestfit_obs-int_1sig_obs[0]))/2,0):.0f} ^\circ$"
    else:
        label_string_stat_only = rf"With Stat. uncert. only: $\alpha^{{H\tau\tau}} = {round(bestfit_obs,0):.0f} ^\circ$"

if args.yr18_systs is not None:
    bestfit_yr18_systs = alpha_yr18_systs[np.argmin(nll_yr18_systs)]
    print(f"Best fit yr18_systs alpha: {bestfit_yr18_systs}")
    # intersections and 1 and 2 sigma
    if int_1sig_yr18_systs_exists:
        label_string_yr18_systs = rf"With YR18 syst. uncert.: $\alpha^{{H\tau\tau}} = {round(np.abs(bestfit_yr18_systs),0):.0f} \pm {round((np.abs(bestfit_yr18_systs-int_1sig_yr18_systs[1])+np.abs(bestfit_yr18_systs-int_1sig_yr18_systs[0]))/2,0):.0f} ^\circ$"
    else:
        label_string_yr18_systs = rf"With YR18 syst. uncert.: $\alpha^{{H\tau\tau}} = {round(np.abs(bestfit_yr18_systs),0):.0f}$ ^\circ"

if args.stat_only is not None:
    ax.plot(alpha_stat_only, nll_obs, linestyle='-', color='red', label=label_string_stat_only)
if args.run3_systs is not None:
    ax.plot(alpha_run3_systs, nll_exp, linestyle='-', color='darkblue', label=label_string_run3_systs)
if args.yr18_systs is not None:
    ax.plot(alpha_yr18_systs, nll_yr18_systs, linestyle='-', color='black', label=label_string_yr18_systs)

ax.set_ylim(0, 40)

ax.set_xlabel(r'$\alpha^{H\tau\tau}$ (degrees)')
ax.set_ylabel(r'-2$\Delta$lnL')

ax.set_xlim(-20, 20)
fig.tight_layout(pad=1.2)
plt.legend(frameon=True, loc='upper center', fontsize=15.5)

#hep.cms.label(ax=ax, label="Preliminary", data=True, lumi="Projection, 3", com='13.6', fontsize=18)
hep.cms.label(
    ax=ax,
    label="Supplementary",
    data=True,
    rlabel=r"$Projection, 3\,\mathrm{ab}^{-1}$ ($13.6$ TeV)",
    fontsize=18,
)

plt.savefig(os.path.join(args.directory, f'alpha_HL_extrap.pdf'))

