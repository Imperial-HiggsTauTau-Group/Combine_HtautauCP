import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import argparse
import os
import math
import scipy.stats

hep.style.use("CMS")
plt.rcParams.update({"font.size": 20})

parser = argparse.ArgumentParser()
parser.add_argument('--directory', type=str, required=True, help='Base directory containing {tauID}_{WP} subdirectories')
parser.add_argument('--tauID', type=str, required=True, help='Tau ID name (e.g. DeepTau, PNet)')
parser.add_argument('--channel', type=str, default='cmb', help='Channel')

args = parser.parse_args()


def find_intersections(x, y, ylevel):
    intersections_x = []
    for i in range(len(x) - 1):
        if (y[i] - ylevel) * (y[i+1] - ylevel) < 0:
            dx = x[i+1] - x[i]
            dy = y[i+1] - y[i]
            x_int = x[i] + (ylevel - y[i]) / (dy / dx)
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
else:
    ch_label = r''

wp_colors = {
    'Medium': 'steelblue',
    'Tight': 'darkorange',
    'VTight': 'forestgreen',
}

wp_name_map ={
    5: 'Medium',
    6: 'Tight',
    7: 'VTight',
}

working_points = [5, 6, 7]  # Corresponding to Medium, Tight, VTight (change if needed)

fig, ax = plt.subplots(figsize=(7.5, 6))
ax.axhline(levels[0], color='grey', linestyle='-', linewidth=0.75)
ax.text(-70, levels[0] + 0.1, '68.3%', color='grey', fontsize=14,
        verticalalignment='bottom', horizontalalignment='right')
ax.axhline(levels[1], color='grey', linestyle='-', linewidth=0.75)
ax.text(-70, levels[1] + 0.1, '95.5%', color='grey', fontsize=14,
        verticalalignment='bottom', horizontalalignment='right')
# find the uncertainties and plot curves
all_nll_max = []


for idx, wp_ID in enumerate(working_points):
    wp = wp_name_map.get(wp_ID)
    color = wp_colors.get(wp)
    exp_path = os.path.join(args.directory, f"{args.tauID}_{wp_ID}", channel, f"alpha_{channel}_EXPECTED.root")

    f = ROOT.TFile.Open(exp_path)
    if not f or f.IsZombie():
        print(f"Warning: could not open {exp_path}, skipping {wp}")
        continue
    print(f">> Opened {exp_path}")

    g = f.GetListOfKeys().At(0).ReadObj()
    n = g.GetN()
    alpha = np.array([g.GetX()[j] for j in range(n)])
    nll   = np.array([g.GetY()[j] for j in range(n)])
    all_nll_max.append(max(nll))

    bestfit = alpha[np.argmin(nll)]
    print(f"Bestfit expected alpha ({wp}): {bestfit:.1f}")

    nll_at_90 = np.interp(90.0, alpha, nll)
    delta_nll_ps = nll_at_90 - min(nll)
    ps_sig = math.sqrt(scipy.stats.chi2.ppf(scipy.stats.chi2.cdf(delta_nll_ps, 1), 1))
    print(f"CP-odd exclusion significance ({wp}): {ps_sig:.2f}σ")

    int_1sig_exists = max(nll) >= levels[0]
    if int_1sig_exists:
        int_1sig = find_intersections(alpha, nll, levels[0])
        print(f">> 1-sigma expected intersections ({wp}): {int_1sig}")
        error_up   = round(np.abs(bestfit - int_1sig[1]), 0)
        error_down = round(np.abs(bestfit - int_1sig[0]), 0)
        if error_up == error_down:
            label = rf"{wp}: $\alpha^{{H\tau\tau}} = (0\pm{error_up:.0f})^\circ$, CP-odd excl. {ps_sig:.2f}$\sigma$"
        else:
            label = rf"{wp}: $\alpha^{{H\tau\tau}} = (0^{{+{error_up:.0f}}}_{{-{error_down:.0f}}})^\circ$, CP-odd excl. {ps_sig:.2f}$\sigma$"
    else:
        label = rf"{wp}: $\alpha^{{H\tau\tau}} = {round(np.abs(bestfit), 0):.0f}^\circ$, CP-odd excl. {ps_sig:.2f}$\sigma$"

    ax.plot(alpha, nll, linestyle='--', color=color, label=label)

ax.set_ylim(0, max(all_nll_max) * 1.4 if all_nll_max else 10)
ax.set_xlim(-90, 90)
ax.set_xlabel(r'$\alpha^{H\tau\tau}$ (degrees)')
ax.set_ylabel(r'-2$\Delta$lnL')

fig.tight_layout(pad=1.2)
plt.legend(frameon=False, loc='upper right', fontsize=14)

algo_label = 'DeepTau' if args.tauID.lower() == 'deeptau2018v2p5' else args.tauID

ax.text(0.03, 0.95, ch_label, transform=ax.transAxes, fontsize=18,
        fontweight='bold', verticalalignment='top', horizontalalignment='left')
ax.text(0.03, 0.88, algo_label, transform=ax.transAxes, fontsize=17,
        fontweight='bold', verticalalignment='top', horizontalalignment='left')

hep.cms.label(ax=ax, label="Work in progress", data=True, lumi='219.66', com='13.6', fontsize=18)

outfile = os.path.join(args.directory, f'alpha_{args.tauID}_{channel}_EXP.pdf')
plt.savefig(outfile, bbox_inches='tight')
print(f"Saved: {outfile}")
