
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

parser = argparse.ArgumentParser(description="Script to plot the correlation matrix of the fit parameters from a MultiDimFit result.")
parser.add_argument("--root_file", "-f", required=True, help="Path to the ROOT file containing the MultiDimFit result.")
parser.add_argument("--output_dir", "-o", required=True, help="Directory where the output plot will be saved.")
args = parser.parse_args()

ch = "tt"

file = args.root_file
fit_name = "fit_mdf"

f = ROOT.TFile.Open(file)
fit_result = f.Get(fit_name)


fontsize_map = {
    'values': 30,
    'cbar_label': 30,
    'cbar_tick_label': 30,
    'x_tick_label': 30,
    'y_tick_label': 30
    # 'x_tick_label': 6,
    # 'y_tick_label': 6
}

# get covariance matrix
cov_matrix = fit_result.covarianceMatrix()
n = cov_matrix.GetNrows()
corr_matrix = np.zeros((n, n))
print(f'Loaded matrix with {n} rows')

# convert to correlation matrix
for i in range(n):
    for j in range(n):
        corr_matrix[i,j]=cov_matrix[i,j]/np.sqrt(cov_matrix[i,i]*cov_matrix[j,j])

print('-> Successfully converted covariance to correlation matrix')


# get nuisance parameters
nuisance_names = [p.GetName() for p in fit_result.floatParsFinal() ]
print(f"Retried {len(nuisance_names)} parameter names")

print(corr_matrix.shape)

indices_to_drop = []
kept_params = []

# Drop nuisancesthat are not of interest (e.g. bbb nuisances) to make the matrix more readable
for i, param in enumerate(nuisance_names):
    if param not in ['alpha', 'muggH', 'muV']:
        print("Dropping nuisance:", param)
        indices_to_drop.append(i)
    else:
        kept_params.append(param)

corr_matrix = np.delete(corr_matrix, indices_to_drop, axis=0)
corr_matrix = np.delete(corr_matrix, indices_to_drop, axis=1)

print("Slimmed matrix to size:", corr_matrix.shape)

# Function to plot numbers on matrix
def plot_numbers_on_matrix(ax,mat,fontsize=16,drop_minor=False):
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            c = mat[j,i]
            if drop_minor:
                if abs(c) > 0.005:
                    ax.text(i,j,'{:.3f}'.format(c),fontdict={'size': fontsize},va='center',ha='center')
            else:
                if c != 0:
                    ax.text(i,j,'{:.3f}'.format(c),fontdict={'size': fontsize},va='center',ha='center')



def plot_corr_matrix( fig, ax, corr_matrix, axis_labels=[], do_text=False, drop_minor=False):
    plt.set_cmap("bwr")
    mat = ax.matshow( corr_matrix, vmin=-1, vmax=1 )
    divider = make_axes_locatable(ax)
    if do_text: plot_numbers_on_matrix(ax, corr_matrix, fontsize=fontsize_map['values'], drop_minor=drop_minor)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(mat, cax=cax)
    cbar.set_label("$\\rho$", fontsize=fontsize_map['cbar_label'])
    cbar.ax.tick_params(labelsize=fontsize_map['cbar_tick_label'])
    ax.set_xticks(np.arange(len(axis_labels)))
    ax.set_yticks(np.arange(len(axis_labels)))
    ax.set_xticklabels( axis_labels, fontsize=fontsize_map['x_tick_label'], rotation=90 )
    ax.set_yticklabels( axis_labels, fontsize=fontsize_map['y_tick_label'] )
    ax.tick_params(axis='x', which='both', size=0, pad=1)
    ax.tick_params(axis='y', which='both', size=0, pad=1)
    ax.xaxis.tick_bottom()


fig, ax = plt.subplots(figsize=(20,20))
plot_corr_matrix( fig, ax, corr_matrix, axis_labels=kept_params, do_text=True, drop_minor=True)
fig.savefig(
    f'{args.output_dir}/correlation_matrix_{ch}.pdf',
    dpi=600)


