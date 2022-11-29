import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import subprocess
import pandas as pd
from numpy import linalg as LA
import lmfit
import shutil
import re
import os
from sys import platform
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Open Sans']
from matplotlib.ticker import FormatStrFormatter

font = {'family': 'normal',
        'weight': 'normal',
        'size': 20}

matplotlib.rc('font', **font)

sns.set_style('white')
sns.set_context('notebook')
sns.set(font='Open Sans')

# Variblen Dictionary --> zentrale änderung von Dateinamen und Konstruktnamen
variables_dict = {
    'name': 'BTL',
    'exp_FRET': 'raw/FRET_40_20.csv',
    'bins': 25
}


def check_os():
    """
    Print's the OS you are worling on
    :return: None
    """
    if platform == "linux" or platform == "linux2":
        print("You are on linux")
    elif platform == "darwin":
        print("Your are on MacOS")


# Farbumrechner um die RGB Farben aus dem Colorbrewer und Paleton direkt einzutragen
def c_c(r, g, b):
    """
    Function to convert rgb color values from e.g. Colobrewer or Paleton in a matplotl
    :param r:
    :param g:
    :param b:
    :return:
    """
    return [round(r / 255, 2), round(g / 255, 2), round(b / 255, 2)]


def make_result_dir(directory_name):
    mkdirpath = f"{os.getcwd()}/{directory_name}"
    try:
        os.mkdir(mkdirpath)
    except FileExistsError:
        print(f"Results can be found in: {mkdirpath}")
    except OSError:
        print(f"Failed to create {mkdirpath}")
    else:
        print(f"Successfully created the directory {mkdirpath}. Results can be found there")


def read_xvg(fname):
    """Read columns of data from file fname

    Returns numpy array of data
    """
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if not line.startswith(('#', '@')):
                skip = i
                break

    return np.genfromtxt(fname, skip_header=skip)


# Tickstyle definitionen
def set_ticksStyle(x_size=30, y_size=30, x_dir='in', y_dir='in'):
    sns.set_style('ticks', {'xtick.major.size': x_size, 'ytick.major.size': y_size, 'xtick.direction': x_dir,
                            'ytick.direction': y_dir})


# Function for kappa2 calculation in python
def get_atoms_coordinates(atom_id, universe):
    universe.trajectory[0]
    acc_c14 = universe.select_atoms(f'id {atom_id}')
    acc_coords = []
    for ts in universe.trajectory:
        coods = list(np.squeeze(acc_c14.positions / 10))
        acc_coords.append(coods)
    return np.array(acc_coords)


def calculate_kappa_2(don_dipol: list, acc_dipol: list) -> list:
    """
    Funktion zur Berechnung von kappa^2 aus den Koordinaten des Dipolvektors, mit zwei gegebenen Atomen

    :param don_dipol: List --> [[[x, v, z],..., [x, v, z]], [[x, v, z],...,[x, v, z]]] Array mit den Koordinaten der Trajektorie die den Dipolvektor des Donorfarbstoffs definieren
    :param acc_dipol: List --> [[[x, v, z],..., [x, v, z]], [[x, v, z],...,[x, v, z]]] Array mit den Koordinaten der Trajektorie die den Dipolvektor des Acceptorfarbstoffs definieren
    :return: Numpy array mit kappa^2 Werten
    """
    # Initilizing vectors woth zero
    dvect = np.zeros([len(don_dipol[1]), 3])
    avect = np.zeros([len(acc_dipol[1]), 3])
    dpos = np.zeros([len(don_dipol[1]), 3])
    apos = np.zeros([len(acc_dipol[1]), 3])
    # Vektoren
    for i in range(0, int(len(don_dipol) / 2)):
        # Vector von einem Atom zum Anderen = Richtungsvektor
        dvect = dvect - don_dipol[2 * i] + don_dipol[2 * i + 1]
        dpos = dpos + don_dipol[2 * i] + don_dipol[2 * i + 1]

    for i in range(0, int(len(acc_dipol) / 2)):
        # Vector von einem Atom zum Anderen = Richtungsvektor
        avect = avect - acc_dipol[2 * i] + acc_dipol[2 * i + 1]
        apos = apos + acc_dipol[2 * i] + acc_dipol[2 * i + 1]
        # Richtungsvektor bestimmen

    # Vektoren Normieren
    # Euklidische Normierung
    dvect = np.divide(dvect, np.expand_dims(LA.norm(dvect, axis=1), axis=1))
    avect = np.divide(avect, np.expand_dims(LA.norm(avect, axis=1), axis=1))

    dpos = 1 / len(don_dipol) * dpos
    apos = 1 / len(acc_dipol) * apos

    # Vektor zwischen den Mittelpunkten der Farbstoffe
    dist = dpos - apos
    distnorm = np.divide(dist, np.expand_dims(LA.norm(dist, axis=1), axis=1))

    a = np.sum(dvect * avect, axis=1)
    b = np.sum(dvect * distnorm, axis=1)
    c = np.sum(distnorm * avect, axis=1)

    kappa = a - 3 * b * c
    kappa = np.around(kappa ** 2, 7)
    return kappa


# Function for dye distance calculation in python
def calculate_inter_dye_distance(mean_donor_atom: list, mean_acceptor_atom: list) -> list:
    """
    Funktion zur Berechnung der Distanz zweier Atome

    :param mean_donor_atom: List --> Trajektorie der Koordinaten des mittleren C Atoms des Donorfarbstoffes
    :param mean_acceptor_atom: List --> Trajektorie der Koordinaten des mittleren C Atoms des Acceptorfarbstoffes
    :return: Liste der Distanzen in Angstrom
    """
    return np.round(np.sqrt(np.sum((np.subtract(mean_donor_atom, mean_acceptor_atom)) ** 2, axis=1)), 7)


# Fit für die Gaussfunktion der FRET Effizienz
def burst_fit(E_FRET, bins_i):
    hist, bins = np.histogram(E_FRET, bins=bins_i, range=(0, 1), weights=np.ones(len(E_FRET)) / len(E_FRET))
    bincenters = binMid = (bins[1:] + bins[:-1]) / 2

    mod = lmfit.models.GaussianModel()
    pars = mod.guess(hist, x=bincenters)
    out = mod.fit(hist, pars, x=bincenters)
    x_fit = np.linspace(0, 1, 100)
    y_fit = mod.func(x_fit, *list(out.params.valuesdict().values())[:3])

    return x_fit, y_fit


# Funktion für die Visualisierung von FRET
def vis_exp_md_acv(md_e, acv_e, exp_e, name, acv_means, show_md=True, show_acv=True, show_exp=True):
    md_fit_x = burst_fit(md_e.FRETefficiencies, variables_dict['bins'])
    acv_fit_x = burst_fit(acv_e.FRETefficiencies, variables_dict['bins'])
    exp_fit_x = burst_fit(exp_e, variables_dict['bins'])
    with sns.axes_style('ticks'):
        set_ticksStyle()
        f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)

        # Plot Experiment und Verteilung Simulation
        if show_md:
            print(f"⟨E_MD⟩= {round(np.mean(md_e.FRETefficiencies), 2)} ± {round(np.std(md_e.FRETefficiencies), 2)}")
            ax[0, 0].hist(md_e.FRETefficiencies, bins=variables_dict['bins'], range=[0, 1], color=c_c(72, 92, 120),
                          edgecolor='w', linewidth=0.5, label='$\mathregular{{\it E}_{MD}}$', density=False,
                          weights=np.ones(len(md_e.FRETefficiencies)) / len(md_e.FRETefficiencies), stacked=False,
                          zorder=3)
            ax[0, 0].plot(md_fit_x[0], md_fit_x[1], color=c_c(52, 72, 100), linewidth=2, zorder=4)
            # ax[0,0].axvline(x=np.mean(md_e.FRETefficiencies), color=c_c(72, 92, 120))

        if show_acv:
            ax[0, 0].hist(acv_e.FRETefficiencies, bins=variables_dict['bins'], range=[0, 1], color=[0.75, 0.51, 0.38],
                          edgecolor='w', linewidth=2, density=False,
                          weights=np.ones(len(acv_e.FRETefficiencies)) / len(acv_e.FRETefficiencies), stacked=False,
                          zorder=1, label='$\mathregular{{\it E}_{MACV}}$')
            ax[0, 0].plot(acv_fit_x[0], acv_fit_x[1], color=c_c(152, 103, 72), linewidth=2, zorder=2)
            # ax[0,0].axvline(x=acv_means, color=c_c(152, 103, 72))
            # ax[0,0].text(x=acv_means, y=1.05, s=r"$⟨\mathregular{{\it E}_{MACV}}⟩$ without shot noise", horizontalalignment='center',verticalalignment='center', transform=ax[0,0].transAxes)

        if show_exp:
            ax[0, 0].hist(exp_e, bins=variables_dict['bins'], range=[0, 1], alpha=0.7, color=c_c(182, 152, 102),
                          edgecolor='w', linewidth=2, label=r"$\mathregular{{\it E}_{exp.}}$", density=False,
                          weights=np.ones(len(exp_e)) / len(exp_e), stacked=False, zorder=5)
            ax[0, 0].plot(exp_fit_x[0], exp_fit_x[1], color=c_c(148, 118, 68), linewidth=2, zorder=6)

        ax[0, 0].set_xlabel('FRET', fontsize=30)
        ax[0, 0].set_ylabel('probability', fontsize=30)
        ax[0, 0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax[0, 0].spines['top'].set_visible(False)
        ax[0, 0].spines['right'].set_visible(False)
        ax[0, 0].spines['bottom'].set_visible(True)
        ax[0, 0].spines['left'].set_visible(True)
        ax[0, 0].locator_params(axis="y", nbins=4)
        ax[0, 0].locator_params(axis="x", nbins=4)
        ax[0, 0].tick_params(axis='both', labelsize=20)
        # plt.xticks(fontsize=14)
        # plt.yticks(fontsize=14)
        ax[0, 0].set_ylim(0, 0.33)
        ax[0, 0].legend(frameon=False, handlelength=0.75, fontsize=20)
        plt.savefig(f"{name}.svg", dpi=300, bbox_inches="tight")
        plt.show()


# Funktion zur Visualisierung von rda
def vis_rda_acv_md(rda_acv: list, rda_md: list, export_file: str, range: list):
    with sns.axes_style('ticks'):
        set_ticksStyle()
        f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)
        ax[0, 0].hist(rda_md, bins=30, range=range, alpha=0.5, color=c_c(52, 72, 100), edgecolor='w', linewidth=0.5,
                      zorder=1, weights=np.ones(len(rda_md)) / len(rda_md))
        ax[0, 0].hist(rda_acv, bins=30, range=range, alpha=1, color=[0.75, 0.51, 0.38], edgecolor='w', linewidth=0.5,
                      zorder=0, weights=np.ones(len(rda_acv)) / len(rda_acv))

        ax[0, 0].set_xlabel('distance $R$ ($\mathregular{\AA}$)', fontsize=30)
        ax[0, 0].set_ylabel('probability', fontsize=30)
        ax[0, 0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        # plt.xticks(fontsize=14)
        # plt.yticks(fontsize=14)
        ax[0, 0].legend(['$\mathregular{{\it R}_{DA-MD}}$', '$\mathregular{{\it R}_{DA-MACV}}$'], frameon=False,
                        handlelength=1, fontsize=20)
        ax[0, 0].spines['top'].set_visible(False)
        ax[0, 0].spines['right'].set_visible(False)
        ax[0, 0].spines['bottom'].set_visible(True)
        ax[0, 0].spines['left'].set_visible(True)
        ax[0, 0].locator_params(axis="y", nbins=4)
        ax[0, 0].locator_params(axis="x", nbins=4)
        ax[0, 0].tick_params(axis='both', labelsize=20)

        plt.savefig(f"{export_file}.svg", dpi=300, bbox_inches="tight")
        plt.show()


# Funktion zur Visualisierung von kapppa2
def vis_kappa_md(kappa2, export_file: str, yp: float = 0.3):
    print(kappa2.mean())
    with sns.axes_style('ticks'):
        set_ticksStyle()
        f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)
        ax[0, 0].hist(kappa2, bins=variables_dict['bins'], range=[0, 4], alpha=0.5, color=c_c(182, 152, 102),
                      edgecolor='w', linewidth=0.5, zorder=1, weights=np.ones(len(kappa2)) / len(kappa2))
        ax[0, 0].axvline(x=kappa2.mean(), ymin=0.9, ymax=1, color=c_c(72, 92, 120))
        ax[0, 0].text(x=kappa2.mean(), y=yp, s=f"{round(kappa2.mean(), 2)}", color=c_c(72, 92, 120))
        ax[0, 0].set_xlabel('$\kappa^2$')
        ax[0, 0].set_ylabel('probability')
        plt.savefig(f"{export_file}.png", dpi=300, bbox_inches="tight")
        plt.show()


# Plots time trajectory and historgram together
def vis_plot_hist(x, y, ylabel, export_file, color, range):
    with sns.axes_style('ticks'):
        set_ticksStyle()
        # start with a square Figure
        fig = plt.figure(figsize=(7, 3))

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        spacing = 0.025

        rect_scatter = [left, bottom, width, height]
        rect_histy = [left + width + spacing, bottom, 0.2, height]
        ax = fig.add_axes(rect_scatter)
        ax_histy = fig.add_axes(rect_histy, sharey=ax)

        # no labels
        ax_histy.tick_params(axis="y", labelleft=False)
        ax_histy.set_xlabel("propability")

        # the scatter plot:
        ax.set_xlabel(r"$time (ns)$")
        ax.set_ylabel(ylabel)
        # ax.set_ylim(range)
        ax.plot(x, y, linewidth=0.5, color=color)

        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
        lim = (int(xymax / binwidth) + 1) * binwidth

        # bins = np.arange(-lim, lim + binwidth, binwidth)
        ax_histy.hist(y, bins=20, range=range, orientation='horizontal', color=color, edgecolor='w', linewidth=0.5,
                      density=False, weights=np.ones(len(y)) / len(y), stacked=False, zorder=1)
        plt.savefig(f"{export_file}.png", dpi=300, bbox_inches="tight")
        plt.show()


def vis_overview(md_e, acv_e, exp_e, acv_means, kappa2, rda_acv, rda_md, export_file_name, yp=0.3):
    md_fit_x = burst_fit(md_e.FRETefficiencies, 30)
    acv_fit_x = burst_fit(acv_e.FRETefficiencies, 30)
    exp_fit_x = burst_fit(exp_e, 30)
    with sns.axes_style('ticks'):
        set_ticksStyle()
        f, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 3), sharex=False, sharey=False, squeeze=False)

        print(f"⟨E_MD⟩= {round(np.mean(md_e.FRETefficiencies), 2)} ± {round(np.std(md_e.FRETefficiencies), 2)}")
        print(f"⟨E_MACV⟩= {round(np.mean(acv_e.FRETefficiencies), 2)} ± {round(np.std(acv_e.FRETefficiencies), 2)}")
        print(f"⟨E_Exp⟩= {round(np.mean(exp_e), 2)} ± {round(np.std(exp_e), 2)}")
        print(f"⟨E_MACV_no_sn⟩= {round(np.mean(acv_means), 2)} ± {round(np.std(acv_means), 2)}")
        print(f"⟨E_ACV_start⟩= {round(acv_means[0], 2)}")
        print(f"⟨E_RDA_MD⟩= {round(np.mean(rda_md), 2)} ± {round(np.std(rda_md), 2)}")
        print(f"⟨E_RDA_MACV⟩= {round(np.mean(rda_acv), 2)} ± {round(np.std(rda_acv), 2)}")
        # Plot Experiment und Verteilung Simulation

        ax[0, 0].hist(kappa2, bins=30, range=[0, 4], alpha=0.5, color=c_c(182, 152, 102), edgecolor='w', linewidth=0.5,
                      zorder=1, weights=np.ones(len(kappa2)) / len(kappa2))
        ax[0, 0].axvline(x=kappa2.mean(), ymin=0.9, ymax=1, color=c_c(72, 92, 120))
        ax[0, 0].text(x=kappa2.mean(), y=yp, s=f"{round(kappa2.mean(), 2)}", color=c_c(72, 92, 120))
        ax[0, 0].set_xlabel('$\kappa^2$')
        ax[0, 0].set_ylabel('probability')

        print(f"⟨E_MD⟩= {round(np.mean(md_e.FRETefficiencies), 2)} ± {round(np.std(md_e.FRETefficiencies), 2)}")
        ax[0, 1].hist(md_e.FRETefficiencies, bins=30, range=[0, 1], color=c_c(72, 92, 120), edgecolor='w',
                      linewidth=0.5, label='$\mathregular{{\it E}_{MD}}$', density=False,
                      weights=np.ones(len(md_e.FRETefficiencies)) / len(md_e.FRETefficiencies), stacked=False, zorder=3)
        ax[0, 1].plot(md_fit_x[0], md_fit_x[1], color=c_c(52, 72, 100), linewidth=2, zorder=4)
        ax[0, 1].axvline(x=np.mean(md_e.FRETefficiencies), color=c_c(72, 92, 120))

        ax[0, 1].hist(acv_e.FRETefficiencies, bins=30, range=[0, 1], color=[0.75, 0.51, 0.38], edgecolor='w',
                      linewidth=0.5, density=False,
                      weights=np.ones(len(acv_e.FRETefficiencies)) / len(acv_e.FRETefficiencies), stacked=False,
                      zorder=1, label='$\mathregular{{\it E}_{MACV}}$')
        ax[0, 1].plot(acv_fit_x[0], acv_fit_x[1], color=c_c(152, 103, 72), linewidth=2, zorder=2)
        ax[0, 1].axvline(x=np.mean(acv_means), color=c_c(152, 103, 72))
        ax[0, 1].text(x=np.mean(acv_means), y=1.05, s="E-MACV ohne sn", horizontalalignment='center',
                      verticalalignment='center', transform=ax[0, 1].transAxes)

        ax[0, 1].hist(exp_e, bins=30, range=[0, 1], alpha=0.7, color=c_c(182, 152, 102), edgecolor='w', linewidth=0.5,
                      label=r"$\mathregular{{\it E}_{Exp}}$", density=False, weights=np.ones(len(exp_e)) / len(exp_e),
                      stacked=False, zorder=5)
        ax[0, 1].plot(exp_fit_x[0], exp_fit_x[1], color=c_c(148, 118, 68), linewidth=2, zorder=6)

        ax[0, 1].set_xlabel('FRET', fontsize='large')
        ax[0, 1].set_ylabel('probability', fontsize='large')
        ax[0, 1].set_ylim(0, 0.4)
        ax[0, 1].legend(frameon=False, handlelength=0.75, fontsize='large')

        ax[0, 2].hist(rda_md, bins=30, range=[20, 120], alpha=0.5, color=c_c(182, 152, 102), edgecolor='w',
                      linewidth=0.5, zorder=1, weights=np.ones(len(rda_md)) / len(rda_md))
        ax[0, 2].hist(rda_acv, bins=30, range=[20, 120], alpha=1, color=c_c(66, 117, 106), edgecolor='w', linewidth=0.5,
                      zorder=0, weights=np.ones(len(rda_acv)) / len(rda_acv))

        ax[0, 2].set_xlabel('distance ($\mathregular{\AA}$)')
        ax[0, 2].set_ylabel('probability')
        ax[0, 2].legend(['$\mathregular{{\it R}_{DA-MD}}$', '$\mathregular{{\it R}_{DA-MACV}}$'], frameon=False,
                        handlelength=1, fontsize='large')

        plt.savefig(f"{export_file_name}.png", dpi=300, bbox_inches="tight")
        plt.show()


# Funktion zur Visualisierung von rda
def vis_ap_md(ap_md: list, export_file: str, range: list):
    with sns.axes_style('ticks'):
        set_ticksStyle()
        f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 4), sharex=False, sharey=False, squeeze=False)
        ax[0, 0].hist(ap_md, bins=30, range=range, alpha=1, color=c_c(52, 72, 100), edgecolor='w', linewidth=0.5,
                      zorder=1, weights=np.ones(len(ap_md)) / len(ap_md))

        ax[0, 0].set_xlabel('distance $R$ ($\mathregular{\AA}$)', fontsize=30)
        ax[0, 0].set_ylabel('probability', fontsize=30)
        ax[0, 0].legend(['$\mathregular{{\it R}_{DA}}$'], frameon=False, handlelength=1, fontsize=20)
        ax[0, 0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax[0, 0].spines['top'].set_visible(False)
        ax[0, 0].spines['right'].set_visible(False)
        ax[0, 0].spines['bottom'].set_visible(True)
        ax[0, 0].spines['left'].set_visible(True)
        # plt.xticks(fontsize=14)
        # plt.yticks(fontsize=14)
        ax[0, 0].locator_params(axis="y", nbins=4)
        ax[0, 0].locator_params(axis="x", nbins=4)
        ax[0, 0].tick_params(axis='both', labelsize=20)
        set_ticksStyle()
        plt.savefig(f"{export_file}.svg", dpi=300, bbox_inches="tight")
        plt.show()


def vis_macv_nsn_md(fret_traj_ACV, ACV_init, export_file: str, range: list):
    with sns.axes_style('ticks'):
        set_ticksStyle()
        f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)
        ax[0, 0].hist(fret_traj_ACV.mean_E_DA, bins=30, range=[0, 1], alpha=1, color=[0.75, 0.51, 0.38], edgecolor='w',
                      linewidth=0.5, zorder=2,
                      weights=np.ones(len(fret_traj_ACV.mean_E_DA)) / len(fret_traj_ACV.mean_E_DA))
        ax[0, 0].axvline(x=ACV_init, ymin=0.8, ymax=0.9, color="#297141", linewidth=3)

        ax[0, 0].set_xlabel('FRET', fontsize=30)
        ax[0, 0].set_ylabel('probability', fontsize=30)
        # ax[0,0].legend(['$\mathregular{{\it E}_{MACV}}$'], frameon=False, handlelength=1, fontsize=20)
        ax[0, 0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax[0, 0].spines['top'].set_visible(False)
        ax[0, 0].spines['right'].set_visible(False)
        ax[0, 0].spines['bottom'].set_visible(True)
        ax[0, 0].spines['left'].set_visible(True)
        ax[0, 0].legend(['$\mathregular{{\it E}_{ACV}}$', '$\mathregular{{\it E}_{MACV}}$'], frameon=False,
                        handlelength=1, fontsize=20)
        # ax[0,0].set_ylim(0,0.25)
        ax[0, 0].locator_params(axis="y", nbins=4)
        ax[0, 0].locator_params(axis="x", nbins=4)
        ax[0, 0].tick_params(axis='both', labelsize=20)
        set_ticksStyle()
        plt.savefig(f"{export_file}.svg", dpi=300, bbox_inches="tight")
        plt.show()


# Funktion zur unterstützung bei der Farbstoffentfernung
def make_ndx_of_rna_without_dyes(gro_file: str, output_file: str):
    """RNA extractor

    This function extracts atom id's belonging to an RNA Molecule and not to dyes and writes an .ndx file for structure extraction in GROMACS

    :param gro_file: path to gro file
    :param output_file: path to the output file where the (should end with .ndx)
    """
    rna_atom_dict = {
        "RU": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N1", "C6",
               "H6", "C5", "H5", "C4", "O4", "N3", "H3", "C2", "O2", "C3'", "H3'", "C2'", "H2'1", "O2'", "HO'2",
               "O3'", ],
        "RG": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N9", "C8",
               "H8", "N7", "C5", "C6", "O6", "N1", "H1", "C2", "N2", "H21", "H22", "N3", "C4", "C3'", "H3'", "C2'",
               "H2'1", "O2'", "HO'2", "O3'"],
        "RA": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N9", "C8",
               "H8", "N7", "C5", "C6", "N6", "H61", "H62", "N1", "C2", "H2", "N3", "C4", "C3'", "H3'", "C2'", "H2'1",
               "O2'", "HO'2", "O3'"],
        "RC": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N1", "C6",
               "H6", "C5", "H5", "C4", "N4", "H41", "H42", "N3", "C2", "O2", "C3'", "H3'", "C2'", "HO'2", "O3'"]
    }
    atoms = []
    with open(gro_file) as f:
        next(f)
        next(f)
        for i, line in enumerate(f):
            if "SOL" not in line:
                atom = line.split()[1]
                base = re.sub(r'\d', '', line.split()[0])
                if any(string in base for string in list(rna_atom_dict.keys())):
                    if atom in rna_atom_dict[base[:2]]:
                        atoms.append(line.split()[2])

    with open(output_file, 'w') as the_file:
        the_file.write('[RNA without Dyes]\n')
        counter = 1
        for element in atoms:
            if counter == 15:
                the_file.write(f"{element.rjust(5)}\n")
                counter = 1
            else:
                the_file.write(f"{element.rjust(5)}")
                counter = counter + 1
        the_file.write('\n')


def run_command_win(command: str, cmd_in: bytes):
    """
    Run a command in bash with user input. Calls python subprocess module

    :param command: bash command
    :param cmd_in: commandline input in bytes
    :return: none
    """
    process = subprocess.Popen(["bash", "-c", command], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    process.communicate(input=cmd_in)
    process.wait()


def run_command(command: str):
    """
    Run a bash command with python subprocess.

    :param command: Bash command as stiring.
    :return: none
    """
    process = subprocess.call(["bash", "-c", command], stdout=subprocess.PIPE)
    # print(process.stdout.readline())


def make_ndx_of_rna(gro_file: str, output_file: str):
    """RNA extractor

    This function extracts atom id's belonging to an RNA Molecule and not to dyes and writes an .ndx file for structure extraction in GROMACS

    :param gro_file: path to gro file
    :param output_file: path to the output file where the (should end with .ndx)
    """
    atoms = []
    with open(gro_file) as f:
        next(f)
        next(f)
        for i, line in enumerate(f):
            if not any(value in line for value in ("SOL", "MG", "K", "CL")):
                atom = line.split()[1]
                base = line.split()[0]
                if any(string in base for string in
                       ["RU", "RG", "RA", "RC", "C3W", "C5W", "RGO", "RUM", "A", "U", "C", "G"]):
                    atoms.append(line.split()[2])

    with open(output_file, 'w') as the_file:
        the_file.write('[RNA]\n')
        counter = 1
        for element in atoms:
            if counter == 15:
                the_file.write(f"{element.rjust(5)}\n")
                counter = 1
            else:
                the_file.write(f"{element.rjust(5)}")
                counter = counter + 1
        the_file.write('\n')


def reduce_sds_file(file):
    dot_bracket = ""
    seq = ""
    with open(f"sds_prediction/{file}") as f:
        next(f)
        for i, line in enumerate(f):
            if i == 1:
                print("1", line)
                dot_bracket = line.split()[0]
            if i == 0:
                print("2", line)
                seq = line

    with open("sds_prediction/dot_bracket.secstruct", "w") as text_file:
        text_file.write(dot_bracket + "\n")
        text_file.write(seq)
    # return [dot_bracket, seq[:-1]]


def uncapitalize_fasta_sequence(sequence_file):
    identifier = ""
    seq = ""
    with open(sequence_file) as f:
        for i, line in enumerate(f):
            if i == 0:
                identifier = line
            if i == 1:
                seq = line.lower()

    with open(sequence_file, "w") as text_file:
        text_file.write(identifier)
        text_file.write(seq)


def read_secondary_structure(sequence_file):
    sec_struct = ""
    seq = ""
    with open(sequence_file) as f:
        for i, line in enumerate(f):
            if i == 0:
                sec_struct = line.strip()
            if i == 1:
                seq = line.strip()
    return [sec_struct, seq]


def predict_2D_structure(sequence_file):
    uncapitalize_fasta_sequence(sequence_file)
    make_result_dir("sds_prediction")
    # subprocess.run(["bash","-c"," mkdir sds_prediction"], capture_output= True)
    run_command(f"RNAfold -i {sequence_file} --noPS > sds_prediction/RNA_fold_output.txt")
    # print(subprocess.run(["bash","-c",f"RNAfold -i {sequence_file} --noPS > sds_prediction/RNA_fold_output.txt"], stdout=subprocess.PIPE, text=True))
    reduce_sds_file("RNA_fold_output.txt")


# Idee: Alle Parameter von allen Tools Sammeln und am Enden eine Datei bauen, die einem ein "Laborbuch" eintrag bastelt
def write_rosetta_parameter(secondary_structure_file, rosetta_parameter):
    make_result_dir("rosetta_results")
    secondary_structure = read_secondary_structure(secondary_structure_file)
    if rosetta_parameter["minimize_rna"]:
        file_content = f"{rosetta_parameter['path_to_rosetta']} -nstruct {rosetta_parameter['nstruct']} -sequence '{secondary_structure[1]}'  -secstruct '{secondary_structure[0]}' -silent silent_out.out -minimize_rna {rosetta_parameter['minimize_rna']} -cycles {rosetta_parameter['cycles']}"
    else:
        file_content = f"{rosetta_parameter['path_to_rosetta']} -nstruct {rosetta_parameter['nstruct']} -sequence '{secondary_structure[1]}'  -secstruct '{secondary_structure[0]}' -silent silent_out.out {rosetta_parameter['minimize_rna']} -cycles {rosetta_parameter['cycles']}"

    with open("rosetta_results/FARFAR2.txt", "w") as text_file:
        text_file.write(file_content)


def predict_3D_structure(secondary_structure_file, parameter):
    write_rosetta_parameter(secondary_structure_file, parameter)
    if platform == "linux" or platform == "linux2":
        # run_command("./scripts/linux/rosetta/submitJobs.sh -i rosetta_results/FARFAR2.txt -d rosetta_results -p 1")
        print(subprocess.run(["bash", "-c",
                              "./scripts/linux/rosetta/submitJobs.sh -i rosetta_results/FARFAR2.txt -d rosetta_results -p 1"],
                             stdout=subprocess.PIPE, text=True))
    elif platform == "darwin":
        # run_command("./scripts/mac_os/rosetta/submitJobs.sh -i rosetta_results/FARFAR2.txt -d rosetta_results -p 1")
        print(subprocess.run(["bash", "-c",
                              "./scripts/mac_os/rosetta/submitJobs.sh -i rosetta_results/FARFAR2.txt -d rosetta_results -p 1"],
                             stdout=subprocess.PIPE, text=True))


def extract_pdb(number_of_pdb):
    if platform == "linux" or platform == "linux2":
        # run_command(f"./scripts/linux/rosetta/extract_pdb.sh -d ./rosetta_results/out/1/ -n {number_of_pdb} -m true -s ./rosetta_results/out/1/silent_out.out")
        print(subprocess.run(["bash", "-c",
                              f"./scripts/linux/rosetta/extract_pdb.sh -d ./rosetta_results/out/1/ -n {number_of_pdb} -m true -s ./rosetta_results/out/1/silent_out.out"],
                             stdout=subprocess.PIPE, text=True))
    elif platform == "darwin":
        # run_command(f"./scripts/mac_os/rosetta/extract_pdb.sh -d ./rosetta_results/out/1/ -n {number_of_pdb} -m true -s ./rosetta_results/out/1/silent_out.out")
        print(subprocess.run(["bash", "-c",
                              f"./scripts/mac_os/rosetta/extract_pdb.sh -d ./rosetta_results/out/1/ -n {number_of_pdb} -m true -s ./rosetta_results/out/1/silent_out.out"],
                             stdout=subprocess.PIPE, text=True))


def change_tamperature_in_nvt(temperature, dir):
    temp_in_K = grad_to_kelvin(temperature)
    content = []
    with open(f"{dir}/mdp/nvt.mdp", 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith(('ref_t', 'gen_temp')):
                # print(re.sub(pattern = "[0-9]+", repl = str(temp_in_K), string=line))
                content.append(re.sub(pattern="[0-9]+", repl=str(temp_in_K), string=line))
            else:
                content.append(line)

    # print(content)
    with open(f"{dir}/mdp/nvt.mdp", 'w') as f:
        for l in content:
            f.write("%s\n" % l)


def change_tamperature_in_npt(temperature, dir):
    temp_in_K = grad_to_kelvin(temperature)
    content = []
    with open(f"{dir}/mdp/npt.mdp", 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith(('ref_t')):
                # print(re.sub(pattern = "[0-9]+", repl = str(temp_in_K), string=line))
                content.append(re.sub(pattern="[0-9]+", repl=str(temp_in_K), string=line))
            else:
                content.append(line)

    # print(content)
    with open(f"{dir}/mdp/npt.mdp", 'w') as f:
        for l in content:
            f.write("%s\n" % l)


def change_sim_time_in_md0(time, dir):
    simulation_steps = sim_time_to_steps(time)
    content = []
    with open(f"{dir}/mdp/md0.mdp", 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith(('nsteps')):
                # print(re.sub(pattern = "[0-9]+", repl = str(simulation_steps), string=line))
                content.append(re.sub(pattern="[0-9]+", repl=str(simulation_steps), string=line))
            else:
                content.append(line)

    print(content)
    with open(f"{dir}/mdp/md0.mdp", 'w') as f:
        for l in content:
            f.write("%s\n" % l)


def sim_time_to_steps(sim_time):
    return int(1000000 * sim_time / 2)


def grad_to_kelvin(grad):
    return grad + 273


def copy_files_to_sim_dir(parameter, md_dir_name):
    working_dir_path = os.getcwd()
    src_folder = working_dir_path + "/scripts/gromacs"
    dst_folder = working_dir_path + f"/{md_dir_name}"

    if os.path.exists(dst_folder) and os.path.isdir(dst_folder):
        shutil.rmtree(dst_folder)

    make_result_dir(md_dir_name)

    shutil.copytree(src_folder + "/amber14sb_OL15.ff", dst_folder + "/amber14sb_OL15.ff")
    shutil.copytree(src_folder + "/mdp", dst_folder + "/mdp")
    shutil.copy(src_folder + "/single_run.sh", dst_folder + "/single_run.sh")


def prepare_simulation(simulation_parameter, md_dir_name):
    # Kopieren der Files von Skripts in den aktuellen Ordner
    copy_files_to_sim_dir(simulation_parameter, md_dir_name)
    # Ändern der Parameter in den mdp files
    change_tamperature_in_nvt(simulation_parameter["temperature[°C]"], md_dir_name)
    change_tamperature_in_npt(simulation_parameter["temperature[°C]"], md_dir_name)
    change_sim_time_in_md0(simulation_parameter["simulation_time[ns]"], md_dir_name)


def copy_input_model(path_to_model_pdb, dir):
    """
    Compieng the modeling result structure to the MD simulation directory. File is renamed to input.pdb

    :param path_to_model_pdb: Path to the modeling result structure
    :param dir: Path to the MD run directory
    :return:
    """
    shutil.copy(os.getcwd() + f"{path_to_model_pdb}", f"{os.getcwd()}/{dir}/input.pdb")


def solvate_molecule(structureFile: str, simulation_parameter: dict, dir: str, working_dir_path: str):
    """
    Running bash commands with python subprocess to solvate MD run with GROMACS. Reference solvate.sh.

    :param structureFile: Name of input structure file
    :param simulation_parameter: simulation parameter dictionary
    :param dir: Path to MD run directory
    :return: none
    """
    os.chdir(working_dir_path + f"/{dir}")
    print(os.getcwd())
    structureName = structureFile[:-4]
    make_result_dir("em")
    run_command_win(
        f"gmx pdb2gmx -f {structureFile} -o em/{structureName}.gro -p {structureName}.top -i em/{structureName}.itp -missing -ignh",
        b"1\n 3\n")

    run_command(f"gmx editconf -f em/{structureName}.gro -o em/{structureName}.gro -bt dodecahedron -d 1")

    run_command(
        f"gmx solvate -cp em/{structureName}.gro -cs tip4p.gro -o em/{structureName}.gro -p {structureName}.top ")

    run_command(
        f"gmx grompp -f mdp/em.mdp -c em/{structureName}.gro -p {structureName}.top -o em/{structureName}.tpr -po em/{structureName}.mdp -maxwarn 2")

    run_command_win(
        f"gmx genion -s em/{structureName}.tpr -o em/{structureName}.gro -p {structureName}.top -nname Cl -pname K -neutral",
        b"3\n")

    run_command(
        f"gmx grompp -f mdp/em.mdp -c em/{structureName}.gro -p {structureName}.top -o em/{structureName}.tpr -po em/{structureName}.mdp -maxwarn 2")

    if (simulation_parameter["c_magnesium_ions[mol/l]"] > 0):
        run_command_win(
            f"gmx genion -s em/{structureName}.tpr -o em/{structureName}.gro -p {structureName}.top -nname Cl -pname MG -pq 2 -conc {simulation_parameter['c_magnesium_ions[mol/l]']}",
            b"4\n")

        run_command(
            f"gmx grompp -f mdp/em.mdp -c em/{structureName}.gro -p {structureName}.top -o em/{structureName}.tpr -po em/{structureName}.mdp -maxwarn 2")

    run_command(
        f"gmx mdrun -v -s em/{structureName}.tpr -c em/{structureName}.gro -o em/{structureName}.trr -e em/{structureName}.edr -g em/{structureName}.log")
    os.chdir(working_dir_path)


def run_simulation_steps(structureFile, dir, working_dir_path):
    """
    Running bash commands with python subprocess to make a single MD run with GROMACS. Reference single_run.sh

    :param structureFile: Name of input structure file
    :param dir: Path of MD run directory
    :return: none
    """
    os.chdir(working_dir_path + f"/{dir}")
    print(os.getcwd())
    structureName = structureFile[:-4]

    make_result_dir("nvt")
    run_command(
        f"gmx grompp -f mdp/nvt.mdp -c em/{structureName}.gro -r em/{structureName}.gro -p {structureName}.top -o nvt/{structureName}.tpr -po nvt/{structureName}.mdp -maxwarn 2")
    run_command(
        f"gmx mdrun -v -s nvt/{structureName}.tpr -c nvt/{structureName}.gro -x nvt/{structureName}.xtc -cpo nvt/{structureName}.cpt -e nvt/{structureName}.edr -g nvt/{structureName}.log")

    make_result_dir("npt")
    run_command(
        f"gmx grompp -f mdp/npt.mdp -c nvt/{structureName}.gro -r nvt/{structureName}.gro -t nvt/{structureName}.cpt -p {structureName}.top -o npt/{structureName}.tpr -po npt/{structureName}.mdp -maxwarn 2")
    run_command(
        f"gmx mdrun -v -s npt/{structureName}.tpr -c npt/{structureName}.gro -x npt/{structureName}.xtc -cpo npt/{structureName}.cpt -e npt/{structureName}.edr -g npt/{structureName}.log")

    make_result_dir("md0")
    run_command(
        f"gmx grompp -f mdp/md0.mdp -c npt/{structureName}.gro -t npt/{structureName}.cpt -p {structureName}.top -o md0/{structureName}.tpr -po md0/{structureName}.mdp  -maxwarn 2")
    run_command(
        f"gmx mdrun -v -s md0/{structureName}.tpr -c md0/{structureName}.gro -x md0/{structureName}.xtc -cpo md0/{structureName}.cpt -e md0/{structureName}.edr -g md0/{structureName}.log")

    os.chdir(working_dir_path)


def reduce_center_xtc(md_dir):
    """
    Reduce the trajectory to the RNA and center it in the simulation Box.

    At first a ndx file of the RNA is created. Here are only atom id's written belonging to RNA molecules. Then two bash commands are called by python subprocess. These two commands using gmx trjconv to convert trajectories and to produce a pdb file of the frst state with the given ndx file.

    :param md_dir: Path of the MD run directory
    :return: none
    """
    make_ndx_of_rna(f"{md_dir}/md0/input.gro", f"{md_dir}/analysis/Index_Files/RNA.ndx")
    run_command(
        f"gmx trjconv -f {md_dir}/md0/input.xtc -s {md_dir}/md0/input.tpr -o {md_dir}/md0/input_centered.xtc -n {md_dir}_analysis/Index_Files/RNA.ndx -pbc mol -center")
    run_command(
        f"gmx trjconv -f {md_dir}/md0/input.xtc -s {md_dir}/md0/input.tpr -o {md_dir}/md0/input_s1.pdb -n {md_dir}_analysis/Index_Files/RNA.ndx -pbc mol -center -b 1 -e 10")


def remove_dyes_from_trajectory(md_analyse_dir):
    """
    Create a trjectory of MD run without dyes.

    This method creates a ndx file where all atom id's of the gro file are listes except of the dyes or linker atoms. With this ndx file the gromacs trjconv tools produces an xtc file without the dyes. A pdb file of the first state without the dyes ist also produced.

    :param md_analyse_dir: name of the MD analysis directory
    :return: none
    """
    make_ndx_of_rna_without_dyes(f"{md_analyse_dir}/raw/input.gro",
                                 f'{md_analyse_dir}/Index_Files/RNA_without_Dyes_python.ndx')
    run_command(
        f"gmx trjconv -f {md_analyse_dir}/raw/input_centered.xtc -s {md_analyse_dir}/raw/input.tpr -o {md_analyse_dir}/raw/input_unlabeled.xtc -n {md_analyse_dir}/Index_Files/RNA_without_Dyes_python.ndx -pbc mol -center")
    run_command(
        f"gmx trjconv -f {md_analyse_dir}/raw/input_centered.xtc -s {md_analyse_dir}/raw/input.tpr -o {md_analyse_dir}/raw/input_unlabeled_s1.pdb -n {md_analyse_dir}/Index_Files/RNA_without_Dyes_python.ndx -pbc mol -center -b 1 -e 10")


def copy_files_to_raw(md_dir):
    """Prepare the directory for MD analsyis procedure"""
    working_dir_path = os.getcwd()
    src_folder = working_dir_path + f"/{md_dir}/md0"
    dst_folder = working_dir_path + f"/{md_dir}/analysis/raw"
    shutil.copy(src_folder + "/input_centered.xtc", dst_folder + "/input_centered.xtc")
    shutil.copy(src_folder + "/input.gro", dst_folder + "/input.gro")
    shutil.copy(src_folder + "/input.tpr", dst_folder + "/input.tpr")
    shutil.copy(src_folder + "/input_s1.pdb", dst_folder + "/input_s1.pdb")


def make_data_analyis_results_dirs(md_dir: str):
    """
    Function to prepare a directory for the MD data analysis.

    A new folder MD_analysis will be created. The xtc file of the MD run is reduced to RNA and a pdb file of the first state from the trajectory ist created. The pdb, xtc, gro and tpr file will be copied to the analysis directory.

    :param md_dir: Name of the MD run folder which should be analyzed
    :return: str --> Path of analysis directory
    """
    analysis_dir = f"{md_dir}/analysis"
    make_result_dir(analysis_dir)
    make_result_dir(f"{analysis_dir}/raw")
    make_result_dir(f"{analysis_dir}/Images")
    make_result_dir(f"{analysis_dir}/Index_Files")
    reduce_center_xtc(md_dir)
    copy_files_to_raw(md_dir)
    return analysis_dir


def get_atom_ids(gro_file, filter):
    """
    Return atom id's from a given residue number and a residue name.

    :param gro_file: gro file with structure informatuion
    :param res_numbers: Number of residue (10)
    :param res_names: Name of residue (C5W or RU)
    :return: Dataframe of resuidue id's with other informations.
    """
    df_list = []
    for set in filter:
        with open(gro_file, 'r') as f:
            next(f)
            next(f)
            for i, line in enumerate(f):
                if not any(value in line for value in ("SOL")):
                    l = line.split()
                    if l[0].startswith(set[0]):
                        if (l[1] == set[1]):
                            df_list.append([l[0], l[1], l[2]])
    return pd.DataFrame(df_list, columns=["Num/Res", "Atom name", "ID"])


def reduce_gro_file(gro_file, ndx_file):
    """
    Reduce a given gro file to atoms of the RNA. This reduced file is necessary for MD-Analysis.

    Reads the ndx file of RNA and saves the atom id's in a list. Then reads the gro file and filter atom ids from the list. The filtered lines are saved into a list and then written into a new reduced .gro file.

    :param gro: path to gro file
    :param ndx: path to ndx file of RNA
    :return: none
    """
    id_list = []
    with open(ndx_file, 'r') as ndx:
        next(ndx)
        for i, line in enumerate(ndx):
            l = line.strip()
            for id in l.split():
                id_list.append(id)
    # print(id_list)

    file_content = []
    iterator = 0
    file_content.append("Text\n")
    file_content.append(f"{len(id_list)}\n")
    with open(gro_file, 'r') as gro:
        next(gro)
        next(gro)
        for i, line in enumerate(gro):
            if (iterator < len(id_list)):
                # if not any(value in line for value in ("SOL")):
                # line = line.strip()
                l = line.split()

                if l[2] in id_list:
                    # print(line.split())
                    file_content.append(line)
                    iterator = iterator + 1

    with open(f'{gro_file[:-4]}_reduced.gro', 'w') as f:
        for line in file_content:
            f.write(f"{line}")


def line_prepender(filename: str, line: str):
    """
    Adds a line at the beginning of a file. Used to simulate xvg files.

    Src: https://stackoverflow.com/questions/5914627/prepend-line-to-beginning-of-a-file
    :param filename: Path to file
    :param line: Line to put at front of file
    :return: none
    """
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


def write_coordinate_file(file_name: str, dipole):
    """
    Creates an file of coordinates for dipoles of dyes. Needed when anisotropy calculcations are performed.

    :param file_name: Path and name of file, where the coordinates should written in.
    :param dipole: Coordinates of dipole: [[[x,y,z]][[x,y,z]]]
    :return: none
    """
    time = np.arange(0, len(dipole[0]) + 9, 10, dtype=int)
    df = pd.DataFrame(list(
        zip(time, dipole[0][:, 0], dipole[0][:, 1], dipole[0][:, 2], dipole[1][:, 0], dipole[1][:, 1],
            dipole[1][:, 1])))
    print(len(df[0]))
    df.to_csv(file_name, sep='\t', header=False, index=False)

    for i in range(0, 13):
        line_prepender(file_name, "#")


def get_rkappa_file_from_dyes(dir, don_dipole, acc_dipole, mean_don_acc):
    """
    Creates an distance and 𝜅^2 file of a trajectory with explicit dyes.

    The 𝜅^2 values are calculated by the coordinates of the donor and acceptor dipoles. Then the mean inter dye distance ist calculated by the coordinates of mean donor and acceptor atoms. Then the tiesteps 𝜅^2 and mean dye distances are combined to a dataframe and written to a file in the fluorburst folder.

    :param dir: Path to MD-analysis directory
    :param don_dipole: Coordinates of donor dipole: [[[x,y,z]][[x,y,z]]]
    :param acc_dipole: Coordinates of acceptor dipole: [[[x,y,z]][[x,y,z]]]
    :param mean_don_acc: Coordinates of donor dipole: [[[x,y,z]][[x,y,z]]]
    :return: none
    """
    kappa_2 = calculate_kappa_2(don_dipole, acc_dipole)
    time = np.arange(0, len(don_dipole[0]) + 9, 10, dtype=int)
    # print(len(time))
    rda_MD = calculate_inter_dye_distance(mean_don_acc[0], mean_don_acc[1])
    # Datei generieren
    df = pd.DataFrame(list(zip(time, rda_MD, kappa_2)))
    print(len(df[0]))
    df.to_csv(f'{dir}/fluorburst/rkappa_neu.dat', sep='\t', header=False, index=False)
    return df
