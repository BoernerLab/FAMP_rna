import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import lmfit
import json
import fretraj as ft
import MDAnalysis as mda
from numpy import linalg as LA
import famp_data_analysis
from matplotlib import rcParams


class Visualizations:
    def __init__(self, vis_paths, bins=20):
        self._bins = bins
        self.working_dir = vis_paths["working_dir"]
        self.path_dye_simulation = vis_paths["dye_simulation"]
        self.path_macv_simulation = vis_paths["macv_simulation"]
        self.path_experimental_fret = vis_paths["experimental_FRET"]
        self.path_experimental_anisotropy = vis_paths["experimental_anisotropy"]
        self.path_burst_parameter = vis_paths["burst_parameter"]
        self.path_analysis_parameter = vis_paths["analysis_parameter"]
        self._burst_parameter = self.parse_burst_parameter()
        self._analysis_parameter = self.parse_analysis_parameter()

    @property
    def bins(self):
        return self._bins

    @bins.setter
    def bins(self, new_bins):
        self._bins = new_bins

    @property
    def burst_parameter(self):
        return self._burst_parameter

    def parse_burst_parameter(self):
        if self.path_burst_parameter is not None:
            with open(self.path_burst_parameter) as file:
                burst_parameter = json.load(file)
                return burst_parameter
        else:
            return None

    @burst_parameter.setter
    def burst_parameter(self, new_burst_parameter: dict):
        self._burst_parameter = new_burst_parameter

    @property
    def analysis_parameter(self):
        return self._analysis_parameter

    def parse_analysis_parameter(self):
        if self.path_burst_parameter is not None:
            with open(self.path_analysis_parameter) as file:
                analysis_parameter = json.load(file)
                return analysis_parameter
        else:
            return None

    @analysis_parameter.setter
    def analysis_parameter(self, new_analysis_parameter: dict):
        self.analysis_parameter = new_analysis_parameter

    @staticmethod
    def c_c(r, g, b):
        """
        Function to convert rgb color values from e.g. Colobrewer or Paleton in a matplotl
        :param r:
        :param g:
        :param b:
        :return:
        """
        return [round(r / 255, 2), round(g / 255, 2), round(b / 255, 2)]

    @staticmethod
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

    @staticmethod
    def burst_fit(e_fret, bins):
        hist, bins = np.histogram(e_fret, bins=bins, range=(0, 1), weights=np.ones(len(e_fret)) / len(e_fret))
        bincenters = binMid = (bins[1:] + bins[:-1]) / 2

        mod = lmfit.models.GaussianModel()
        pars = mod.guess(hist, x=bincenters)
        out = mod.fit(hist, pars, x=bincenters)
        x_fit = np.linspace(0, 1, 100)
        y_fit = mod.func(x_fit, *list(out.params.valuesdict().values())[:3])

        return x_fit, y_fit

    @staticmethod
    def calculate_kappa_2(don_dipol: list, acc_dipol: list) -> np.ndarray:
        """
        Funktion zur Berechnung von kappa^2 aus den Koordinaten des Dipolvektors, mit zwei gegebenen Atomen

        :param don_dipol: List --> [[[x, v, z],..., [x, v, z]], [[x, v, z],...,[x, v, z]]] Array mit den Koordinaten der Trajektorie die den Dipolvektor des Donorfarbstoffs definieren
        :param acc_dipol: List --> [[[x, v, z],..., [x, v, z]], [[x, v, z],...,[x, v, z]]] Array mit den Koordinaten der Trajektorie die den Dipolvektor des Acceptorfarbstoffs definieren
        :return: Numpy array mit kappa^2 Werten
        """
        # Initializing vectors with zero
        d_vect = np.zeros([len(don_dipol[1]), 3])
        a_vect = np.zeros([len(acc_dipol[1]), 3])
        d_pos = np.zeros([len(don_dipol[1]), 3])
        a_pos = np.zeros([len(acc_dipol[1]), 3])
        # Vektoren
        for i in range(0, int(len(don_dipol) / 2)):
            # Vector von einem Atom zum Anderen = Richtungsvektor
            d_vect = d_vect - don_dipol[2 * i] + don_dipol[2 * i + 1]
            d_pos = d_pos + don_dipol[2 * i] + don_dipol[2 * i + 1]

        for i in range(0, int(len(acc_dipol) / 2)):
            # Vector von einem Atom zum Anderen = Richtungsvektor
            a_vect = a_vect - acc_dipol[2 * i] + acc_dipol[2 * i + 1]
            a_pos = a_pos + acc_dipol[2 * i] + acc_dipol[2 * i + 1]
            # Richtungsvektor bestimmen

        # Vektoren Normieren
        # Euklidische Normierung
        d_vect = np.divide(d_vect, np.expand_dims(LA.norm(d_vect, axis=1), axis=1))
        a_vect = np.divide(a_vect, np.expand_dims(LA.norm(a_vect, axis=1), axis=1))

        d_pos = 1 / len(don_dipol) * d_pos
        a_pos = 1 / len(acc_dipol) * a_pos

        # Vektor zwischen den Mittelpunkten der Farbstoffe
        dist = d_pos - a_pos
        dist_norm = np.divide(dist, np.expand_dims(LA.norm(dist, axis=1), axis=1))

        a = np.sum(d_vect * a_vect, axis=1)
        b = np.sum(d_vect * dist_norm, axis=1)
        c = np.sum(dist_norm * a_vect, axis=1)

        kappa = a - 3 * b * c
        kappa = np.around(kappa ** 2, 7)
        return kappa

    @staticmethod
    def calculate_inter_dye_distance(mean_donor_atom: list, mean_acceptor_atom: list) -> list:
        """
        Function to calculate the distance between to atoms.

        :param mean_donor_atom: List → Trajektorie der Koordinaten des mittleren C Atoms des Donorfarbstoffes
        :param mean_acceptor_atom: List → Trajektorie der Koordinaten des mittleren C Atoms des Acceptorfarbstoffes
        :return: Liste der Distanzen in Angstrom
        """
        return np.round(np.sqrt(np.sum((np.subtract(mean_donor_atom, mean_acceptor_atom)) ** 2, axis=1)), 7)

    def calculate_burst_dye(self):
        if self.burst_parameter["species"]["unix_pattern_don_coords"] is not None and self.burst_parameter["species"]["unix_pattern_acc_coords"] is not None:
            return ft.burst.Experiment(self.path_dye_simulation, self.burst_parameter, compute_anisotropy=True, units='nm')
        elif self.burst_parameter["species"]["unix_pattern_don_coords"] is None and self.burst_parameter["species"]["unix_pattern_acc_coords"] is None:
            return ft.burst.Experiment(self.path_dye_simulation, self.burst_parameter, compute_anisotropy=False, units='nm')
        else:
            print("Anisotropy will not be calculated due to an Error in the Parameter file. "
                  "Check the unix pattern for don and acc cords in the burst parameter file.")
            return ft.burst.Experiment(self.path_dye_simulation, self.burst_parameter, compute_anisotropy=False, units='nm')

    def calculate_burst_macv(self):
        return ft.burst.Experiment(self.path_macv_simulation, self.burst_parameter, compute_anisotropy=False, units='nm')

    def read_exp_fret_dist(self):
        try:
            exp_data = np.genfromtxt(self.path_experimental_fret, delimiter=',', skip_header=True)
            print(f"{exp_data.mean()} ± {exp_data.std()}")
            return exp_data
        except FileNotFoundError:
            print("No experimental data can be found.")

    def set_ticksStyle(self, x_size=30, y_size=30, x_dir='in', y_dir='in'):
        sns.set_style('ticks', {'xtick.major.size': x_size, 'ytick.major.size': y_size, 'xtick.direction': x_dir,
                                'ytick.direction': y_dir})

    def vis_exp_md_acv(self, md_e, acv_e, exp_e, show_md=True, show_acv=True, show_exp=True, export_name="exp_dye_macv"):
        md_fit_x = self.burst_fit(md_e.FRETefficiencies, self.bins)
        acv_fit_x = self.burst_fit(acv_e.FRETefficiencies, self.bins)
        exp_fit_x = self.burst_fit(exp_e, self.bins)
        with sns.axes_style('ticks'):
            self.set_ticksStyle()
            f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)

            # Plot Experiment und Verteilung Simulation
            if show_md:
                print(f"⟨E_MD⟩= {round(np.mean(md_e.FRETefficiencies), 2)} ± {round(np.std(md_e.FRETefficiencies), 2)}")
                ax[0, 0].hist(md_e.FRETefficiencies, bins=self.bins, range=[0, 1], color=self.c_c(72, 92, 120),
                              edgecolor='w', linewidth=0.5, label='$\mathregular{{\it E}_{MD}}$', density=False,
                              weights=np.ones(len(md_e.FRETefficiencies)) / len(md_e.FRETefficiencies), stacked=False,
                              zorder=3)
                ax[0, 0].plot(md_fit_x[0], md_fit_x[1], color=self.c_c(52, 72, 100), linewidth=2, zorder=4)
                # ax[0,0].axvline(x=np.mean(md_e.FRETefficiencies), color=c_c(72, 92, 120))

            if show_acv:
                ax[0, 0].hist(acv_e.FRETefficiencies, bins=self.bins, range=[0, 1],
                              color=[0.75, 0.51, 0.38],
                              edgecolor='w', linewidth=0.5, density=False,
                              weights=np.ones(len(acv_e.FRETefficiencies)) / len(acv_e.FRETefficiencies), stacked=False,
                              zorder=1, label='$\mathregular{{\it E}_{MACV}}$')
                ax[0, 0].plot(acv_fit_x[0], acv_fit_x[1], color=self.c_c(152, 103, 72), linewidth=2, zorder=2)
                # ax[0,0].axvline(x=acv_means, color=c_c(152, 103, 72))
                # ax[0,0].text(x=acv_means, y=1.05, s=r"$⟨\mathregular{{\it E}_{MACV}}⟩$ without shot noise", horizontalalignment='center',verticalalignment='center', transform=ax[0,0].transAxes)

            if show_exp:
                ax[0, 0].hist(exp_e, bins=self.bins, range=[0, 1], alpha=0.7, color=self.c_c(182, 152, 102),
                              edgecolor='w', linewidth=0.5, label=r"$\mathregular{{\it E}_{exp.}}$", density=False,
                              weights=np.ones(len(exp_e)) / len(exp_e), stacked=False, zorder=5)
                ax[0, 0].plot(exp_fit_x[0], exp_fit_x[1], color=self.c_c(148, 118, 68), linewidth=2, zorder=6)

            ax[0, 0].set_xlabel('FRET', fontsize=30)
            ax[0, 0].set_ylabel('probability', fontsize=30)
            ax[0, 0].yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
            ax[0, 0].spines['top'].set_visible(False)
            ax[0, 0].spines['right'].set_visible(False)
            ax[0, 0].spines['bottom'].set_visible(True)
            ax[0, 0].spines['left'].set_visible(True)
            ax[0, 0].locator_params(axis="y", nbins=4)
            ax[0, 0].locator_params(axis="x", nbins=4)
            ax[0, 0].tick_params(axis='both', labelsize=20)
            # plt.xticks(fontsize=14)
            # plt.yticks(fontsize=14)
            ax[0, 0].set_ylim(0, 0.5)
            ax[0, 0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(f"{self.working_dir}/Images/{export_name}.svg", dpi=300, bbox_inches="tight")
            plt.savefig(f"{self.working_dir}/Images/{export_name}.png", dpi=300, bbox_inches="tight")
            plt.show()

    def vis_exp_md_acv_grid(self, md_e, acv_e, exp_e, export_name="exp_dye_macv_grid"):
        md_fit_x = self.burst_fit(md_e.FRETefficiencies, self.bins)
        acv_fit_x = self.burst_fit(acv_e.FRETefficiencies, self.bins)
        exp_fit_x = self.burst_fit(exp_e, self.bins)
        with sns.axes_style('ticks'):
            self.set_ticksStyle()
            f, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 3), sharex=False, sharey=True, squeeze=False)

            # Plot Experiment und Verteilung Simulation
            print(f"⟨E_MD⟩= {round(np.mean(md_e.FRETefficiencies), 2)} ± {round(np.std(md_e.FRETefficiencies), 2)}")
            for i in np.arange(0,3,2):
                ax[0, i].hist(md_e.FRETefficiencies, bins=self.bins, range=[0, 1], color=self.c_c(72, 92, 120),
                              edgecolor='w', linewidth=0.5, label='$\mathregular{{\it E}_{MD}}$', density=False,
                              weights=np.ones(len(md_e.FRETefficiencies)) / len(md_e.FRETefficiencies), stacked=False,
                              zorder=3)
                ax[0, i].plot(md_fit_x[0], md_fit_x[1], color=self.c_c(52, 72, 100), linewidth=2, zorder=4)
                # ax[0,0].axvline(x=np.mean(md_e.FRETefficiencies), color=c_c(72, 92, 120))

            for i in range(1,3):
                ax[0, i].hist(acv_e.FRETefficiencies, bins=self.bins, range=[0, 1],
                              color=[0.75, 0.51, 0.38],
                              edgecolor='w', linewidth=0.5, density=False,
                              weights=np.ones(len(acv_e.FRETefficiencies)) / len(acv_e.FRETefficiencies), stacked=False,
                              zorder=1, label='$\mathregular{{\it E}_{MACV}}$')
                ax[0, i].plot(acv_fit_x[0], acv_fit_x[1], color=self.c_c(152, 103, 72), linewidth=2, zorder=2)
                # ax[0,0].axvline(x=acv_means, color=c_c(152, 103, 72))
                # ax[0,0].text(x=acv_means, y=1.05, s=r"$⟨\mathregular{{\it E}_{MACV}}⟩$ without shot noise", horizontalalignment='center',verticalalignment='center', transform=ax[0,0].transAxes)

            for i in range(0,2):
                ax[0, i].hist(exp_e, bins=self.bins, range=[0, 1], alpha=0.7, color=self.c_c(182, 152, 102),
                              edgecolor='w', linewidth=0.5, label=r"$\mathregular{{\it E}_{exp.}}$", density=False,
                              weights=np.ones(len(exp_e)) / len(exp_e), stacked=False, zorder=5)
                ax[0, i].plot(exp_fit_x[0], exp_fit_x[1], color=self.c_c(148, 118, 68), linewidth=2, zorder=6)

            ax[0, 0].set_xlabel('FRET', fontsize=30)
            ax[0, 0].set_ylabel('probability', fontsize=30)
            ax[0, 0].yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
            ax[0, 0].spines['top'].set_visible(False)
            ax[0, 0].spines['right'].set_visible(False)
            ax[0, 0].spines['bottom'].set_visible(True)
            ax[0, 0].spines['left'].set_visible(True)
            ax[0, 0].locator_params(axis="y", nbins=4)
            ax[0, 0].locator_params(axis="x", nbins=4)
            ax[0, 0].tick_params(axis='both', labelsize=20)
            # plt.xticks(fontsize=14)
            # plt.yticks(fontsize=14)
            #ax[0, 0].set_ylim(0, 0.33)
            #ax[0, 0].legend(loc='center left', bbox_to_anchor=(1, 0.5))

            f.legend(["TEst","TEst","bla"], loc='outside upper right')
            plt.savefig(f"{self.working_dir}/Images/{export_name}.svg", dpi=300, bbox_inches="tight")
            plt.savefig(f"{self.working_dir}/Images/{export_name}.png", dpi=300, bbox_inches="tight")
            plt.show()

    def vis_burst_size_dist(self, burstsizes, export_name="burstzie"):
        with sns.axes_style('ticks'):
            self.set_ticksStyle()
            f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)
            ax[0, 0].hist(burstsizes, bins=20, color="gray")
            #ax[0, 0].set_ylim(0, 200)
            ax[0, 0].set_xlabel('burstsize', fontsize=30)
            ax[0, 0].set_ylabel('counts', fontsize=30)
            ax[0, 0].set_yscale("log")
            ax[0, 0].set_xscale("log")
            ax[0, 0].yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
            ax[0, 0].spines['top'].set_visible(False)
            ax[0, 0].spines['right'].set_visible(False)
            ax[0, 0].spines['bottom'].set_visible(True)
            ax[0, 0].spines['left'].set_visible(True)
            plt.savefig(f"{self.working_dir}/Images/{export_name}.png", dpi=300, bbox_inches="tight")
            plt.show()

    ########################################################################################################################
    # Here old Code #
    ########################################################################################################################

    # Funktion zur Visualisierung von rda
    def vis_rda_acv_md(self, rda_acv: list, rda_md: list, export_file: str, range: list):
        with sns.axes_style('ticks'):
            set_ticksStyle()
            f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)
            ax[0, 0].hist(rda_md, bins=30, range=range, alpha=0.5, color=c_c(52, 72, 100), edgecolor='w', linewidth=0.5,
                          zorder=1, weights=np.ones(len(rda_md)) / len(rda_md))
            ax[0, 0].hist(rda_acv, bins=30, range=range, alpha=1, color=[0.75, 0.51, 0.38], edgecolor='w',
                          linewidth=0.5,
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
    def vis_kappa_md(self, kappa2, export_file: str, yp: float = 0.3):
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
    def vis_plot_hist(self, x, y, ylabel, export_file, color, range):
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

    def vis_overview(self, md_e, acv_e, exp_e, acv_means, kappa2, rda_acv, rda_md, export_file_name, yp=0.3):
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

            ax[0, 0].hist(kappa2, bins=30, range=[0, 4], alpha=0.5, color=c_c(182, 152, 102), edgecolor='w',
                          linewidth=0.5,
                          zorder=1, weights=np.ones(len(kappa2)) / len(kappa2))
            ax[0, 0].axvline(x=kappa2.mean(), ymin=0.9, ymax=1, color=c_c(72, 92, 120))
            ax[0, 0].text(x=kappa2.mean(), y=yp, s=f"{round(kappa2.mean(), 2)}", color=c_c(72, 92, 120))
            ax[0, 0].set_xlabel('$\kappa^2$')
            ax[0, 0].set_ylabel('probability')

            print(f"⟨E_MD⟩= {round(np.mean(md_e.FRETefficiencies), 2)} ± {round(np.std(md_e.FRETefficiencies), 2)}")
            ax[0, 1].hist(md_e.FRETefficiencies, bins=30, range=[0, 1], color=c_c(72, 92, 120), edgecolor='w',
                          linewidth=0.5, label='$\mathregular{{\it E}_{MD}}$', density=False,
                          weights=np.ones(len(md_e.FRETefficiencies)) / len(md_e.FRETefficiencies), stacked=False,
                          zorder=3)
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

            ax[0, 1].hist(exp_e, bins=30, range=[0, 1], alpha=0.7, color=c_c(182, 152, 102), edgecolor='w',
                          linewidth=0.5,
                          label=r"$\mathregular{{\it E}_{Exp}}$", density=False,
                          weights=np.ones(len(exp_e)) / len(exp_e),
                          stacked=False, zorder=5)
            ax[0, 1].plot(exp_fit_x[0], exp_fit_x[1], color=c_c(148, 118, 68), linewidth=2, zorder=6)

            ax[0, 1].set_xlabel('FRET', fontsize='large')
            ax[0, 1].set_ylabel('probability', fontsize='large')
            ax[0, 1].set_ylim(0, 0.4)
            ax[0, 1].legend(frameon=False, handlelength=0.75, fontsize='large')

            ax[0, 2].hist(rda_md, bins=30, range=[20, 120], alpha=0.5, color=c_c(182, 152, 102), edgecolor='w',
                          linewidth=0.5, zorder=1, weights=np.ones(len(rda_md)) / len(rda_md))
            ax[0, 2].hist(rda_acv, bins=30, range=[20, 120], alpha=1, color=c_c(66, 117, 106), edgecolor='w',
                          linewidth=0.5,
                          zorder=0, weights=np.ones(len(rda_acv)) / len(rda_acv))

            ax[0, 2].set_xlabel('distance ($\mathregular{\AA}$)')
            ax[0, 2].set_ylabel('probability')
            ax[0, 2].legend(['$\mathregular{{\it R}_{DA-MD}}$', '$\mathregular{{\it R}_{DA-MACV}}$'], frameon=False,
                            handlelength=1, fontsize='large')

            plt.savefig(f"{export_file_name}.png", dpi=300, bbox_inches="tight")
            plt.show()

    # Funktion zur Visualisierung von rda
    def vis_ap_md(self, ap_md: list, export_file: str, range: list):
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

    def vis_macv_nsn_md(self, fret_traj_ACV, ACV_init, export_file: str, range: list):
        with sns.axes_style('ticks'):
            set_ticksStyle()
            f, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3), sharex=False, sharey=False, squeeze=False)
            ax[0, 0].hist(fret_traj_ACV.mean_E_DA, bins=30, range=[0, 1], alpha=1, color=[0.75, 0.51, 0.38],
                          edgecolor='w',
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


if __name__ == '__main__':

    paths = {
        "working_dir": "/Users/felixerichson/Documents/Simulationen/Simulations_KLTL_complete/BTL_Paper",
        "dye_simulation": "/Users/felixerichson/Documents/Simulationen/Simulations_KLTL_complete/BTL_Paper/fluorburst/BTL_labeled",
        "macv_simulation": "/Users/felixerichson/Documents/Simulationen/Simulations_KLTL_complete/BTL_Paper/MACV/BTL_labeled",
        "experimental_FRET": "/Users/felixerichson/Documents/Simulationen/Simulations_KLTL_complete/BTL_CEM_Dist_Rest_unrestraint_2/analysis/raw/FRET_40_20.csv",
        "experimental_anisotropy": "/Users/felixerichson/Documents/Simulationen/Simulations_KLTL_complete/BTL_CEM_Dist_Rest_unrestraint_2/analysis/raw/data_sample.txt",
        "burst_parameter": "standard_parameter_files/burst_parameter.json",
        "analysis_parameter": "standard_parameter_files/analysis_parameter.json"
    }

    test = Visualizations(paths, bins=15)
    print(test.burst_parameter)
    e_macv = test.calculate_burst_macv()
    e_dye = test.calculate_burst_dye()
    #e_dye_2 = test.calculate_burst_dye()
    #e_dye_3 = test.calculate_burst_dye()
    e_exp = test.read_exp_fret_dist()

    #print(e_dye.decaytimes_AA["A_f"])
    #print(e_dye.bursts)

    #test.vis_exp_md_acv(e_dye,e_macv,e_exp, show_acv= True, show_exp=True, show_md=False,export_name="corrected exp_with_est_burstsize_fret_sim")

    test.vis_burst_size_dist(e_macv.burstsizes, export_name="burstsize_from_exp")
    print(len(np.unique(e_dye.FRETefficiencies)))
    print(len(e_dye.FRETefficiencies))

    #test.vis_exp_md_acv_grid(e_dye,e_macv,e_exp,export_name="unrestraint_dist_QA_03_QD_0074_gamma_false")
    #test.vis_burst_size_dist()



