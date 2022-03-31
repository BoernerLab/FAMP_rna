import sys

import sys
#import fretraj as ft
import numpy as np
import seaborn as sns
import pandas as pd
import mdtraj as md
import lmfit
import burst as ft
from PyQt5.QtWidgets import (
    QApplication, QDialog, QMainWindow, QMessageBox, QPushButton, QVBoxLayout, QFileDialog
)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from PyQt5.uic import loadUi

from Fluorburst_Designed import Ui_Main


class Window(QMainWindow, Ui_Main):

    def __init__(self, parent=None):
        self.parameters = {}
        self.experiment = []
        self.trajectory_list = {}
        self.xtc_filename = ""
        self.pdb_file_name = ""
        super().__init__(parent)
        self.setupUi(self)
        # self.connectSignalsSlots()
        #print(self.i_qa.text())

        self.clicked = ""

        self.b_load_xtc_pdb.clicked.connect(self.load_trajectory)
        #self.b_load_pdb.clicked.connect(self.clicked_pdb)
        self.b_calculate.clicked.connect(self.print_pars)
        self.b_FRET.clicked.connect(self.vis_exp_md_acv)
        self.b_kappa2.clicked.connect(self.vis_kappa_md)
        self.b_rda.clicked.connect(self.vis_rda_md)
        self.r_ap.clicked.connect(self.vis_rap_md)

    def load_trajectory(self):
        xtc_file_name = QFileDialog.getOpenFileName(self, "Load trajectory","","xtc files (*.xtc)")
        if xtc_file_name:
            pdb_file_name = QFileDialog.getOpenFileName(self, "Open Topology File", "", "PDB files (*.pdb)")
            if pdb_file_name:
                self.xtc_filename = xtc_file_name[0]
                self.pdb_file_name = pdb_file_name[0]
                print('check')
                #self.label_xtc_file.setText("  loading ...")
                print('Check2')
                QApplication.processEvents()

                traj = md.load(self.xtc_filename, top=self.pdb_file_name)

                self.trajectory_list[self.xtc_filename] = traj
                traj_df = traj.top.to_dataframe()[0]

                # TODO Die Farbstoffe sollten über die GUI Spezifiziert werden zumindest Farbstoff und die Residue Nummer --> Dipolatome sollten für Farbstoffe vordefiniert werden
                D14 = traj_df.iloc[traj_df[((traj_df['name'] == 'C14') & (
                        traj_df['resName'] == 'C3W'))].index.values[0]]['serial']
                D2 = traj_df.iloc[traj_df[((traj_df['name'] == 'C2') & (
                        traj_df['resName'] == 'C3W'))].index.values[0]]['serial']
                A14 = traj_df.iloc[traj_df[((traj_df['name'] == 'C14') & (
                        traj_df['resName'] == 'C5W'))].index.values[0]]['serial']
                A2 = traj_df.iloc[traj_df[((traj_df['name'] == 'C2') & (
                        traj_df['resName'] == 'C5W'))].index.values[0]]['serial']
                MD = traj_df.iloc[traj_df[((traj_df['name'] == 'C11') & (
                        traj_df['resName'] == 'C3W'))].index.values[0]]['serial']
                MA = traj_df.iloc[traj_df[((traj_df['name'] == 'C32') & (
                        traj_df['resName'] == 'C5W'))].index.values[0]]['serial']
                print('Check3')
                donor = [traj.xyz[:, D14], traj.xyz[:, D2]]
                acceptor = [traj.xyz[:, A14], traj.xyz[:, A2]]
                mean_donor = traj.xyz[:, MD]
                mean_acceptor = traj.xyz[:, MA]
                print('Check4')
                kappa2 = self.calculate_kappa2(donor, acceptor)
                R_DA = self.calculate_inter_dye_distance(mean_donor, mean_acceptor)
                time = np.arange(0, 1000010, 10)

                print(len(kappa2), len(R_DA), len(time))

                df = pd.DataFrame(list(zip(time, R_DA, kappa2)))
                print("Check1")
                df.to_csv('Data/rkappa.dat', sep='\t', header=False, index=False)
                print("Check 2")


    """
    def clicked_pdb(self):
        # Spezifizieren Dateiname
        fname = QFileDialog.getOpenFileName(self, "Open Topology File", "", "PDB files (*.pdb)")
        if fname:
            self.label_pdb_file.setText(fname[0].split("/")[-1])
            self.pdb_file_name = fname[0]
            traj = md.load(self.xtc_filename, top=self.pdb_file_name)
            traj_df = traj.top.to_dataframe()[0]
            # TODO Die Farbstoffe sollten über die GUI Spezifiziert werden zumindest Farbstoff und die Residue Nummer --> Dipolatome sollten für Farbstoffe vordefiniert werden
            D14 = traj_df.iloc[traj_df[((traj_df['name'] == 'C14') & (
                    traj_df['resName'] == 'C3W'))].index.values[0]]['serial']
            D2 = traj_df.iloc[traj_df[((traj_df['name'] == 'C2') & (
                    traj_df['resName'] == 'C3W'))].index.values[0]]['serial']
            A14 = traj_df.iloc[traj_df[((traj_df['name'] == 'C14') & (
                    traj_df['resName'] == 'C5W'))].index.values[0]]['serial']
            A2 = traj_df.iloc[traj_df[((traj_df['name'] == 'C2') & (
                    traj_df['resName'] == 'C5W'))].index.values[0]]['serial']
            MD = traj_df.iloc[traj_df[((traj_df['name'] == 'C11') & (
                    traj_df['resName'] == 'C3W'))].index.values[0]]['serial']
            MA = traj_df.iloc[traj_df[((traj_df['name'] == 'C32') & (
                    traj_df['resName'] == 'C5W'))].index.values[0]]['serial']

            donor = [traj.xyz[:, D14], traj.xyz[:, D2]]
            acceptor = [traj.xyz[:, A14], traj.xyz[:, A2]]
            mean_donor = traj.xyz[:, MD]
            mean_acceptor = traj.xyz[:, MA]

            kappa2 = self.calculate_kappa2(donor, acceptor)
            R_DA = self.calculate_inter_dye_distance(mean_donor, mean_acceptor)
            time = np.arange(0, 1000010, 10)

            print(len(kappa2), len(R_DA), len(time))

            df = pd.DataFrame(list(zip(time, R_DA, kappa2)))
            print("Check1")
            df.to_csv('Data/rkappa.dat', sep='\t', header=False, index=False)
            print("Check 2")
    """
    def set_ticksStyle(self, x_size=4, y_size=4, x_dir='in', y_dir='in'):
        sns.set_style('ticks', {'xtick.major.size': x_size, 'ytick.major.size': y_size, 'xtick.direction': x_dir,
                                'ytick.direction': y_dir})

    def calculate_inter_dye_distance(self, mean_donor_atom: list, mean_acceptor_atom: list) -> list:
        """
        Funktion zur Berechnung der Distanz zweier Atome

        :param mean_donor_atom: List --> Trajektorie der Koordinaten des mittleren C Atoms des Donorfarbstoffes
        :param mean_acceptor_atom: List --> Trajektorie der Koordinaten des mittleren C Atoms des Acceptorfarbstoffes
        :return: Liste der Distanzen in Angstrom
        """
        return np.sqrt(np.sum((np.subtract(mean_donor_atom, mean_acceptor_atom)) ** 2, axis=1))

    def calculate_kappa2(self, don_dipol: list, acc_dipol: list) -> list:
        """
        Funktion zur Berechnung von kappa^2 aus den Koordinaten des Dipolvektors, mit zwei gegebenen Atomen

        :param don_dipol: List --> [[[x, v, z],..., [x, v, z]], [[x, v, z],...,[x, v, z]]] Arry mit den Koordinaten der Trajektorie die den Dipolvektor des Donorfarbstoffs definieren
        :param acc_dipol: List --> [[[x, v, z],..., [x, v, z]], [[x, v, z],...,[x, v, z]]] Arry mit den Koordinaten der Trajektorie die den Dipolvektor des Acceptorfarbstoffs definieren
        :return: Numpy array mit kappa^2 Werten
        """
        kappas = []
        A14 = acc_dipol[0]
        D14 = don_dipol[0]
        A2 = acc_dipol[1]
        D2 = don_dipol[1]

        # Richtungsvektor bestimmen
        dvect = 0 - D2 + D14
        avect = 0 - A2 + A14
        # Vektoren Normieren
        dvect = np.divide(dvect, np.expand_dims(np.linalg.norm(dvect, axis=1), axis=1))
        avect = np.divide(avect, np.expand_dims(np.linalg.norm(avect, axis=1), axis=1))

        # Mittelpunkt des Farbstoffes bestimmen!
        # TODO Zeilen 27 - 31 zusammenführen!
        # Endpostion der Ortsvektoren
        dpos = 0 + D2 + D14
        apos = 0 + A2 + A14
        # Die Nenner beim Bruch ist die Anzahl der Atome?
        # Mittelpunkt des Farbstoffes!
        dpos = 1 / 2 * dpos
        apos = 1 / 2 * apos

        # Vektor zwischen den Mittelpunkten der Farbstoffe
        dist = dpos - apos
        distnorm = np.divide(dist, np.expand_dims(np.linalg.norm(dist, axis=1), axis=1))

        kappa = np.sum(dvect * avect, axis=1) - (
                    3 * (np.sum(dvect * distnorm, axis=1)) * np.sum(distnorm * avect, axis=1))
        kappa = kappa ** 2
        return kappa

    def c_c(self,r, g, b):
        return [round(r / 255, 2), round(g / 255, 2), round(b / 255, 2)]

    def set_ticksStyle(self, x_size=4, y_size=4, x_dir='in', y_dir='in'):
        sns.set_style('ticks', {'xtick.major.size': x_size, 'ytick.major.size': y_size, 'xtick.direction': x_dir,
                                'ytick.direction': y_dir})

    def burst_fit(self, E_FRET, bins_i):
        hist, bins = np.histogram(E_FRET, bins=bins_i, range=(0, 1), weights=np.ones(len(E_FRET)) / len(E_FRET))
        bincenters = binMid = (bins[1:] + bins[:-1]) / 2

        mod = lmfit.models.GaussianModel()
        pars = mod.guess(hist, x=bincenters)
        out = mod.fit(hist, pars, x=bincenters)
        x_fit = np.linspace(0, 1, 100)
        y_fit = mod.func(x_fit, *list(out.params.valuesdict().values())[:3])

        return x_fit, y_fit

    def print_pars(self):
        self.parameters = {
            "dyes": {
                "tauD": float(self.i_taud.text()),
                "tauA": float(self.i_tauA.text()),
                "QD": float(self.i_qd.text()),
                "QA": float(self.i_qa.text()),
                "dipole_angle_abs_em": 10.5
            },
            "sampling": {
                "nbursts": int(self.i_bursts.text()),
                "skipframesatstart": int(self.i_sfstart.text()),
                "skipframesatend": int(self.i_sfend.text()),
                "multiprocessing": True
            },
            "fret": {
                "R0": float(self.I_r0.text()),
                "kappasquare": float(self.i_kappa.text()),
                "no_gamma": self.r_gamma.isChecked(),
                "quenching_radius": float(self.i_quench.text())
            },
            "species": {
                "name": ["all"],
                "unix_pattern_rkappa": ["*.dat"],
                "unix_pattern_don_coords": ["*Cy3*.xvg"],
                "unix_pattern_acc_coords": ["*Cy5*.xvg"],
                "probability": [1],
                "n_trajectory_splits": None
            },
            "bursts": {
                "lower_limit": int(self.i_llimit.text()),
                "upper_limit": int(self.i_ulimit.text()),
                "lambda": float(self.i_lambda.text()),
                "QY_correction": False,
                "averaging": "all",
                "burst_size_file": None
            }
        }
        #print(parameters)

        experiment_fd = ft.Experiment('Data', self.parameters, compute_anisotropy=False)
        self.experiment = experiment_fd.FRETefficiencies

    def change_button_style(self):
        clicked = self.clicked
        print("Huh")
        if clicked == "FRET":
            FRET_bg = "#ffe278"
        else:
            FRET_bg = "#ffffff"

        if clicked == "kappa":
            kappa_bg = "#ffe278"
        else:
            kappa_bg = "#ffffff"

        if clicked == "rda":
            rda_bg = "#ffe278"
        else:
            rda_bg = "#ffffff"

        if clicked == "ap":
            r_ap_bg = "#ffe278"
        else:
            r_ap_bg = "#ffffff"

        self.b_FRET.setStyleSheet(f"QPushButton{{background-color: {FRET_bg}; border-top-left-radius: 6px; border-left: 2px solid #000000; border-right: 2px solid #000000; border-top: 2px solid #000000;}} QPushButton:hover{{ background-color: #ffe278;}}")

        self.b_kappa2.setStyleSheet(f"QPushButton{{ background-color: {kappa_bg}; border-right: 2px solid #000000; border-top: 2px solid #000000;}} QPushButton:hover{{background-color: #ffe278;}}")

        self.b_rda.setStyleSheet(f"QPushButton{{background-color: {rda_bg}; border-right: 2px solid #000000;border-top: 2px solid #000000;}} QPushButton:hover{{background-color: #ffe278;}}")

        self.r_ap.setStyleSheet(f"QPushButton{{ background-color: {r_ap_bg}; border-top-right-radius: 6px;border-right: 2px solid #000000;border-top: 2px solid #000000;}} QPushButton:hover{{background-color: #ffe278;}}")


    def vis_exp_md_acv(self):
        print("öhm")
        self.clicked = "FRET"
        self.change_button_style()
        # md_e = self.experiment_fd
        md_fit_x = self.burst_fit(self.experiment, 30)
        with sns.axes_style('ticks'):
            self.set_ticksStyle()
            self.figure.clear()

            ax = self.figure.add_subplot(111)

            # Plot Experiment und Verteilung Simulation
            ax.hist(self.experiment, bins=30, range=[0, 1], color=self.c_c(72, 92, 120),
                               edgecolor='w', linewidth=0.5, label='Fluorburst', density=False,
                               weights=np.ones(len(self.experiment)) / len(self.experiment),
                               stacked=False,
                               zorder=3)
            ax.plot(md_fit_x[0], md_fit_x[1], color=self.c_c(52, 72, 100), linewidth=2, zorder=4)

            ax.axvline(x=np.mean(self.experiment), color=self.c_c(72, 92, 120))

            ax.set_xlabel('FRET', fontsize='large')
            ax.set_ylabel('probability', fontsize='large')
            ax.set_ylim(0, 0.3)
            ax.legend(['$\mathregular{{\it E}_{MD}}$', '$\mathregular{{\it E}_{MACV}}$',
                                  '$\mathregular{{\it E}_{Experiment}}$'], frameon=False, handlelength=0.75,
                                 fontsize='large')
            self.canvas.draw()

    def vis_kappa_md(self):
        self.clicked = "kappa"
        self.change_button_style()
        read_kappa = pd.read_csv(f'Data/rkappa.dat', delim_whitespace=True, header=None,
                                 names=['step', 'RDA', 'kappa'], index_col=False)
        print(read_kappa['kappa'].mean())
        with sns.axes_style('ticks'):
            self.figure.clear()
            self.set_ticksStyle()
            ax = self.figure.add_subplot(111)
            ax.hist(read_kappa['kappa'], bins=30, range=[0, 4], alpha=0.5,
                          color=self.c_c(182, 152, 102), edgecolor='w', linewidth=0.5, zorder=1,
                          weights=np.ones(len(read_kappa['kappa'])) / len(read_kappa['kappa']))
            ax.axvline(x=read_kappa['kappa'].mean(), ymin=0.9, ymax=1, color=self.c_c(72, 92, 120))
            ax.text(x=(read_kappa['kappa'].mean()+0.1)/4, y=0.85, ha='center', va='center', transform=ax.transAxes, s=f"{round(read_kappa['kappa'].mean(), 2)}",
                         color=self.c_c(72, 92, 120))
            ax.set_xlabel('$\kappa^2$', fontsize="large")
            ax.set_ylabel('probability', fontsize="large")
            self.canvas.draw()

    def vis_rda_md(self):
        self.clicked = "rda"
        self.change_button_style()

    def vis_rap_md(self):
        self.clicked = "ap"
        self.change_button_style()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Window()
    win.show()
    sys.exit(app.exec())
