import subprocess
import os
import MDAnalysis as mda
import numpy as np
import pandas
from numpy import linalg as LA
import pathlib
import re
import fileinput
import shutil
import pandas as pd
import mdtraj as md
import fretraj as ft
from famp.pdb_cleaner import *


class Dye:
    def __init__(self, dye_parameter: tuple) -> None:
        self.source_path = pathlib.Path(__file__).parent.resolve()
        self.ff_abbreviation = dye_parameter[0]
        self.attachment_residue = dye_parameter[1]
        self.attach_chain_id = dye_parameter[2]
        self.dye_name = None
        self.central_c = None
        self.dipole_names = None
        self.dipole_ids = None
        self.central_c_id = None
        self.attechmentpoint_id = None
        self.get_attributes_from_file()

    @staticmethod
    def read_pdb_file(pdb_file):
        column_names = []
        content = []
        if os.path.isfile(pdb_file):
            with open(pdb_file, "r") as file:
                for line in file:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        splitted_line = line.split()
                        if len(splitted_line) == 11:
                            column_names = ["type", "atom_id", "atom_name", "residue_name", "residue_number", "X", "Y",
                                            "Z", "occupancy", "temperature_factor", "element"]
                        elif len(splitted_line) == 12:
                            column_names = ["type", "atom_id", "atom_name", "residue_name", "chain_id",
                                            "residue_number",
                                            "X", "Y", "Z", "occupancy", "temperature_factor", "element"]
                        else:
                            print("The PDB file has more columns than expected.")
                        content.append(splitted_line)
        else:
            print(f"Cant find file: {pdb_file}"
                  f" Please check if this file and folder exist. Rename the file if necessary.")

        pdb_df = pd.DataFrame(content, columns=column_names)
        return pdb_df

    def get_attechment_id_from_pdb(self, pdb_file):
        pdb_df = self.read_pdb_file(pdb_file)
        #print(pdb_file)
        residue_numbers = pdb_df["residue_number"].unique()
        attechment_atom = ""
        if self.attachment_residue == 1:
            atom_id = pdb_df[(pdb_df["residue_number"] == str(1)) & (pdb_df["atom_name"] == "O5'") & (pdb_df["chain_id"] == self.attach_chain_id)]
            attechment_atom = "O5'"
        elif self.attachment_residue == int(residue_numbers[-1]):
            atom_id = pdb_df[(pdb_df["residue_number"] == str(residue_numbers[-1])) & (pdb_df["atom_name"] == "O3'") & (pdb_df["chain_id"] == self.attach_chain_id)]
            attechment_atom = "O3'"
        else:
            atom_id = pdb_df[(pdb_df["residue_number"] == str(self.attachment_residue)) & (pdb_df["atom_name"] == "C5") & (pdb_df["chain_id"] == self.attach_chain_id)]
            attechment_atom = "C5"

        atom_id = int(atom_id["atom_id"])
        print(f"Attechmentpoint for {self.attachment_residue} {self.ff_abbreviation}: {atom_id}, {attechment_atom}, Chain: {self.attach_chain_id} ")
        return atom_id, attechment_atom

    def get_attributes_from_file(self):
        print(self.ff_abbreviation)
        dye_table = pandas.read_csv(f"{self.source_path}/scripts/dye_properties.csv", header=0)
        dye_row = dye_table.loc[dye_table["dye_abbreviation"] == self.ff_abbreviation]
        self.dye_name = dye_row["dye_name"].to_string(index=False)
        self.central_c = dye_row["central_c"].to_string(index=False)
        self.dipole_names = (dye_row["dipol_1"].to_string(index=False), dye_row["dipol_2"].to_string(index=False))

    def get_ids_from_gro(self, gro_file):
        dipole_id_1 = 0
        dipole_id_2 = 0
        if os.path.isfile(gro_file):
            with open(gro_file, 'r') as f:
                next(f)
                next(f)
                for i, line in enumerate(f):
                    if not any(value in line for value in "SOL"):
                        split_line = line.split()
                        if split_line[0].startswith(f"{self.attachment_residue}{self.ff_abbreviation}"):
                            if split_line[1] == self.central_c:
                                self.central_c_id = int(split_line[2])
                            if split_line[1] == self.dipole_names[0]:
                                dipole_id_1 = int(split_line[2])
                            if split_line[1] == self.dipole_names[1]:
                                dipole_id_2 = int(split_line[2])
                            if dipole_id_1 > 0 and dipole_id_2 > 0:
                                self.dipole_ids = (dipole_id_1, dipole_id_2)
        else:
            print(
                f"Cant find file: {gro_file}"
                f" Please check if this file and folder exist. Rename the file if necessary.")


class DataAnalysis:
    def __init__(self, working_dir, path_sim_results: str, analysis_parameter: dict, macv_label_pars: dict) -> None:
        self.source_path = pathlib.Path(__file__).parent.resolve()
        self.working_dir = working_dir
        self.path_sim_results = path_sim_results
        self.analysis_parameter = analysis_parameter
        self.macv_label_pars = macv_label_pars
        self.analysis_dir = f"{self.path_sim_results}/analysis"
        self.input_structure_name = self.analysis_parameter["input_structure_name"]
        self.md_traj = None
        self.fret_macv = None
        self.simulation_type = None
        self.parameter_result_file_checker()
        self.acceptor_dye = Dye(dye_parameter=self.analysis_parameter["Acceptor_residue_name_number"])
        self.donor_dye = Dye(dye_parameter=self.analysis_parameter["Donor_residue_name_number"])

    @staticmethod
    def run_command(command: str):
        process = subprocess.Popen(["bash", "-c", command], stdout=subprocess.PIPE, text=True)
        while process.stdout.readable():
            line = process.stdout.readline()

            if not line:
                break

            print(line.strip())

    @staticmethod
    def make_dir(directory_name):
        """
        Creates a directory

        :param: directory_name - Name of the new directory
        :return: None
        """
        result_dir = directory_name
        try:
            os.mkdir(result_dir)
        except FileExistsError:
            print(f"Results can be found in: {result_dir}")
        except OSError:
            print(f"Failed to create {result_dir}")
        else:
            print(f"Successfully created the directory {result_dir}. Results can be found there")

    @staticmethod
    def make_ndx_of_rna(gro_file: str, output_file: str):
        """RNA extractor

        This function extracts atom id's belonging to an RNA Molecule and not to dyes and writes an .ndx file for
        structure extraction in GROMACS

        :param gro_file: path to gro file
        :param output_file: path to the output file where the (should end with .ndx)
        """

        atoms = []
        rna_residues = [
            "RU", "RG", "RA", "RC", "C3W", "C5W", "C5K",
            "RGO", "RUM", "A", "U", "C", "G", "RG5", "RA5", "RC5", "RU5", "RG3", "RA3", "RC3", "RU3",
            "DT", "DA", "DG", "DC", "DGP", "DC3", "DGK","DG5", "C3N", "DTM", "C3B", "DCI", "DT3", "C5I", "C3P"
        ]
        exclude_exact = ["SOL", "MG", "K", "CL"]

        with open(gro_file) as f:
            next(f)  # Überspringe Titel
            next(f)  # Überspringe Atomanzahl
            for line in f:
                if len(line.strip()) < 20:
                    continue

                residue_name = line[5:10].strip()
                atom_number = line[15:20].strip()

                if residue_name in exclude_exact:
                    continue

                if residue_name in rna_residues:
                    atoms.append(atom_number)



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

    @staticmethod
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

    @staticmethod
    def make_ndx_of_rna_without_dyes(gro_file: str, output_file: str):
        # TODO: Extend for 3 and 5 prime residues
        """RNA extractor

        This function extracts atom id's belonging to an RNA Molecule and not to dyes and writes an .ndx file
        for structure extraction in GROMACS

        :param gro_file: path to gro file
        :param output_file: path to the output file where the (should end with .ndx)
        """
        rna_atom_dict = {
            "RU": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N1",
                   "C6",
                   "H6", "C5", "H5", "C4", "O4", "N3", "H3", "C2", "O2", "C3'", "H3'", "C2'", "H2'1", "O2'", "HO'2",
                   "O3'", ],
            "RG": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N9",
                   "C8",
                   "H8", "N7", "C5", "C6", "O6", "N1", "H1", "C2", "N2", "H21", "H22", "N3", "C4", "C3'", "H3'", "C2'",
                   "H2'1", "O2'", "HO'2", "O3'"],
            "RA": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N9",
                   "C8",
                   "H8", "N7", "C5", "C6", "N6", "H61", "H62", "N1", "C2", "H2", "N3", "C4", "C3'", "H3'", "C2'",
                   "H2'1",
                   "O2'", "HO'2", "O3'"],
            "RC": ["P", "O1P", "O2P", "O5'", "H5T", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'", "C1'", "H1'", "N1",
                   "C6",
                   "H6", "C5", "H5", "C4", "N4", "H41", "H42", "N3", "C2", "O2", "C3'", "H3'", "C2'", "HO'2", "O3'", "O2'"]
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

    def print_dye_informations(self, residue_numbers, atom_names):
        if self.md_traj is not None:
            traj_df = self.md_traj.top.to_dataframe()[0]
            return traj_df.loc[((traj_df['resSeq'] == residue_numbers[0]) & (traj_df['name'] == atom_names[0])) | (
                    (traj_df['resSeq'] == residue_numbers[1]) & (traj_df['name'] == [1]))]
        else:
            self.set_md_traj()
            traj_df = self.md_traj.top.to_dataframe()[0]
            return traj_df.loc[((traj_df['resSeq'] == residue_numbers[0]) & (traj_df['name'] == atom_names[0])) | (
                    (traj_df['resSeq'] == residue_numbers[1]) & (traj_df['name'] == [1]))]

    def set_md_traj(self):
        traj = md.load(f'{self.analysis_dir}/raw/{self.input_structure_name}_unlabeled.xtc',
                       top=f'{self.analysis_dir}/raw/{self.input_structure_name}_unlabeled_s1.pdb')
        self.md_traj = traj

    # ---------------------Preparing Files for MD Analysis procedure----------------------------------------------------

    def parameter_result_file_checker(self):
        dye_list = ["C3P", "C5I", "C3W", "C5W", "C7N", "C55", "C75", "A35", "A48", "A53", "A56", "A59", "A64", "T39", "T42", "T46",
                    "T48", "T49", "T51", "T52", "T61"]

        par_acceptor = self.analysis_parameter['Acceptor_residue_name_number']
        par_donor = self.analysis_parameter['Donor_residue_name_number']
        if os.path.isdir(self.working_dir):
            pass
        else:
            print(f"The workring directory {self.working_dir} does not exist")

        if os.path.isdir(self.path_sim_results):
            pass
        else:
            print(f"The directory of MD results {self.working_dir} does not exist")

        if os.path.isdir(f"{self.path_sim_results}/md0"):
            results_file_path = f"{self.path_sim_results}/md0/{self.input_structure_name}"
            if os.path.isfile(f"{results_file_path}.gro") and os.path.isfile(f"{results_file_path}.xtc"):
                if par_acceptor[0] in dye_list and par_donor[0] in dye_list:
                    print(f"{results_file_path}.gro")
                    with open(f"{results_file_path}.gro") as file:
                        if f"{par_acceptor[1]}{par_acceptor[0]}" and f"{par_donor[1]}{par_donor[0]}" in file.read():
                            print("Dyes found: Parameter are ok. Dye and MACV calculations can be performed")
                            self.simulation_type = "DYE"
                        else:
                            print("Dyes not found. Please check your parameter input if you have performed a dye "
                                  "simulation. Otherwise only MACV calculation can be performed")
                            self.simulation_type = "MACV"
                else:
                    print("Specified dyes are not known. Please check your entries of the dye names in the "
                          "analyses parameter.")
            else:
                print(f"Cannot find the file {results_file_path}.gro or {results_file_path}.xtc."
                      f"Please make sure that the paths and inputsrtructure name are correct.")
        else:
            print(f"The directory {self.path_sim_results}/md0 expected but not found")

    def reduce_center_xtc(self, pbc_method="mol"):
        """
        Reduce the trajectory to the RNA and center it in the simulation Box.

        At first a ndx file of the RNA is created. Here are only atom id's written belonging to RNA molecules.
        Then two bash commands are called by python subprocess. These two commands using gmx trjconv to convert
        trajectories and to produce a pdb file of the frst state with the given ndx file.

        :return: none
        """
        sim_name = self.analysis_parameter['input_structure_name']
        md_dir = self.path_sim_results
        self.make_ndx_of_rna(f"{md_dir}/md0/{sim_name}.gro", f"{self.analysis_dir}/Index_Files/RNA.ndx")
        self.run_command(
            f"gmx trjconv -f {md_dir}/md0/{sim_name}.xtc -s {md_dir}/md0/{sim_name}.tpr -o {md_dir}/md0/{sim_name}_centered.xtc -n {self.analysis_dir}/Index_Files/RNA.ndx -pbc {pbc_method} -center")
        self.run_command(
            f"gmx trjconv -f {md_dir}/md0/{sim_name}.xtc -s {md_dir}/md0/{sim_name}.tpr -o {md_dir}/md0/{sim_name}_s1.pdb -n {self.analysis_dir}/Index_Files/RNA.ndx -pbc mol -center -b 1 -e 10")
        clean_pdb(f"{md_dir}/md0/{sim_name}_s1.pdb", chain_id=True)

    def export_pdb_trajectory(self, time_steps, unlabeled_traj=False):
        """
        Exports a pdb trajectory from a GROMACS simulation file.

        :params: time_steps: step width for exporting states.
        :return: none
        """
        sim_name = self.analysis_parameter['input_structure_name']
        if unlabeled_traj:
            self.run_command(
                f"gmx trjconv -f {self.analysis_dir}/raw/{sim_name}_centered.xtc -s {self.analysis_dir}/raw/{sim_name}.tpr -o {self.analysis_dir}/raw/{sim_name}_unlabeled_traj.pdb -n {self.analysis_dir}/Index_Files/RNA_without_Dyes_python.ndx -pbc mol -dt {time_steps} -center")

        else:
            self.run_command(
                f"gmx trjconv -f {self.analysis_dir}/raw/{sim_name}_centered.xtc -s {self.analysis_dir}/raw/{sim_name}.tpr -o {self.analysis_dir}/raw/{sim_name}_traj.pdb -n {self.analysis_dir}/Index_Files/RNA.ndx -pbc mol -dt {time_steps} -center")

    def export_range_pdb_trajectory(self, time_steps: int, traj_range: list):
        """
        Exports a pdb trajectory from a GROMACS simulation file.

        :params: time_steps: step width for exporting states.
        :return: none
        """
        sim_name = self.analysis_parameter['input_structure_name']
        self.run_command(
            f"gmx trjconv -f {self.analysis_dir}/raw/{sim_name}_centered.xtc -s {self.analysis_dir}/raw/{sim_name}.tpr -o {self.analysis_dir}/raw/{sim_name}_traj_{traj_range[0]}_{traj_range[1]}.pdb -n {self.analysis_dir}/Index_Files/RNA.ndx -pbc mol -dt {time_steps} -b {traj_range[0]} -e {traj_range[1]} -center")

    def copy_files_to_raw(self):
        """Prepare the directory for MD analysis procedure"""
        src_folder = f"{self.path_sim_results}/md0"
        dst_folder = f"{self.analysis_dir}/raw"
        sim_name = self.analysis_parameter["input_structure_name"]
        shutil.copy(src_folder + f"/{sim_name}_centered.xtc", dst_folder + f"/{sim_name}_centered.xtc")
        shutil.copy(src_folder + f"/{sim_name}.gro", dst_folder + f"/{sim_name}.gro")
        shutil.copy(src_folder + f"/{sim_name}.tpr", dst_folder + f"/{sim_name}.tpr")
        shutil.copy(src_folder + f"/{sim_name}_s1.pdb", dst_folder + f"/{sim_name}_s1.pdb")

    def reduce_gro_file(self):
        """
        Reduce a given gro file to atoms of the RNA. This reduced file is necessary for MD-Analysis.

        Reads the ndx file of RNA and saves the atom id's in a list. Then reads the gro file and filter atom ids from
        the list. The filtered lines are saved into a list and then written into a new reduced .gro file.

        :return: none
        """
        path_to_gro_file = f"{self.analysis_dir}/raw/{self.input_structure_name}.gro"
        path_to_ndx_file = f"{self.analysis_dir}/Index_Files/RNA.ndx"
        if path_to_gro_file and path_to_ndx_file:
            id_list = []
            with open(path_to_ndx_file, 'r') as ndx:
                next(ndx)
                for i, line in enumerate(ndx):
                    l = line.strip()
                    for atom_id in l.split():
                        id_list.append(atom_id)
            # print(id_list)

            file_content = []
            iterator = 0
            file_content.append("Text\n")
            file_content.append(f"{len(id_list)}\n")
            with open(path_to_gro_file, 'r') as gro:
                next(gro)
                next(gro)
                for i, line in enumerate(gro):
                    if iterator < len(id_list):
                        # if not any(value in line for value in ("SOL")):
                        # line = line.strip()
                        l = line.split()

                        if l[2] in id_list:
                            # print(line.split())
                            file_content.append(line)
                            iterator = iterator + 1

            with open(f'{path_to_gro_file[:-4]}_reduced.gro', 'w') as f:
                for line in file_content:
                    f.write(f"{line}")
        else:
            print(f"Please check if these both files {path_to_gro_file} and {path_to_ndx_file} exist.")

    def make_data_analysis_results_dirs(self, pbc_method="mol"):
        """
        Function to prepare a directory for the MD data analysis.

        A new folder MD_analysis will be created. The xtc file of the MD run is reduced to RNA and a pdb file of the
        first state from the trajectory ist created. The pdb, xtc, gro and tpr file will be copied to the analysis
        directory.

        """
        analysis_dir = f"{self.analysis_dir}"
        self.make_dir(analysis_dir)
        self.make_dir(f"{analysis_dir}/raw")
        self.make_dir(f"{analysis_dir}/Images")
        self.make_dir(f"{analysis_dir}/Index_Files")
        self.reduce_center_xtc(pbc_method=pbc_method)
        self.copy_files_to_raw()
        self.reduce_gro_file()

    # ----------------------MACV calculations---------------------------------------------------------------------------
    def remove_dyes_from_trajectory(self):
        """
        Create a trajectory of MD run without dyes.

        This method creates a ndx file where all atom id's of the gro file are listed except of the dyes or linker atoms
        . With this ndx file the gromacs trjconv tools produces a xtc file without the dyes. A pdb file of the first
        state without the dyes ist also produced.

        :return: none
        """
        sim_name = self.analysis_parameter["input_structure_name"]
        self.make_ndx_of_rna_without_dyes(f"{self.analysis_dir}/raw/{sim_name}.gro",
                                          f'{self.analysis_dir}/Index_Files/RNA_without_Dyes_python.ndx')
        self.run_command(
            f"gmx trjconv -f {self.analysis_dir}/raw/{sim_name}_centered.xtc -s {self.analysis_dir}/raw/{sim_name}.tpr -o {self.analysis_dir}/raw/{sim_name}_unlabeled.xtc -n {self.analysis_dir}/Index_Files/RNA_without_Dyes_python.ndx -pbc mol -center")
        self.run_command(
            f"gmx trjconv -f {self.analysis_dir}/raw/{sim_name}_centered.xtc -s {self.analysis_dir}/raw/{sim_name}.tpr -o {self.analysis_dir}/raw/{sim_name}_unlabeled_s1.pdb -n {self.analysis_dir}/Index_Files/RNA_without_Dyes_python.ndx -pbc mol -center -b 1 -e 10")
        clean_pdb(f"{self.analysis_dir}/raw/{sim_name}_unlabeled_s1.pdb", chain_id=True)
        self.rewrite_atoms_after_unlabeling()

    def rewrite_atoms_after_unlabeling(self):
        """
        Function to rename dye specific bases to normal bases in the unlabeled pdb file. Example: RUM --> RU
        :return: none
        """
        path = f'{self.analysis_dir}/raw/{self.input_structure_name}_unlabeled_s1.pdb'
        with fileinput.FileInput(path, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("RUM", " RU"), end='')

        with fileinput.FileInput(path, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("RGO", " RG"), end='')

        with fileinput.FileInput(path, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("RCO", " RC"), end='')

        with fileinput.FileInput(path, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("RUO", " RU"), end='')

        with fileinput.FileInput(path, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace("RAO", " RA"), end='')

    def build_acv_parameter(self) -> dict:

        standart_acv_parameter = {
            "use_LabelLib": False,
            "grid_spacing": 1.0,
            "simulation_type": "AV3",
            "state": 1,
            "frame_mdtraj": 0,
            "contour_level_AV": 0,
            "contour_level_CV": 0.7,
            "b_factor": 100,
            "gaussian_resolution": 2,
            "grid_buffer": 2.0,
            "transparent_AV": True,
            "mol_selection": "all"
        }

        donor_analysis_pars = self.analysis_parameter["Donor_residue_name_number"]
        acceptor_analysis_pars = self.analysis_parameter["Acceptor_residue_name_number"]

        path_to_pdb = f"{self.analysis_dir}/raw/{self.analysis_parameter['input_structure_name']}_unlabeled_s1.pdb"
        donor_ap = self.donor_dye.get_attechment_id_from_pdb(path_to_pdb)
        acceptor_ap = self.acceptor_dye.get_attechment_id_from_pdb(path_to_pdb)

        donor_position = f"{self.macv_label_pars['Donor']['name']}-{donor_analysis_pars[1]}-{donor_ap[1]}"
        acceptor_position = f"{self.macv_label_pars['Acceptor']['name']}-{acceptor_analysis_pars[1]}-{acceptor_ap[1]}"

        donor_params = {}
        donor_params = self.macv_label_pars["Donor"]
        donor_params["attach_id"] = donor_ap[0]
        donor_params.update(standart_acv_parameter)

        acceptor_params = {}
        acceptor_params = self.macv_label_pars["Acceptor"]
        acceptor_params["attach_id"] = acceptor_ap[0]
        acceptor_params.update(standart_acv_parameter)

        labels = {"Position": {}, "Distance": {}}
        labels["Position"][donor_position] = donor_params
        labels["Position"][acceptor_position] = acceptor_params
        labels["Distance"].update(self.macv_label_pars["Distance"])

        return labels

    def get_selected_frames(self):
        time = self.md_traj.time[-1]
        step = self.md_traj.timestep
        max_time = round(time / step, 0)
        time_step = 1
        ts = max_time / 10000
        if ts >= 1:
            time_step = round(ts, 0)
        else:
            time_step = 1

        s_frames = [int(max_time + 1), int(time_step)]
        return s_frames

    def calculate_macv(self, macv_parameter, pkl_file_name="macv_calculation", frame_factor=1):
        s_frames = self.get_selected_frames()
        selected_frames = range(0, s_frames[0], s_frames[1]*frame_factor)
        print(s_frames, selected_frames)

        donor_site = list(macv_parameter["Position"].keys())[0]
        acceptor_site = list(macv_parameter["Position"].keys())[1]
        fret_pair = list(macv_parameter["Distance"].keys())[0]

        print(donor_site,acceptor_site,fret_pair)

        fret = ft.cloud.pipeline_frames(self.md_traj, donor_site, acceptor_site, macv_parameter, selected_frames,
                                        fret_pair)
        ft.cloud.save_obj(f'{self.analysis_dir}/macv/{self.input_structure_name}{pkl_file_name}.pkl', fret)
        return fret

    def load_macv(self, pkl_file_name="macv_calculation"):
        """
        Loads a macv calculation from a pkl file
        :param pkl_file_name: name of the pkl
        """
        fret = ft.cloud.load_obj(f'{self.analysis_dir}/macv/{pkl_file_name}.pkl')
        self.fret_macv = fret
        #print(len(fret))
        return fret

    def write_rkappa_file_from_macv(self, file_name, frame_factor=1):
        s_frames = self.get_selected_frames()
        fret_traj = ft.cloud.Trajectory(self.fret_macv, timestep=(s_frames[1]*frame_factor) * int(self.md_traj.timestep),
                                        kappasquare=0.66)
        fret_traj.save_traj(f'{self.analysis_dir}/macv/R_kappa_{file_name}.dat', format='txt', R_kappa_only=True, units='nm',
                            header=False)
        fret_traj.dataframe.head()

    def genarate_rkappa_file_from_macv(self, calculate_macv=True, pkl_file_name="macv_calculation", frame_factor=1):
        self.make_dir(f"{self.analysis_dir}/macv")
        self.remove_dyes_from_trajectory()
        self.rewrite_atoms_after_unlabeling()
        macv_parameter = self.build_acv_parameter()
        self.set_md_traj()

        if calculate_macv:
            self.fret_macv = self.calculate_macv(macv_parameter, frame_factor=frame_factor)
        else:
            self.fret_macv = self.load_macv()

        self.write_rkappa_file_from_macv(pkl_file_name, frame_factor=frame_factor)



    # -------------------------explicit dye simulation handling---------------------------------------------------------

    def get_atoms_coordinates(self, atom_id, universe):
        universe.trajectory[0]
        acc_c14 = universe.select_atoms(f'id {atom_id}')
        acc_coords = []
        for ts in universe.trajectory:
            coords = list(np.squeeze(acc_c14.positions / 10))
            acc_coords.append(coords)
        return np.array(acc_coords)

    def write_coordinate_file(self, file_name: str, dipole):
        """
        Creates a file of coordinates for dipoles of dyes. Needed when anisotropy calculations are performed.

        :param file_name: Path and name of file, where the coordinates should be written in.
        :param dipole: Coordinates of dipole: [[[x,y,z]][[x,y,z]]]
        :return: none
        """
        time = np.arange(0, len(dipole[0]) + 9, 1, dtype=int)
        df = pd.DataFrame(list(
            zip(time * 10, dipole[0][:, 0], dipole[0][:, 1], dipole[0][:, 2], dipole[1][:, 0], dipole[1][:, 1],
                dipole[1][:, 2])))
        print(len(df[0]))
        df.to_csv(f"{self.analysis_dir}/explicit_dyes/{file_name}", sep='\t', header=False, index=False)

        for i in range(0, 13):
            self.line_prepender(f"{self.analysis_dir}/explicit_dyes/{file_name}", "#")

    def write_rkappa_file_from_dyes(self, don_dipole, acc_dipole, mean_don_acc):
        """
        Creates a distance and 𝜅^2 file of a trajectory with explicit dyes.

        The 𝜅^2 values are calculated by the coordinates of the donor and acceptor dipoles. Then the mean inter dye
        distance ist calculated by the coordinates of mean donor and acceptor atoms. Then the time steps 𝜅^2 and mean
        dye distances are combined to a dataframe and written to a file in the explicit_dyes folder.

        :param don_dipole: Coordinates of donor dipole: [[[x,y,z]][[x,y,z]]]
        :param acc_dipole: Coordinates of acceptor dipole: [[[x,y,z]][[x,y,z]]]
        :param mean_don_acc: Coordinates of donor dipole: [[[x,y,z]][[x,y,z]]]
        :return: none
        """
        kappa_2 = self.calculate_kappa_2(don_dipole, acc_dipole)
        time = np.arange(0, len(don_dipole[0]) + 9, 1, dtype=int)
        rda_md = self.calculate_inter_dye_distance(mean_don_acc[0], mean_don_acc[1])
        # generate file
        df = pd.DataFrame(list(zip(time * 10, rda_md, kappa_2)))
        print(len(df[0]))
        df.to_csv(f'{self.analysis_dir}/explicit_dyes/rkappa.dat', sep='\t', header=False, index=False)
        return df

    def generate_r_kappa_from_dyes(self):

        self.make_dir(f"{self.analysis_dir}/explicit_dyes")

        self.acceptor_dye.get_ids_from_gro(f"{self.analysis_dir}/raw/{self.input_structure_name}_reduced.gro")
        self.donor_dye.get_ids_from_gro(f"{self.analysis_dir}/raw/{self.input_structure_name}_reduced.gro")

        u = mda.Universe(f"{self.analysis_dir}/raw/{self.input_structure_name}_reduced.gro",
                         f"{self.analysis_dir}/raw/{self.input_structure_name}_centered.xtc")

        donor_dipole_coords = [self.get_atoms_coordinates(str(self.donor_dye.dipole_ids[0]), u),
                               self.get_atoms_coordinates(str(self.donor_dye.dipole_ids[1]), u)]
        acceptor_dipole_coords = [self.get_atoms_coordinates(str(self.acceptor_dye.dipole_ids[0]), u),
                                  self.get_atoms_coordinates(str(self.acceptor_dye.dipole_ids[1]), u)]

        print(f"Donor ID:\ncentral c: {self.donor_dye.central_c_id}\ndipole: {self.donor_dye.dipole_ids}")
        print(f"Acceptor ID:\ncentral c: {self.acceptor_dye.central_c_id}\ndipole: {self.acceptor_dye.dipole_ids}")

        central_donor = str(self.donor_dye.central_c_id)
        central_acceptor = str(self.acceptor_dye.central_c_id)

        mean_don_acc = [self.get_atoms_coordinates(central_donor, u), self.get_atoms_coordinates(central_acceptor, u)]
        self.write_rkappa_file_from_dyes(donor_dipole_coords, acceptor_dipole_coords, mean_don_acc)
        self.write_coordinate_file(f"Donor_coords.txt", donor_dipole_coords)
        self.write_coordinate_file(f"Acceptor_coords.txt", acceptor_dipole_coords)


if __name__ == '__main__':
    analysis_paras = {
        "simulation_name": "lukas_poly_u_test_2021er_roteded_restraint",
        "input_structure_name": "lukas_polyu_rotated",
        "Donor_residue_name_number": ["C3W", 10, "A"],
        "Acceptor_residue_name_number": ["C5W", 45, "B"],
    }

    dye_acv_parameter = {
        "Acceptor": {
            "name": "sCy5",
            "linker_length": 20,
            "linker_width": 3.5,
            "dye_radius1": 9.5,
            "dye_radius2": 3,
            "dye_radius3": 1.5,
            "cv_fraction": 0.99,
            "cv_thickness": 3
        },
        "Donor": {
            "name": "sCy3",
            "linker_length": 20,
            "linker_width": 3.5,
            "dye_radius1": 8.0,
            "dye_radius2": 3,
            "dye_radius3": 1.5,
            "cv_fraction": 0.99,
            "cv_thickness": 3,
        },
        "Distance": {"sCy3-sCy5":
                         {"R0": 54,
                          "n_dist": 10 ** 6}
                     }

    }

    print(os.getcwd())

    md_analysis = DataAnalysis(working_dir="/home/felix/Documents/md_lukas_test_2021er",
                               path_sim_results="/home/felix/Documents/md_lukas_test_2021er/lukas_poly_u_test_2021er_roteded_restraint",
                               analysis_parameter=analysis_paras, macv_label_pars=dye_acv_parameter)

    # 1. get all files ready in new analysis folder
    md_analysis.make_data_analysis_results_dirs()
    md_analysis.export_pdb_trajectory(5)
    # 2. calculate r_kappa from explicit dyes
    #md_analysis.generate_r_kappa_from_dyes()
    # 3. calculate r_kappa from macv
    #md_analysis.make_dir(f"{md_analysis.analysis_dir}/macv")
    #md_analysis.remove_dyes_from_trajectory()
    #md_analysis.rewrite_atoms_after_unlabeling()
    #md_analysis.genarate_rkappa_file_from_macv()
    # dye = Dye(("C3W", 63))
    # pdb = dye.get_attechment_id_from_pdb(f"/home/felix/Documents/md_pipeline_testfolder/m_tlr_ub/analysis/raw/m_tlr_ub_1_s1.pdb")
    # print(pdb)
