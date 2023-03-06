import subprocess
import os
import MDAnalysis as mda
import numpy as np
from numpy import linalg as LA
import re
import fileinput
import shutil
import pandas as pd
import mdtraj as md
import fretraj as ft


class DataAnalysis:
    def __init__(self, working_dir, path_sim_results: str, analysis_parameter: dict, macv_label_pars: dict) -> None:
        self.working_dir = working_dir
        self.path_sim_results = path_sim_results
        self.analysis_parameter = analysis_parameter
        self.macv_label_pars = macv_label_pars
        self.analysis_dir = f"{self.path_sim_results}/analysis"
        self.input_structure_name = self.analysis_parameter["input_structure_name"]
        self.md_traj = None
        self.fret_macv = None

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

    def make_result_dir(self, directory_name):
        """
        Creates a directory where the results for operations are stored.
        The directory is created in the defined working directory

        :param: directory_name - Name of the new directory
        :return: None
        """
        result_dir = f"{self.working_dir}/{directory_name}"
        try:
            os.mkdir(result_dir)
        except FileExistsError:
            print(f"Results can be found in: {result_dir}")
        except OSError:
            print(f"Failed to create {result_dir}")
        else:
            print(f"Successfully created the directory {result_dir}. Results can be found there")

    def set_md_traj(self):
        traj = md.load(f'{self.analysis_dir}/raw/{self.input_structure_name}_unlabeled.xtc',
                       top=f'{self.analysis_dir}/raw/{self.input_structure_name}_unlabeled_s1.pdb')
        self.md_traj = traj

    # -----------------------------------------------------------------------------------------------------------------------

    def reduce_center_xtc(self):
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
            f"gmx trjconv -f {md_dir}/md0/{sim_name}.xtc -s {md_dir}/md0/{sim_name}.tpr -o {md_dir}/md0/{sim_name}_centered.xtc -n {self.analysis_dir}/Index_Files/RNA.ndx -pbc mol -center")
        self.run_command(
            f"gmx trjconv -f {md_dir}/md0/{sim_name}.xtc -s {md_dir}/md0/{sim_name}.tpr -o {md_dir}/md0/{sim_name}_s1.pdb -n {self.analysis_dir}/Index_Files/RNA.ndx -pbc mol -center -b 1 -e 10")

    def export_pdb_trajectory(self, time_steps):
        """
        Exports a pdb trajectory from a GROMACS simulation file.

        :params: time_steps: step width for exporting states.
        :return: none
        """
        sim_name = self.analysis_parameter['input_structure_name']
        self.run_command(
            f"gmx trjconv -f {self.analysis_dir}/raw/{sim_name}_centered.xtc -s {self.analysis_dir}/raw/{sim_name}.tpr -o {self.analysis_dir}/raw/{sim_name}_traj.pdb -n {self.analysis_dir}/Index_Files/RNA.ndx -pbc mol -dt {time_steps} -center")

    def make_data_analysis_results_dirs(self):
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
        self.reduce_center_xtc()
        self.copy_files_to_raw()

    def copy_files_to_raw(self):
        """Prepare the directory for MD analysis procedure"""
        src_folder = f"{self.path_sim_results}/md0"
        dst_folder = f"{self.analysis_dir}/raw"
        sim_name = self.analysis_parameter["input_structure_name"]
        shutil.copy(src_folder + f"/{sim_name}_centered.xtc", dst_folder + f"/{sim_name}_centered.xtc")
        shutil.copy(src_folder + f"/{sim_name}.gro", dst_folder + f"/{sim_name}.gro")
        shutil.copy(src_folder + f"/{sim_name}.tpr", dst_folder + f"/{sim_name}.tpr")
        shutil.copy(src_folder + f"/{sim_name}_s1.pdb", dst_folder + f"/{sim_name}_s1.pdb")

    # -----------------------------------------------------------------------------------------------------------------------

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

        self.rewrite_atoms_after_unlabeling()

    def get_selected_frames(self):
        time = self.md_traj.time[-1]
        step = self.md_traj.timestep
        max_time = round(time / step, 0)
        time_step = 1
        ts = max_time / 1000
        if ts >= 1:
            time_step = round(ts, 0)
        else:
            time_step = 1

        s_frames = [int(max_time + 1), int(time_step)]
        return s_frames

    def calculate_macv(self):
        s_frames = self.get_selected_frames()
        selected_frames = range(0, s_frames[0], s_frames[1])
        print(s_frames)
        fret = ft.cloud.pipeline_frames(self.md_traj, 'Cy5-10-C5', "Cy3-65-O3'", self.macv_label_pars, selected_frames, 'Cy3-Cy5')
        ft.cloud.save_obj(f'{self.analysis_dir}/macv/{self.input_structure_name}_macv.pkl', fret)
        return fret

    def load_macv(self):
        fret = ft.cloud.load_obj(f'{self.analysis_dir}/macv/{self.input_structure_name}_macv.pkl')
        return fret

    def write_rkappa_file_from_macv(self):
        fret_traj = ft.cloud.Trajectory(self.fret_macv, timestep=self.fret_macv.timestep, kappasquare=0.66)
        fret_traj.save_traj(f'{self.analysis_dir}/macv/R_kappa_ACV.dat', format='txt', R_kappa_only=True, units='nm',
                            header=False)
        fret_traj.dataframe.head()

    def genarate_rkappa_file_from_macv(self, calculate_macv=True):
        self.make_dir(f"{self.analysis_dir}/macv")
        self.remove_dyes_from_trajectory()
        self.rewrite_atoms_after_unlabeling()
        self.set_md_traj()

        if calculate_macv:
            self.fret_macv = self.calculate_macv()
        else:
            self.fret_macv = self.load_macv()

        self.write_rkappa_file_from_macv()

    # -----------------------------------------------------------------------------------------------------------------------

    def get_atoms_coordinates(self, atom_id, universe):
        universe.trajectory[0]
        acc_c14 = universe.select_atoms(f'id {atom_id}')
        acc_coords = []
        for ts in universe.trajectory:
            coords = list(np.squeeze(acc_c14.positions / 10))
            acc_coords.append(coords)
        return np.array(acc_coords)

    def get_atom_ids(self, dye_filter: list) -> pd.DataFrame:
        """
        Return atom id's from a given residue number and a residue name.

        :param dye_filter: List of Number of residue (10) and Name of residue (C5W or RU)
        :return: Dataframe of residue id's with other information.
        """

        if os.path.isfile(f"{self.analysis_dir}/raw/{self.input_structure_name}.gro"):
            df_list = []
            for filter_set in dye_filter:
                with open(f"{self.analysis_dir}/raw/{self.input_structure_name}.gro", 'r') as f:
                    next(f)
                    next(f)
                    for i, line in enumerate(f):
                        if not any(value in line for value in "SOL"):
                            split_line = line.split()
                            if split_line[0].startswith(filter_set[0]):
                                if split_line[1] == filter_set[1]:
                                    df_list.append([split_line[0], split_line[1], split_line[2]])
            return pd.DataFrame(df_list, columns=["Num/Res", "Atom name", "ID"])
        else:
            print(
                f"There should be a file named {self.input_structure_name}.gro in the folder {self.analysis_dir}/raw/."
                f" Please check if this file and folder exist. Rename the file if necessary.")

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

    def write_coordinate_file(self, file_name: str, dipole):
        """
        Creates a file of coordinates for dipoles of dyes. Needed when anisotropy calculations are performed.

        :param file_name: Path and name of file, where the coordinates should be written in.
        :param dipole: Coordinates of dipole: [[[x,y,z]][[x,y,z]]]
        :return: none
        """
        time = np.arange(0, len(dipole[0]) + 9, 10, dtype=int)
        df = pd.DataFrame(list(
            zip(time, dipole[0][:, 0], dipole[0][:, 1], dipole[0][:, 2], dipole[1][:, 0], dipole[1][:, 1],
                dipole[1][:, 1])))
        print(len(df[0]))
        df.to_csv(f"{self.analysis_dir}/fluorburst/{file_name}", sep='\t', header=False, index=False)

        for i in range(0, 13):
            self.line_prepender(f"{self.analysis_dir}/fluorburst/{file_name}", "#")

    def write_rkappa_file_from_dyes(self, don_dipole, acc_dipole, mean_don_acc):
        """
        Creates a distance and 𝜅^2 file of a trajectory with explicit dyes.

        The 𝜅^2 values are calculated by the coordinates of the donor and acceptor dipoles. Then the mean inter dye
        distance ist calculated by the coordinates of mean donor and acceptor atoms. Then the time steps 𝜅^2 and mean
        dye distances are combined to a dataframe and written to a file in the fluorburst folder.

        :param don_dipole: Coordinates of donor dipole: [[[x,y,z]][[x,y,z]]]
        :param acc_dipole: Coordinates of acceptor dipole: [[[x,y,z]][[x,y,z]]]
        :param mean_don_acc: Coordinates of donor dipole: [[[x,y,z]][[x,y,z]]]
        :return: none
        """
        kappa_2 = self.calculate_kappa_2(don_dipole, acc_dipole)
        time = np.arange(0, len(don_dipole[0]) + 9, 10, dtype=int)
        # print(len(time))
        rda_md = self.calculate_inter_dye_distance(mean_don_acc[0], mean_don_acc[1])
        # Datei generieren
        df = pd.DataFrame(list(zip(time, rda_md, kappa_2)))
        print(len(df[0]))
        df.to_csv(f'{self.analysis_dir}/fluorburst/rkappa.dat', sep='\t', header=False, index=False)
        return df

    def generate_r_kappa_from_dyes(self):
        self.reduce_gro_file()
        self.make_dir(f"{self.analysis_dir}/fluorburst")
        u = mda.Universe(f"{self.analysis_dir}/raw/{self.input_structure_name}_reduced.gro",
                         f"{self.analysis_dir}/raw/{self.input_structure_name}_centered.xtc")

        donor_dipol = self.analysis_parameter["donor_dipole"]
        acceptor_dipol = self.analysis_parameter["acceptor_dipole"]

        donor_dipole_coords = [self.get_atoms_coordinates(str(donor_dipol[0]), u),
                               self.get_atoms_coordinates(str(donor_dipol[1]), u)]
        acceptor_dipole_coords = [self.get_atoms_coordinates(str(acceptor_dipol[0]), u),
                                  self.get_atoms_coordinates(str(acceptor_dipol[1]), u)]

        central_donor = str(self.analysis_parameter["mean_donor_atom"])
        central_acceptor = str(self.analysis_parameter["mean_acceptor_atom"])

        mean_don_acc = [self.get_atoms_coordinates(central_donor, u), self.get_atoms_coordinates(central_acceptor, u)]
        self.write_rkappa_file_from_dyes(donor_dipole_coords, acceptor_dipole_coords, mean_don_acc)
        self.write_coordinate_file(f"Donor_coords.txt", donor_dipole_coords)
        self.write_coordinate_file(f"Acceptor_coords.txt", acceptor_dipole_coords)


if __name__ == '__main__':
    analysis_paras = {
        "simulation_name": "cryo_em_model_labeled",
        "input_structure_name": "input",
        "mean_donor_atom": 2126,
        "donor_dipole": [2123, 2155],
        "mean_acceptor_atom": 304,
        "acceptor_dipole": [311, 334]
    }

    labels = {"Position":
                  {"Cy5-10-C5":
                       {"attach_id": 310,
                        "mol_selection": "all",
                        "linker_length": 21,
                        "linker_width": 3.5,
                        "dye_radius1": 9.5,
                        "dye_radius2": 3,
                        "dye_radius3": 1.5,
                        "cv_fraction": 0.25,
                        "cv_thickness": 3,
                        "use_LabelLib": False,
                        "grid_spacing": 1.0,
                        "simulation_type": "AV3",
                        "state": (int, 1),
                        "frame_mdtraj": (int, 0),
                        "contour_level_AV": ((int, float), 0),
                        "contour_level_CV": ((int, float), 0.7),
                        "b_factor": (int, 100),
                        "gaussian_resolution": (int, 2),
                        "grid_buffer": ((int, float), 2.0),
                        "transparent_AV": (bool, True)
                        },
                   "Cy3-65-O3'":
                       {"attach_id": 2052,
                        "mol_selection": "all",
                        "linker_length": 20.5,
                        "linker_width": 3.5,
                        "dye_radius1": 8,
                        "dye_radius2": 3,
                        "dye_radius3": 1.5,
                        "cv_fraction": 0.25,
                        "cv_thickness": 3,
                         "use_LabelLib": False,
                        "grid_spacing": 1.0,
                        "simulation_type": "AV3",
                        "state": (int, 1),
                        "frame_mdtraj": (int, 0),
                        "contour_level_AV": ((int, float), 0),
                        "contour_level_CV": ((int, float), 0.7),
                        "b_factor": (int, 100),
                        "gaussian_resolution": (int, 2),
                        "grid_buffer": ((int, float), 2.0),
                        "transparent_AV": (bool, True)},
                   },
              "Distance": {"Cy3-Cy5":
                               {"R0": 54,
                                "n_dist": (int, 10 ** 6)}
                           }
              }

    filter_par = [["10", "C14"], ["10", "C2"], ["10", "C32"], ["65", "C14"], ["65", "C2"], ["65", "C11"]]

    print(os.getcwd())
    md_analysis = DataAnalysis(
        working_dir=f"/home/felix/Documents/md_BTL_Ros_PyM_04_04/md_CryoEM_without_restraints_labeled",
        path_sim_results=f"/home/felix/Documents/md_BTL_Ros_PyM_04_04/md_CryoEM_without_restraints_labeled",
        analysis_parameter=analysis_paras,
        macv_label_pars=labels)

    # md_analysis.make_data_analysis_results_dirs()
    # md_analysis.reduce_center_xtc()
    # md_analysis.export_pdb_trajectory(1000)
    # md_analysis.generate_r_kappa_from_dyes()
    md_analysis.genarate_rkappa_file_from_macv()

