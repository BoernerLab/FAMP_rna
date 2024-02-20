from sys import platform
import subprocess
import os
import pathlib


class Modeling:
    def __init__(self, working_dir: str, file_path_sequence: str, modeling_parameter: dict) -> None:
        self.working_dir = working_dir
        self.file_path_sequence = file_path_sequence
        self.modeling_parameter = modeling_parameter
        self.source_path = pathlib.Path(__file__).parent.resolve()
        self.sequence = self.read_fasta_file()

    @staticmethod
    def check_os():
        """
        Print's the OS you are working on
        :return: None
        """
        if platform == "linux" or platform == "linux2":
            print("You are working under linux")
        elif platform == "darwin":
            print("Your are working under MacOS")
        elif platform == "":
            print("You are working under Windows. The pipeline has not yet "
                  "been developed for this operating system.")

    @staticmethod
    def run_command_win(command: str, cmd_in: bytes):
        """
        Run a command in bash with user input. Calls python subprocess module.

        :param command: bash command
        :param cmd_in: commandline input in bytes
        :return: none
        """
        process = subprocess.Popen(["bash", "-c", command], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        process.communicate(input=cmd_in)
        process.wait()

    @staticmethod
    def run_command(command: str):
        try:
            output = subprocess.check_output(
                command, stderr=subprocess.STDOUT, shell=True,
                universal_newlines=True)
        except subprocess.CalledProcessError as exc:
            print("Status : FAIL", exc.returncode, exc.output)


        #process = subprocess.Popen(["bash", "-c", command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, text=True)
        #while process.stdout.readable():
        #    line = process.stdout.readline()

        #    if not line:
        #        break

        #    print(line)


    def run_command_rosetta(command: str):
        counter = 1
        line_finders = ["--------------", "Scores", "fa_atr", "fa_rep", "fa_intra_rep", "lk_nonpolar", "fa_elec_rna_phos_phos", "rna_torsion", "suiteness_bonus", "rna_sugar_close", "fa_stack", "stack_elec", "geom_sol_fast", "hbond_sr_bb_sc", "hbond_lr_bb_sc", "hbond_sc", "ref", "free_suite", "free_2HOprime", "intermol", "other_pose", "loop_close", "linear_chainbreak"]

        process = subprocess.Popen(["bash", "-c", command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, text=True)
        print(f"modeling structure {counter}")
        while process.stdout.readable():
            line = process.stdout.readline()

            if not line:
                break
            if any(s in line for s in line_finders):
                print(line)
            if "S_" in line:
                print(line)
                counter=+1
                print(f"modeling structure {counter}")

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

    def read_fasta_file(self) -> str:
        """
        Reads in the sequence of a fasta file and make the latter's lowe case.
        Writes the lowercase letters into the file

        :return: RNA Sequence in lower case letters
        """
        identifier = ""
        seq = ""
        with open(self.file_path_sequence) as f:
            for i, line in enumerate(f):
                if i == 0:
                    identifier = line
                if i == 1:
                    seq = line.lower()
        print(f"Read in Sequence: {identifier}")
        with open(self.file_path_sequence, "w") as text_file:
            text_file.write(identifier)
            text_file.write(seq)
        return seq

    def reduce_sds_file(self, path_result_file: str):
        """
        Reduces a secondary structure file to only the sequence and the dot bracket formatted 2D structure
        Parameters
        ----------
        path_result_file (str) Path to the file into which the results should be written.

        Returns
        -------
        None
        """

        dot_bracket = ""
        seq = ""
        with open(f"{self.working_dir}/secondary_prediction/{path_result_file}") as f:
            next(f)
            for i, line in enumerate(f):
                if i == 1:
                    print("1", line)
                    dot_bracket = line.split()[0]
                if i == 0:
                    print("2", line)
                    seq = line

        with open(f"{self.working_dir}/secondary_prediction/dot_bracket.secstruct", "w") as text_file:
            text_file.write(dot_bracket + "\n")
            text_file.write(seq)
        # return [dot_bracket, seq[:-1]]

    def predict_2d_structure(self):
        self.make_result_dir("secondary_prediction")
        # subprocess.run(["bash","-c"," mkdir sds_prediction"], capture_output= True)
        self.run_command(f"RNAfold -i {self.file_path_sequence} --noPS "
                         f"> {self.working_dir}/secondary_prediction/RNA_fold_output.txt")
        self.reduce_sds_file("RNA_fold_output.txt")

    @staticmethod
    def read_secondary_structure(path_2d_structure_file):
        sec_struct = ""
        seq = ""
        with open(path_2d_structure_file) as f:
            for i, line in enumerate(f):
                if i == 0:
                    sec_struct = line.strip()
                if i == 1:
                    seq = line.strip()
        return [sec_struct, seq]

    def write_rna_denovo_parameter(self, path_2d_structure_file: str) -> str:
        """
        Writes a txt file with the rosetta parameter to run rosetta scripts
        Parameters
        ----------
        path_2d_structure_file

        Returns
        -------

        """
        self.make_result_dir("rosetta_results")
        secondary_structure = self.read_secondary_structure(path_2d_structure_file)
        flags = f"{self.modeling_parameter['path_to_rosetta']} " \
                f"-nstruct {self.modeling_parameter['nstruct']} " \
                f"-sequence '{secondary_structure[1]}'  " \
                f"-secstruct '{secondary_structure[0]}' " \
                f"-silent {self.working_dir}/rosetta_results/silent_out.out " \
                f"-minimize_rna {self.modeling_parameter['minimize_rna']} " \
                f"-cycles {self.modeling_parameter['cycles']}" \

        if self.modeling_parameter["overwrite"]:
            flags = f"{flags} -overwrite"

        return flags

    def predict_3d_structure(self, path_2d_structure_file: str) -> None:
        """
        Uses the rosetta modul rna_denovo to predict tertiary structures

        :param path_2d_structure_file: Path to the secondary structure file.

        :return: None
        """
        command = self.write_rna_denovo_parameter(path_2d_structure_file)
        self.run_command(command)

    def extract_pdb(self, number_of_pdb: int) -> None:
        """
        Exports structures from the .out file to a PDB file.

        :param: number_of_pdb: Number of PDB structures to be exported from the ensemble

        :return: None

        """
        #TODO not check the parameters. Check the number of structures in .out
        if number_of_pdb > self.modeling_parameter["nstruct"]:
            number_of_pdb = self.modeling_parameter["nstruct"]

        script_source_dir_linux = f"{self.source_path}/scripts/linux/rosetta/"
        script_source_dir_mac = f"{self.source_path}/scripts/mac_os/rosetta/"

        if platform == "linux" or platform == "linux2":
            self.run_command(f"{script_source_dir_linux}linux_extract_pdb.sh "
                             f"-d {self.working_dir}/rosetta_results/ "
                             f"-n {number_of_pdb}"
                             f" -m true "
                             f"-s {self.working_dir}/rosetta_results/silent_out.out")

        elif platform == "darwin":
            self.run_command(f"{script_source_dir_mac}mac_extract_pdb.sh "
                             f"-d {self.working_dir}/rosetta_results/ "
                             f"-n {number_of_pdb} "
                             f"-m true "
                             f"-s {self.working_dir}/rosetta_results/silent_out.out")


if __name__ == '__main__':
    rosetta_parameter = {
        "path_to_rosetta": "rna_denovo.default.linuxgccrelease",
        "nstruct": 500,
        "fasta": "BTL.fasta",
        "minimize_rna": True,
        "cycles": 10000,
        "overwrite": True,
    }

    print(os.getcwd())
    params = {}
    BTL_modeling = Modeling(working_dir=f"/home/felix/Documents/Rosetta_KLTL_ensemble_unbound",
                    file_path_sequence=f"/home/felix/Documents/Rosetta_KLTL_ensemble_unbound/BTL.fasta",
                    modeling_parameter=rosetta_parameter)
    #BTL_modeling.predict_2d_structure()
    BTL_modeling.predict_3d_structure(path_2d_structure_file=f"/home/felix/Documents/Rosetta_KLTL_ensemble_unbound/secondary_prediction/dot_bracket.secstruct")
    BTL_modeling.extract_pdb(500)
