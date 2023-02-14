from sys import platform
import subprocess
import os


class Modeling:
    def __init__(self, working_dir: str, file_path_sequence: str, modeling_parameter: dict) -> None:
        self.working_dir = working_dir
        self.file_path_sequence = file_path_sequence
        self.modeling_parameter = modeling_parameter
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
        """
        Run a bash command with python subprocess.

        :param command: Bash command as string.
        :return: none
        """
        subprocess.call(["bash", "-c", command], stdout=subprocess.PIPE)

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

    def read_secondary_structure(self):
        sec_struct = ""
        seq = ""
        with open(self.file_path_sequence) as f:
            for i, line in enumerate(f):
                if i == 0:
                    sec_struct = line.strip()
                if i == 1:
                    seq = line.strip()
        return [sec_struct, seq]

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

    def write_rosetta_parameter(self, secondary_structure_file, rosetta_parameter):
        self.make_result_dir("rosetta_results")
        secondary_structure = self.read_secondary_structure()
        if rosetta_parameter["minimize_rna"]:
            file_content = f"{rosetta_parameter['path_to_rosetta']} -nstruct {rosetta_parameter['nstruct']} -sequence '{secondary_structure[1]}'  -secstruct '{secondary_structure[0]}' -silent silent_out.out -minimize_rna {rosetta_parameter['minimize_rna']} -cycles {rosetta_parameter['cycles']}"
        else:
            file_content = f"{rosetta_parameter['path_to_rosetta']} -nstruct {rosetta_parameter['nstruct']} -sequence '{secondary_structure[1]}'  -secstruct '{secondary_structure[0]}' -silent silent_out.out {rosetta_parameter['minimize_rna']} -cycles {rosetta_parameter['cycles']}"

        with open("rosetta_results/FARFAR2.txt", "w") as text_file:
            text_file.write(file_content)

    def predict_3D_structure(self, secondary_structure_file, parameter):
        self.write_rosetta_parameter(secondary_structure_file, parameter)
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

    def extract_pdb(self, number_of_pdb: int) -> None:
        """
        Exports structures from the .out file to a PDB file.

        Parameters
        ----------
        number_of_pdb: Number of PDB structures to be exported from the ensemble

        Returns
        -------

        """

        if platform == "linux" or platform == "linux2":
            self.run_command(f"./scripts/linux/rosetta/extract_pdb.sh -d ./rosetta_results/out/1/ -n {number_of_pdb}"
                             f" -m true -s ./rosetta_results/out/1/silent_out.out")

        elif platform == "darwin":
            self.run_command(f"./scripts/mac_os/rosetta/extract_pdb.sh -d ./rosetta_results/out/1/ -n {number_of_pdb} "
                             f"-m true -s ./rosetta_results/out/1/silent_out.out")

if __name__ == '__main__':

    print(os.getcwd())
    params = {}
    test = Modeling(f"{os.getcwd()}/data", f"{os.getcwd()}/data/RNA_Hairpin.fasta", params)
    test.predict_2d_structure()
