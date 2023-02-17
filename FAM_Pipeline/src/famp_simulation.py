import subprocess
import re
import os
import shutil


class MDSimulation:
    def __init__(self, working_dir: str, file_path_input: str, md_parameter: dict) -> None:
        self.working_dir = self.create_working_dir(working_dir)
        self.file_path_input = file_path_input
        self.md_parameter = md_parameter
        self.path_simulation_folder = self.get_simulation_path()
        self.input_structure_name = self.get_input_structure_name()

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
        process = subprocess.Popen(["bash", "-c", command], stdout=subprocess.PIPE, text=True)
        while process.stdout.readable():
            line = process.stdout.readline()

            if not line:
                break

            print(line.strip())

    @staticmethod
    def sim_time_to_steps(sim_time):
        return int(1000000 * sim_time / 2)

    @staticmethod
    def degree_to_kelvin(degree):
        return degree + 273

    @staticmethod
    def create_working_dir(working_dir):
        if os.path.isdir(working_dir):
            print("The specified folder exists.")
        else:
            os.mkdir(working_dir)
            print("The specified folder does not exist but was created.")
        return working_dir

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

    def get_simulation_path(self):
        return f"{self.working_dir}/{self.md_parameter['simulation_name']}"

    def get_input_structure_name(self):
        file_name = os.path.basename(os.path.normpath(self.file_path_input))
        structure_name = file_name[:-4]
        return structure_name

    def change_temperature_in_nvt(self, temperature):
        temp_in_k = self.degree_to_kelvin(temperature)
        for source_dir in ["nvt", "md0"]:
            content = []
            with open(f"{self.path_simulation_folder}/mdp/{source_dir}.mdp", 'r') as f:
                for i, line in enumerate(f):
                    line = line.strip()
                    if line.startswith(('ref_t', 'gen_temp')):
                        # print(re.sub(pattern = "[0-9]+", repl = str(temp_in_K), string=line))
                        content.append(re.sub(pattern="[0-9]+", repl=str(temp_in_k), string=line))
                    else:
                        content.append(line)

            # print(content)
            with open(f"{self.path_simulation_folder}/mdp/{source_dir}.mdp", 'w') as f:
                for line in content:
                    f.write("%s\n" % line)

    def change_temperature_in_npt(self, temperature):
        temp_in_k = self.degree_to_kelvin(temperature)
        content = []
        with open(f"{self.path_simulation_folder}/mdp/npt.mdp", 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if line.startswith('ref_t'):
                    # print(re.sub(pattern = "[0-9]+", repl = str(temp_in_K), string=line))
                    content.append(re.sub(pattern="[0-9]+", repl=str(temp_in_k), string=line))
                else:
                    content.append(line)

        # print(content)
        with open(f"{self.path_simulation_folder}/mdp/npt.mdp", 'w') as f:
            for line in content:
                f.write("%s\n" % line)

    def change_sim_time_in_md0(self, time):
        simulation_steps = self.sim_time_to_steps(time)
        content = []
        with open(f"{self.path_simulation_folder}/mdp/md0.mdp", 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if line.startswith('nsteps'):
                    # print(re.sub(pattern = "[0-9]+", repl = str(simulation_steps), string=line))
                    content.append(re.sub(pattern="[0-9]+", repl=str(simulation_steps), string=line))
                else:
                    content.append(line)

        with open(f"{self.path_simulation_folder}/mdp/md0.mdp", 'w') as f:
            for line in content:
                f.write("%s\n" % line)

    def copy_files_to_sim_dir(self):
        src_folder = "./scripts/gromacs"
        dst_folder = self.working_dir + f"/{self.md_parameter['simulation_name']}"

        if os.path.exists(dst_folder) and os.path.isdir(dst_folder):
            print(
                "MD run already exists. To make a new Simulation change the Name od the Simulation in the MD parameter")
        else:
            self.make_result_dir(self.md_parameter['simulation_name'])

        if os.path.exists(dst_folder + "/amber14sb_OL15.ff") and os.path.isdir(dst_folder + "/amber14sb_OL15.ff"):
            pass
        else:
            shutil.copytree(src_folder + "/amber14sb_OL15.ff", dst_folder + "/amber14sb_OL15.ff")

        if os.path.exists(dst_folder + "/mdp") and os.path.isdir(dst_folder + "/mdp"):
            pass
        else:
            shutil.copytree(src_folder + "/mdp", dst_folder + "/mdp")

    def prepare_new_md_run(self):
        self.copy_files_to_sim_dir()
        self.copy_input_model()
        # Ändern der Parameter in den mdp files
        self.change_temperature_in_nvt(simulation_parameter["temperature[°C]"])
        self.change_temperature_in_npt(simulation_parameter["temperature[°C]"])
        self.change_sim_time_in_md0(simulation_parameter["simulation_time[ns]"])

    def update_parameter(self):
        self.change_temperature_in_nvt(simulation_parameter["temperature[°C]"])
        self.change_temperature_in_npt(simulation_parameter["temperature[°C]"])
        self.change_sim_time_in_md0(simulation_parameter["simulation_time[ns]"])

    def copy_input_model(self) -> None:
        """
        Copies the modeling result structure to the MD simulation directory.
        """
        self.make_result_dir(self.md_parameter["simulation_name"])
        input_file_name = os.path.basename(os.path.normpath(self.file_path_input))
        shutil.copy(f"{self.file_path_input}",
                    f"{self.working_dir}/{self.md_parameter['simulation_name']}/{input_file_name}")

    def solvate_molecule(self) -> None:
        """
        Running bash commands with python subprocess to solvate MD run with GROMACS. Reference solvate.sh.

        :return: none
        """
        # os.chdir(working_dir_path + f"/{dir}")
        # print(os.getcwd())
        self.make_result_dir(f"{self.md_parameter['simulation_name']}/em")
        self.run_command_win(
            f"gmx pdb2gmx "
            f"-f {self.path_simulation_folder}/{self.input_structure_name}.pdb "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-i {self.path_simulation_folder}/em/{self.input_structure_name}.itp "
            f"-missing "
            f"-ignh",
            b"1\n 3\n")

        self.run_command(f"gmx editconf "
                         f"-f {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                         f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                         f"-bt dodecahedron "
                         f"-d {self.md_parameter['dist_to_box[nm]']}")

        self.run_command(
            f"gmx solvate "
            f"-cp {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-cs tip4p.gro "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top ")

        self.run_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/em.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        self.run_command_win(
            f"gmx genion "
            f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-nname Cl "
            f"-pname K "
            f"-neutral",
            b"3\n")

        self.run_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/em.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        if simulation_parameter["c_magnesium_ions[mol/l]"] > 0:
            self.run_command_win(
                f"gmx genion "
                f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
                f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
                f"-nname Cl "
                f"-pname MG "
                f"-pq 2 "
                f"-conc {simulation_parameter['c_magnesium_ions[mol/l]']}",
                b"4\n")

            self.run_command(
                f"gmx grompp "
                f"-f {self.path_simulation_folder}/mdp/em.mdp "
                f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
                f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
                f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
                f"-maxwarn 2")

        self.run_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.trr "
            f"-e {self.path_simulation_folder}/em/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/em/{self.input_structure_name}.log")

    def run_simulation_steps(self) -> None:
        """
        Running bash commands with python subprocess to make a single MD run with GROMACS. Reference single_run.sh

        :return: None
        """
        self.make_result_dir(f"{self.md_parameter['simulation_name']}/nvt")

        self.run_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/nvt.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-r {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/nvt/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/nvt/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        self.run_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/nvt/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/nvt/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/nvt/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/nvt/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/nvt/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/nvt/{self.input_structure_name}.log")

        self.make_result_dir(f"{self.md_parameter['simulation_name']}/npt")

        self.run_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/npt.mdp "
            f"-c {self.path_simulation_folder}/nvt/{self.input_structure_name}.gro "
            f"-r {self.path_simulation_folder}/nvt/{self.input_structure_name}.gro "
            f"-t {self.path_simulation_folder}/nvt/{self.input_structure_name}.cpt "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/npt/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/npt/{self.input_structure_name}.mdp "
            f"-maxwarn 2")
        self.run_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/npt/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/npt/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/npt/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/npt/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/npt/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/npt/{self.input_structure_name}.log")

        self.make_result_dir(f"{self.md_parameter['simulation_name']}/md0")

        self.run_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/md0.mdp "
            f"-c {self.path_simulation_folder}/npt/{self.input_structure_name}.gro "
            f"-t {self.path_simulation_folder}/npt/{self.input_structure_name}.cpt "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/md0/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/md0/{self.input_structure_name}.mdp  "
            f"-maxwarn 2")
        self.run_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/md0/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/md0/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/md0/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/md0/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/md0/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/md0/{self.input_structure_name}.log")

        # os.chdir(working_dir_path)


if __name__ == '__main__':
    simulation_parameter = {
        "simulation_name": "hairpin_labeled",
        "c_magnesium_ions[mol/l]": 0.00,
        "simulation_time[ns]": 0.3,
        "temperature[°C]": 25,
        "dist_to_box[nm]": "1",
    }
    print(os.getcwd())
    hairpin_labeled = MDSimulation(working_dir=f"{os.getcwd()}/data",
                                   file_path_input=f"{os.getcwd()}/data/rosetta_results/silent_out.pdb",
                                   md_parameter=simulation_parameter)

    hairpin_labeled.prepare_new_md_run()
    hairpin_labeled.update_parameter()
    hairpin_labeled.solvate_molecule()
    hairpin_labeled.run_simulation_steps()
