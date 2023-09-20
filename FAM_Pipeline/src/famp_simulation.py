import subprocess
import re
import os
import shutil
from exceptions import GromacsMdrunError, GromacsEditconfError, GromacsGenionError, GromacsGromppError, GromacsSolvateError, GromacsPdb2gmxError


class MDSimulation:
    def __init__(self, working_dir: str, file_path_input: str, md_parameter: dict) -> None:
        self.working_dir = self.create_working_dir(working_dir)
        self.file_path_input = file_path_input
        self.md_parameter = md_parameter
        self.path_simulation_folder = self.get_simulation_path()
        self.input_structure_name = self.get_input_structure_name()

    @staticmethod
    def run_gromacs_command(command):
        module_name = command.split()[1]
        error_massege = f"{module_name} failed! \n" \
                        f"Command: {command} \n" \
                        f"Please read the GROMACS error massage for further trouble shooting!"

        try:
            output = subprocess.check_output(
                command, stderr=subprocess.STDOUT, shell=True,
                universal_newlines=True)
        except subprocess.CalledProcessError as exc:
            print("Status : FAIL", exc.returncode, exc.output)
            if module_name == "pdb2gmx":
                print(error_massege)
                raise GromacsPdb2gmxError
            elif module_name == "editconf":
                print(error_massege)
                raise GromacsEditconfError
            elif module_name == "solvate":
                print(error_massege)
                raise GromacsSolvateError
            elif module_name == "grompp":
                print(error_massege)
                raise GromacsGromppError
            elif module_name == "genion":
                print(error_massege)
                raise GromacsGenionError
            elif module_name == "mdrun":
                print(error_massege)
                raise GromacsMdrunError
            else:
                print(error_massege)
                raise Exception
        else:
            print("Output: \n{}\n".format(output))

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

    @staticmethod
    def make_ndx_of_SOL(gro_file: str, output_file: str):
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
                if "SOL" in line:
                    if len(line.split()) == 6:
                        atom = line.split()[2]
                        atoms.append(atom)
                    else:
                        atom = line.split()[1]
                        number = re.split("OW|HW1|HW2|MW", atom)
                        atoms.append(number[1])

        with open(output_file, 'w') as the_file:
            the_file.write('[SOL]\n')
            counter = 1
            for element in atoms:
                if counter == 15:
                    the_file.write(f"{element.rjust(6)}\n")
                    counter = 1
                else:
                    the_file.write(f"{element.rjust(6)}")
                    counter = counter + 1
            the_file.write('\n')

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

        water_file = ""
        if self.md_parameter["water_model"] == "tip3p":
            water_file = "spc216.gro"
        elif self.md_parameter["water_model"] == "tip4p":
            water_file = "tip4p.gro"
        else:
            raise ValueError("Please use tip3p or tip4p as water model. Other water models are currently not implemented within the pipeline.")

        # os.chdir(working_dir_path + f"/{dir}")
        # print(os.getcwd())
        self.make_result_dir(f"{self.md_parameter['simulation_name']}/em")
        self.run_gromacs_command(
            f"gmx pdb2gmx "
            f"-f {self.path_simulation_folder}/{self.input_structure_name}.pdb "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-i {self.path_simulation_folder}/em/{self.input_structure_name}.itp "
            f"-missing "
            f"-ignh "
            f"-ff amber14sb_OL15 "
            f"-water {self.md_parameter['water_model']}")

        self.run_gromacs_command(
            f"gmx editconf "
            f"-f {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-bt dodecahedron "
            f"-d {self.md_parameter['dist_to_box[nm]']}")

        self.run_gromacs_command(
            f"gmx solvate "
            f"-cp {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-cs {water_file} "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top ")

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/em.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        self.make_ndx_of_SOL(f"{self.path_simulation_folder}/em/{self.input_structure_name}.gro",
                             f"{self.path_simulation_folder}/em/SOL.ndx")

        self.run_gromacs_command(
            f"gmx genion "
            f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-nname Cl "
            f"-pname K "
            f"-neutral "
            f"-n {self.path_simulation_folder}/em/SOL.ndx")

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/em.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        if simulation_parameter["c_magnesium_ions[mol/l]"] > 0:

            self.make_ndx_of_SOL(f"{self.path_simulation_folder}/em/{self.input_structure_name}.gro",
                                 f"{self.path_simulation_folder}/em/SOL.ndx")

            self.run_gromacs_command(
                f"gmx genion "
                f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
                f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
                f"-nname Cl "
                f"-pname MG "
                f"-pq 2 "
                f"-conc {simulation_parameter['c_magnesium_ions[mol/l]']} "
                f"-n {self.path_simulation_folder}/em/SOL.ndx")

            self.run_gromacs_command(
                f"gmx grompp "
                f"-f {self.path_simulation_folder}/mdp/em.mdp "
                f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
                f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
                f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
                f"-maxwarn 2")

        self.run_gromacs_command(
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

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/nvt.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-r {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/nvt/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/nvt/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        self.run_gromacs_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/nvt/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/nvt/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/nvt/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/nvt/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/nvt/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/nvt/{self.input_structure_name}.log")

        self.make_result_dir(f"{self.md_parameter['simulation_name']}/npt")

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/npt.mdp "
            f"-c {self.path_simulation_folder}/nvt/{self.input_structure_name}.gro "
            f"-r {self.path_simulation_folder}/nvt/{self.input_structure_name}.gro "
            f"-t {self.path_simulation_folder}/nvt/{self.input_structure_name}.cpt "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/npt/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/npt/{self.input_structure_name}.mdp "
            f"-maxwarn 2")
        self.run_gromacs_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/npt/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/npt/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/npt/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/npt/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/npt/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/npt/{self.input_structure_name}.log")

        self.make_result_dir(f"{self.md_parameter['simulation_name']}/md0")

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/md0.mdp "
            f"-c {self.path_simulation_folder}/npt/{self.input_structure_name}.gro "
            f"-t {self.path_simulation_folder}/npt/{self.input_structure_name}.cpt "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/md0/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/md0/{self.input_structure_name}.mdp  "
            f"-maxwarn 2")
        self.run_gromacs_command(
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
        "simulation_name": "m_tlr_ub",
        "c_magnesium_ions[mol/l]": 0.02,
        "simulation_time[ns]": 0.1,
        "temperature[°C]": 25,
        "dist_to_box[nm]": "1.25",
        "water_model": "tip3p"
    }
    print(os.getcwd())
    hairpin_labeled = MDSimulation(working_dir=f"/home/felix/Documents/md_pipeline_testfolder",
                                   file_path_input=f"/home/felix/Documents/md_pipeline_testfolder/m_tlr_ub_1.pdb",
                                   md_parameter=simulation_parameter)

    hairpin_labeled.prepare_new_md_run()
    hairpin_labeled.update_parameter()
    hairpin_labeled.solvate_molecule()
    hairpin_labeled.run_simulation_steps()
