import subprocess
import re
import os
import shutil


class MDSimulation:
    def __init__(self, working_dir: str, file_path_input: str, md_parameter: dict) -> None:
        self.working_dir = working_dir
        self.file_path_input = file_path_input
        self.md_parameter = md_parameter

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
    def grad_to_kelvin(grad):
        return grad + 273

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

    def change_tamperature_in_nvt(self, temperature, dir):
        temp_in_K = self.grad_to_kelvin(temperature)
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

    def change_tamperature_in_npt(self, temperature, dir):
        temp_in_K = self.grad_to_kelvin(temperature)
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

    def change_sim_time_in_md0(self, time, dir):
        simulation_steps = self.sim_time_to_steps(time)
        content = []
        with open(f"{dir}/mdp/md0.mdp", 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if line.startswith('nsteps'):
                    # print(re.sub(pattern = "[0-9]+", repl = str(simulation_steps), string=line))
                    content.append(re.sub(pattern="[0-9]+", repl=str(simulation_steps), string=line))
                else:
                    content.append(line)

        print(content)
        with open(f"{dir}/mdp/md0.mdp", 'w') as f:
            for l in content:
                f.write("%s\n" % l)

    def copy_files_to_sim_dir(self, parameter, md_dir_name):
        src_folder = self.working_dir + "./scripts/gromacs"
        dst_folder = self.working_dir + f"/{md_dir_name}"

        if os.path.exists(dst_folder) and os.path.isdir(dst_folder):
            shutil.rmtree(dst_folder)

        self.make_result_dir(md_dir_name)

        shutil.copytree(src_folder + "/amber14sb_OL15.ff", dst_folder + "/amber14sb_OL15.ff")
        shutil.copytree(src_folder + "/mdp", dst_folder + "/mdp")
        shutil.copy(src_folder + "/single_run.sh", dst_folder + "/single_run.sh")

    def prepare_simulation(self, simulation_parameter, md_dir_name):
        self.make_result_dir(self.md_parameter["simulation_name"])
        self.copy_files_to_sim_dir(simulation_parameter, md_dir_name)
        # Ändern der Parameter in den mdp files
        self.change_tamperature_in_nvt(simulation_parameter["temperature[°C]"], md_dir_name)
        self.change_tamperature_in_npt(simulation_parameter["temperature[°C]"], md_dir_name)
        self.change_sim_time_in_md0(simulation_parameter["simulation_time[ns]"], md_dir_name)

    def copy_input_model(self, path_to_model_pdb, dir):
        """
        Compieng the modeling result structure to the MD simulation directory. File is renamed to input.pdb

        :param path_to_model_pdb: Path to the modeling result structure
        :param dir: Path to the MD run directory
        :return:
        """
        self.make_result_dir(self.md_parameter["simulation_name"])
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

    def run_simulation_steps(self, structureFile, dir, working_dir_path):
        """
        Running bash commands with python subprocess to make a single MD run with GROMACS. Reference single_run.sh

        :param structureFile: Name of input structure file
        :param dir: Path of MD run directory
        :return: none
        """
        os.chdir(working_dir_path + f"/{dir}")
        print(os.getcwd())
        structureName = structureFile[:-4]

        self.make_result_dir("nvt")
        self.run_command(
            f"gmx grompp -f mdp/nvt.mdp -c em/{structureName}.gro -r em/{structureName}.gro -p {structureName}.top -o nvt/{structureName}.tpr -po nvt/{structureName}.mdp -maxwarn 2")
        self.run_command(
            f"gmx mdrun -v -s nvt/{structureName}.tpr -c nvt/{structureName}.gro -x nvt/{structureName}.xtc -cpo nvt/{structureName}.cpt -e nvt/{structureName}.edr -g nvt/{structureName}.log")

        self.make_result_dir("npt")
        self.run_command(
            f"gmx grompp -f mdp/npt.mdp -c nvt/{structureName}.gro -r nvt/{structureName}.gro -t nvt/{structureName}.cpt -p {structureName}.top -o npt/{structureName}.tpr -po npt/{structureName}.mdp -maxwarn 2")
        self.run_command(
            f"gmx mdrun -v -s npt/{structureName}.tpr -c npt/{structureName}.gro -x npt/{structureName}.xtc -cpo npt/{structureName}.cpt -e npt/{structureName}.edr -g npt/{structureName}.log")

        self.make_result_dir("md0")
        self.run_command(
            f"gmx grompp -f mdp/md0.mdp -c npt/{structureName}.gro -t npt/{structureName}.cpt -p {structureName}.top -o md0/{structureName}.tpr -po md0/{structureName}.mdp  -maxwarn 2")
        self.run_command(
            f"gmx mdrun -v -s md0/{structureName}.tpr -c md0/{structureName}.gro -x md0/{structureName}.xtc -cpo md0/{structureName}.cpt -e md0/{structureName}.edr -g md0/{structureName}.log")

        os.chdir(working_dir_path)

    def reduce_center_xtc(self, md_dir):
        """
        Reduce the trajectory to the RNA and center it in the simulation Box.

        At first a ndx file of the RNA is created. Here are only atom id's written belonging to RNA molecules. Then two bash commands are called by python subprocess. These two commands using gmx trjconv to convert trajectories and to produce a pdb file of the frst state with the given ndx file.

        :param md_dir: Path of the MD run directory
        :return: none
        """
        self.make_ndx_of_rna(f"{md_dir}/md0/input.gro", f"{md_dir}/analysis/Index_Files/RNA.ndx")
        md_dir = f"{os.getcwd()}/{md_dir}"
        self.run_command(
            f"gmx trjconv -f {md_dir}/md0/input.xtc -s {md_dir}/md0/input.tpr -o {md_dir}/md0/input_centered.xtc -n {md_dir}_analysis/Index_Files/RNA.ndx -pbc mol -center")
        self.run_command(
            f"gmx trjconv -f {md_dir}/md0/input.xtc -s {md_dir}/md0/input.tpr -o {md_dir}/md0/input_s1.pdb -n {md_dir}_analysis/Index_Files/RNA.ndx -pbc mol -center -b 1 -e 10")


if __name__ == '__main__':
    simulation_parameter = {
        "simulation_name": "Hairpin_labeled",
        "c_magnesium_ions[mol/l]": 0.00,
        "simulation_time[ns]": 0.3,
        "temperature[°C]": 25,
        "dist_to_box[nm]": "1",
    }
    hairpin_labeled = MDSimulation(f"{os.getcwd()}/data",
                                   f"{os.getcwd()}/data/rosetta_results/silent_out.pdb",
                                   simulation_parameter)

