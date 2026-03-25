import subprocess
import re
import os
import shutil
import pathlib
import numpy as np
import mdtraj as md
import fretraj as ft
from famp.exceptions import GromacsMdrunError, GromacsEditconfError, GromacsGenionError, GromacsGromppError, \
    GromacsSolvateError, GromacsPdb2gmxError


class MDSimulation:
    def __init__(self, working_dir: str, file_path_input: str, md_parameter: dict) -> None:
        self.source_path = pathlib.Path(__file__).parent.resolve()
        self.working_dir = self.create_working_dir(working_dir)
        self.file_path_input = file_path_input
        self.md_parameter = md_parameter
        self.path_simulation_folder = self.get_simulation_path()
        self.input_structure_name = self.get_input_structure_name()
        self.temp_pseudo_restraint_data = []
        self.restraints = []

    @staticmethod
    def run_gromacs_command(command: str):
        """Uses subprocess to run gromacs commands.

        This function executes subprocesses by commands. In particular, GROMACS bash commands are executed here and
        their errors are checked. Specific exceptions for different modules are thrown for different GROMACS modules.
        Command stdout is printed when the process is finished.

        :param command: Simple string with a bash command
        """
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
    def run_mdrun_command(command: str):
        """Uses subprocess to run gromacs mdrun with live output.

        Utilized solely for executing the mdrun command, displaying real-time output in the terminal.
        Any command errors are not captured.
        :param command: Simple string with a bash command
        """
        process = subprocess.Popen(["bash", "-c", command], stdout=subprocess.PIPE, universal_newlines=True)

        process.wait()

    @staticmethod
    def sim_time_to_steps(sim_time: float) -> int:
        """Converting simulation time to simulation steps.

        The function calculates the number of simulation steps needed based on a given total simulation time in
        nanoseconds for an integration time step of 2 fs.

        :param sim_time:

        """
        return int(1000000 * sim_time / 2)

    @staticmethod
    def degree_to_kelvin(degree: int) -> int:
        """Converts °C to Kelvin

        :param degree: °C

        """
        return degree + 273

    @staticmethod
    def create_working_dir(working_dir: str) -> str:
        """Creates a (working) directory

        Creates a directory if it does not already exist. If the directory exists, only a print is output.

        :param working_dir: Path of the directory to be created

        """
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

    def get_simulation_path(self) -> str:
        return f"{self.working_dir}/{self.md_parameter['simulation_name']}"

    def get_input_structure_name(self) -> str:
        file_name = os.path.basename(os.path.normpath(self.file_path_input))
        structure_name = file_name[:-4]
        return structure_name

    #______________ Pseudoatom Section__________________________________________________________________________________

    def prepare_structure_with_pseudoatoms(self, acv_parameter_file: str):
        """
        Calculates ACVs based on labels, adds Pseudo-Atoms (PS) to the structure,
        and saves a new PDB. It also stores restraint data for later remapping.

        :param acv_parameter_file: Path to the JSON file containing label definitions.
        """
        print("--- Preparing Pseudo Atoms ---")

        labels = ft.cloud.labeling_params(acv_parameter_file, verbose=False)

        structure = md.load(self.file_path_input)

        volumes = []
        if 'Position' in labels:
            for dye_name in labels['Position'].keys():
                print(f"Calculating Volume for {dye_name}...")
                vol = ft.cloud.Volume(structure, dye_name, labels)
                volumes.append(vol)

        if not volumes:
            print("No label positions found in JSON.")
            return

        new_xyz = structure.xyz.copy()
        top = structure.top

        # get last chain
        last_chain = top.chain(top.n_chains - 1)
        last_residue = last_chain.residue(last_chain.n_residues - 1)
        current_resSeq = last_residue.resSeq + 1

        self.temp_pseudo_restraint_data = []  # Reset data

        for volume in volumes:
            # add new residue in pdb
            res = top.add_residue("PS", last_chain, resSeq=current_resSeq)
            current_resSeq += 1

            # add new atom in pdb
            top.add_atom("PS", md.element.bromine, res)

            # calcaulate pseudo atom coords
            pseudo_pos = np.array(volume.acv.mp, dtype=float) / 10.0

            if pseudo_pos.shape != (3,):
                raise ValueError(f"Pseudoatom position shape error: {pseudo_pos.shape}")

            new_xyz = np.concatenate([new_xyz, pseudo_pos.reshape(1, 1, 3)], axis=1)

            # Restraints berechnen (Geometrie vor pdb2gmx)
            res_data = self.calculate_distance_restraints_geometry(volume)
            self.temp_pseudo_restraint_data.append(res_data)

        # 5. save new structure
        new_filename = "pseudo_atoms_added.pdb"
        new_path = os.path.join(self.working_dir, new_filename)

        new_structure = md.Trajectory(new_xyz, top)
        new_structure.save_pdb(new_path)

        print(f"New structure with pseudo atoms saved to: {new_path}")

        # 6. Update Input File Path
        self.file_path_input = new_path
        self.input_structure_name = "pseudo_atoms_added"

    def calculate_distance_restraints_geometry(self, volume, selection="name P", cutoff=15):
        """Internal helper to calculate distances and mapping info."""
        current_structure = volume.structure
        top = current_structure.top
        attach_id_mdtraj = volume.attach_id - 1

        # calculate neighbor atoms
        neighbor_idx = md.compute_neighbors(current_structure, cutoff / 10, [attach_id_mdtraj])[0]
        atomsele_idx = top.select(selection)
        atoms_idx = np.array([idx for idx in neighbor_idx if idx in atomsele_idx])

        if len(atoms_idx) != 0:
            mapping_info = []
            for idx in atoms_idx:
                atom = top.atom(idx)
                mapping_info.append({
                    "resSeq": atom.residue.resSeq,
                    "resName": atom.residue.name,
                    "name": atom.name
                })

            atoms_xyz = current_structure.xyz[0, atoms_idx, :] * 10  # nm to Angstrom calc logic from notebook
            mp_angstrom = volume.acv.mp

            # Distance in nm
            distances = np.sqrt(np.sum((mp_angstrom - atoms_xyz) ** 2, axis=1)) / 10

            return {
                "attach_atoms": {
                    "xyz": current_structure.xyz[0, atoms_idx, :],
                    "mapping": mapping_info
                },
                "distances": distances,
            }
        else:
            print("No selected atom types within cutoff for restraints.")
            return None

    def create_pseudo_restraints(self, gro_file_path: str, itp_filename="distance_restraints.itp"):
        """
        Called after pdb2gmx. Maps original atoms to new GROMACS indices,
        identifies pseudo atoms, and writes the .itp file.
        """
        if not self.temp_pseudo_restraint_data:
            print("No pseudo atom restraint data found. Skipping .itp generation.")
            return

        print(f"Mapping restraints to {gro_file_path}...")

        traj_gro = md.load(gro_file_path)
        top_gro = traj_gro.top

        # 1. create lookup table
        lookup_map = {}
        for atom in top_gro.atoms:
            key = (atom.residue.resSeq, atom.name)
            lookup_map[key] = atom.index + 1  # GROMACS is 1-based

        # 2. find pseudo atoms (PS)
        ps_indices = top_gro.select("resname PS")
        if len(ps_indices) == 0:
            print("Warning: No residues named 'PS' found in .gro file.")
            return

        pseudo_atom_ids_gro = [idx + 1 for idx in ps_indices]

        # 3. remap data
        remapped_restraints = []
        for group in self.temp_pseudo_restraint_data:
            if group is None: continue

            new_group = group.copy()
            new_serials = []

            for atom_info in group['attach_atoms']['mapping']:
                key = (atom_info['resSeq'], atom_info['name'])
                if key in lookup_map:
                    new_serials.append(lookup_map[key])
                else:
                    raise ValueError(f"Atom {atom_info} not found in {gro_file_path} after pdb2gmx.")

            new_group['attach_atoms']['serial'] = np.array(new_serials)
            remapped_restraints.append(new_group)

        # 4. write .itp file
        itp_path = os.path.join(self.path_simulation_folder, "em", itp_filename)  # Save in em folder first

        margin = 0.3  # nm deviation
        content = [
            "; GROMACS distance restraints for pseudo atoms",
            "[ distance_restraints ]",
            ";  ai    aj   type   index   type'      low     up1     up2     fac"
        ]

        counter = 0
        for group_data, p_id in zip(remapped_restraints, pseudo_atom_ids_gro):
            atom_id_1 = int(p_id)
            attached = group_data['attach_atoms']['serial']
            dists = group_data['distances']

            for atom_id_2, d_val_np in zip(attached, dists):
                d_val = float(d_val_np)
                lower = max(0.0, d_val - margin)
                upper = d_val + margin

                line = f"{atom_id_1:6d} {int(atom_id_2):6d} 1 {counter:4d} 1 {lower:.4f} {d_val:.4f} {upper:.4f} 1.0"
                content.append(line)
                counter += 1

        with open(itp_path, "w") as f:
            f.write("\n".join(content) + "\n")

        print(f"Generated {itp_filename} with {counter} restraints.")

        # 5. In Topology einbinden
        top_file = os.path.join(self.path_simulation_folder, f"{self.input_structure_name}.top")
        self._include_itp_in_top(top_file, itp_filename)

    def _include_itp_in_top(self, top_file, itp_file):
        """Injects #include statement into topol.top"""
        with open(top_file, 'r') as f:
            lines = f.readlines()

        include_line = f'#include "{itp_file}"\n'
        if any(include_line.strip() in l for l in lines):
            return

            # Insert before [ molecules ] or after forcefield includes
        # Strategy: Find [ molecules ] and insert before
        insert_idx = -1
        for i, line in enumerate(lines):
            if "[ molecules ]" in line:
                insert_idx = i
                break

        if insert_idx != -1:
            lines.insert(insert_idx, f"\n{include_line}\n")
            with open(top_file, 'w') as f:
                f.writelines(lines)
            print(f"Included {itp_file} in topology.")

    #___________________________________________________________________________________________________________________

    def change_temperature_in_mdp_files(self, temperature: int):
        """Changes the tempreature in the mdp files

        :param temperature: integer value
        """
        temp_in_k = self.degree_to_kelvin(temperature)
        for source_dir in ["nvt", "npt", "md0"]:
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

    def change_sim_time_in_md0(self, time):
        """Changing the simulation time

        This function changes the simulation steps in a mdp file. The entered simulation time (ns) is converted
        into the simulation time with integration step size of 2 fs with the function sim_time_to_step().

        :param time: Specification of the time to be simulated in ns.
        """
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
        """Copies files to the simulaiton folder.

        A new folder for the simulation run is created here. The standard parameter files (mdp) and the force field are
        copied from the scripts folder to the newly created directory.

        """
        src_folder = f"{self.source_path}/scripts/gromacs"
        dst_folder = self.working_dir + f"/{self.md_parameter['simulation_name']}"

        if os.path.exists(dst_folder) and os.path.isdir(dst_folder):
            print(
                "MD run already exists. To make a new Simulation change the Name or the Simulation in the MD parameter")
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

    def copy_input_model(self) -> None:
        """
        Copies the modeling result structure to the MD simulation directory.
        """
        self.make_result_dir(self.md_parameter["simulation_name"])
        input_file_name = os.path.basename(os.path.normpath(self.file_path_input))
        shutil.copy(f"{self.file_path_input}",
                    f"{self.working_dir}/{self.md_parameter['simulation_name']}/{input_file_name}")

    def prepare_new_md_run(self):
        """Prepares a new MD run.

        A folder for the simulation run is created and all necessary files like parameters, input structure and force
        field are added to the directory. Then the temperature and the simulation time, defined in the parameters,
        are updated within the parameter files.

        """

        self.copy_files_to_sim_dir()
        self.copy_input_model()
        self.change_temperature_in_mdp_files(self.md_parameter["temperature[C]"])
        self.change_sim_time_in_md0(self.md_parameter["simulation_time[ns]"])

    def update_parameter(self):
        """
        Function to update the parameters simulation time and temperature within the mdp files.
        """
        self.change_temperature_in_mdp_files(self.md_parameter["temperature[C]"])
        self.change_sim_time_in_md0(self.md_parameter["simulation_time[ns]"])

    def solvate_molecule(self) -> None:
        """
        Running bash commands with python subprocess to solvate MD run with GROMACS. Reference solvate.sh.
        Distinguishes between tip3p and tip4p water. Always uses the amber14sb_OL15 force field for the simulation,
        which must be located in the installation directory of GROMACS. With the dockerization of GROMACS this
        will be updated in the future.


        :return: none
        """

        if self.md_parameter["water_model"] == "tip3p":
            water_file = "spc216.gro"
        elif self.md_parameter["water_model"] == "tip4p":
            water_file = "tip4p.gro"
        else:
            raise ValueError(
                "Please use tip3p or tip4p as water model. Other water models are "
                "currently not implemented within the pipeline.")

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

        gro_file = f"{self.path_simulation_folder}/em/{self.input_structure_name}.gro"
        if self.temp_pseudo_restraint_data:
            self.create_pseudo_restraints(gro_file, itp_filename="distance_restraints.itp")

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

        self.run_gromacs_command(
            f"gmx select "
            f"-f {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-on {self.path_simulation_folder}/em/SOL_K.ndx "
            f"-select 'SOL'")

        self.run_gromacs_command(
            f"gmx genion "
            f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-nname Cl "
            f"-pname K "
            f"-neutral "
            f"-n {self.path_simulation_folder}/em/SOL_K.ndx")

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/em.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        if self.md_parameter["c_magnesium_ions[mol/l]"] > 0:
            self.run_gromacs_command(
                f"gmx select "
                f"-f {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
                f"-on {self.path_simulation_folder}/em/SOL_MG.ndx "
                f"-select 'SOL'")

            self.run_gromacs_command(
                f"gmx genion "
                f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
                f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
                f"-nname Cl "
                f"-pname MG "
                f"-pq 2 "
                f"-conc {self.md_parameter['c_magnesium_ions[mol/l]']} "
                f"-n {self.path_simulation_folder}/em/SOL_MG.ndx")

            self.run_gromacs_command(
                f"gmx grompp "
                f"-f {self.path_simulation_folder}/mdp/em.mdp "
                f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
                f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
                f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
                f"-po {self.path_simulation_folder}/em/{self.input_structure_name}.mdp "
                f"-maxwarn 2")

        self.run_mdrun_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/em/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-o {self.path_simulation_folder}/em/{self.input_structure_name}.trr "
            f"-e {self.path_simulation_folder}/em/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/em/{self.input_structure_name}.log")

    def run_simulation_steps(self) -> None:
        """
        Running bash commands with python subprocess to make a single MD run with GROMACS. Reference single_run.sh
        If the tip3p water model is selected, the calculation is preferred to the GPU, as this improves the performance
        of the simulation.

        :return: None
        """
        self.make_result_dir(f"{self.md_parameter['simulation_name']}/nvt")

        if self.md_parameter["water_model"] == "tip3p":
            update_parameter = "gpu"
        else:
            update_parameter = "auto"

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/nvt.mdp "
            f"-c {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-r {self.path_simulation_folder}/em/{self.input_structure_name}.gro "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/nvt/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/nvt/{self.input_structure_name}.mdp "
            f"-maxwarn 2")

        self.run_mdrun_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/nvt/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/nvt/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/nvt/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/nvt/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/nvt/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/nvt/{self.input_structure_name}.log "
            f"-update {update_parameter} "
            f"-nt {self.md_parameter['number_of_cores']}")

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

        self.run_mdrun_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/npt/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/npt/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/npt/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/npt/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/npt/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/npt/{self.input_structure_name}.log "
            f"-update {update_parameter} "
            f"-nt {self.md_parameter['number_of_cores']}")

        self.make_result_dir(f"{self.md_parameter['simulation_name']}/md0")

        self.apply_restraints()

        self.run_gromacs_command(
            f"gmx grompp "
            f"-f {self.path_simulation_folder}/mdp/md0.mdp "
            f"-c {self.path_simulation_folder}/npt/{self.input_structure_name}.gro "
            f"-t {self.path_simulation_folder}/npt/{self.input_structure_name}.cpt "
            f"-p {self.path_simulation_folder}/{self.input_structure_name}.top "
            f"-o {self.path_simulation_folder}/md0/{self.input_structure_name}.tpr "
            f"-po {self.path_simulation_folder}/md0/{self.input_structure_name}.mdp  "
            f"-maxwarn 2")
        self.run_mdrun_command(
            f"gmx mdrun -v "
            f"-s {self.path_simulation_folder}/md0/{self.input_structure_name}.tpr "
            f"-c {self.path_simulation_folder}/md0/{self.input_structure_name}.gro "
            f"-x {self.path_simulation_folder}/md0/{self.input_structure_name}.xtc "
            f"-cpo {self.path_simulation_folder}/md0/{self.input_structure_name}.cpt "
            f"-e {self.path_simulation_folder}/md0/{self.input_structure_name}.edr "
            f"-g {self.path_simulation_folder}/md0/{self.input_structure_name}.log "
            f"-update {update_parameter} "
            f"-nt {self.md_parameter['number_of_cores']}")

        # os.chdir(working_dir_path)

    def add_restraint(self, restraint: dict):
        """Adds a restraint to the restraint attribute

        :param restraint:dictionary with the parameters describing a restraint.
        """
        self.restraints.append(restraint)

    def remove_all_restraints(self):
        """initializes the restraint attribute

        """
        self.restraints = []

    def generate_text_for_topology(self) -> str:
        """Generates the text to be added to the .top file to define restraints.
        The header for the restraints is created. The restraint attribute (type: list) is then iterated through to
        list the restraints in a formatted way.

        :return: Formatted text for the topology file. Type:str
        """
        formatted_text = f"\n[ distance_restraints ] \n; ai   aj   type   index   type'      low     up1     up2     fac\n"

        for i, restraint in enumerate(self.restraints):
            formatted_text += f"{restraint['atom_id_1']}\t{restraint['atom_id_2']}\t1\t{i}\t1\t{restraint['lower_distance_limit']}\t{restraint['atoms_distance']}\t{restraint['upper_distance_limit']}\t{restraint['force_constant_fraction']}\n"

        return formatted_text

    def add_restraints_to_topology(self):
        """Adds the formatted text that defines the restraints to the topology file.

        If a section starting with [distance restraints] exists, it is removed,
        including everything up to "; Include water topology".
        The string ("; Include water topology") must exist in the file.
        The new restraints are added directly before "; Include water topology".
        """
        file_content = []
        start_restraints = None
        end_restraints = None
        water_topology_position = None

        # Read the topology file and identify the sections to modify
        with open(f"{self.path_simulation_folder}/{self.input_structure_name}.top", 'r') as f:
            for i, line in enumerate(f):
                stripped_line = line.strip()
                if "[distance restraints]" in stripped_line and start_restraints is None:
                    start_restraints = i
                if "; Include water topology" in stripped_line:
                    water_topology_position = i
                    if start_restraints is not None and end_restraints is None:
                        end_restraints = i
                file_content.append(line.rstrip())  # Strip trailing whitespace to avoid empty lines

        # Remove existing [distance restraints] section if it exists
        if start_restraints is not None and end_restraints is not None:
            del file_content[start_restraints:end_restraints]

        # Generate the new restraint text
        restraint_text = self.generate_text_for_topology().split("\n")
        restraint_text = [line.strip() for line in restraint_text if line.strip()]  # Remove empty lines from new text

        # Insert the new restraint text before "; Include water topology"
        if water_topology_position is not None:
            for i, text in enumerate(reversed(restraint_text)):
                file_content.insert(water_topology_position, text)

        # Write the updated content back to the file
        with open(f"{self.path_simulation_folder}/{self.input_structure_name}.top", 'w') as f:
            f.write("\n".join(file_content) + "\n")

    def activate_restraints_mdp(self):
        """Activates the parameter for distance restraints.

        The function activates the parameter "disre = simple" in the mdp file md0. This causes the restraints to be applied in
        the simulation long run. The parameter is added if the parameter does not exist in the md0 file.


        """

        content = []
        with open(f"{self.path_simulation_folder}/mdp/md0.mdp", 'r') as f:
            if "disre" in f.read():
                for i, line in enumerate(f):
                    line = line.strip()
                    if line.startswith('disre'):
                        if line.contains("simple"):
                            content.append(line)
                        elif line.contains("no"):
                            content.append("disre\t= simple")
                    else:
                        content.append(line)
            else:
                f.seek(0)
                for i, line in enumerate(f):
                    line = line.strip()
                    content.append(line)
                content.append("disre\t= simple")


        with open(f"{self.path_simulation_folder}/mdp/md0.mdp", 'w') as f:
            for line in content:
                f.write("%s\n" % line)

    def deactivate_restraints_mdp(self):
        """Disables distnace retraints in the md0 file.

        If the parameter "disre" exists in the md0 file it will be set to "disre = no".

        """
        content = []
        with open(f"{self.path_simulation_folder}/mdp/md0.mdp", 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if line.startswith('disre'):
                    if line.contains("simple"):
                        content.append("disre\t= no")
                    elif line.contains("no"):
                        content.append(line)
                else:
                    content.append(line)

        # print(content)
        with open(f"{self.path_simulation_folder}/mdp/md0.mdp", 'w') as f:
            for line in content:
                f.write("%s\n" % line)

    def apply_restraints(self):
        """
        Allows the application of restraints.

        The function checks the simulation parameters whether distance restraints should be activated or
        deactivated and executes the corresponding functions.
        """
        if self.md_parameter['distance_restraints']:
            if not self.restraints:
                print('There are no restraints to add. Restraints are deactivated')
                self.deactivate_restraints_mdp()
            else:
                self.add_restraints_to_topology()
                self.activate_restraints_mdp()
        else:
            self.deactivate_restraints_mdp()


if __name__ == '__main__':
    simulation_parameter = {
        "simulation_name": "lukas_poly_u_test_2021er_roteded_restraint",
        "c_magnesium_ions[mol/l]": 0.02,
        "simulation_time[ns]": 0.2,
        "temperature[C]": 25,
        "dist_to_box[nm]": "1",
        "water_model": "tip4p",
        "distance_restraints": True,
        "number_of_cores": "24"
    }
    print(os.getcwd())
    hairpin_labeled = MDSimulation(working_dir=f"/home/felix/Documents/md_lukas_test_2021er",
                                   file_path_input=f"/home/felix/Documents/md_lukas_test_2021er/lukas_polyu_rotated.pdb",
                                   md_parameter=simulation_parameter)

    hairpin_labeled.prepare_new_md_run()
    hairpin_labeled.update_parameter()
    hairpin_labeled.solvate_molecule()

    restraint_1 = {
        "atom_id_1": 193,
        "atom_id_2": 1791,
        "lower_distance_limit": 0.0,
        "atoms_distance": 0.32,
        "upper_distance_limit": 0.52,
        "force_constant_fraction": 1.0,

    }

    restraint_2 = {
        "atom_id_1": 181,
        "atom_id_2": 973,
        "lower_distance_limit": 0.0,
        "atoms_distance": 0.37,
        "upper_distance_limit": 0.57,
        "force_constant_fraction": 1.0,

    }

    restraint_3 = {
        "atom_id_1": 214,
        "atom_id_2": 973,
        "lower_distance_limit": 0.0,
        "atoms_distance": 0.34,
        "upper_distance_limit": 0.54,
        "force_constant_fraction": 1.0,

    }

    restraint_4 = {
        "atom_id_1": 220,
        "atom_id_2": 1758,
        "lower_distance_limit": 0.0,
        "atoms_distance": 0.31,
        "upper_distance_limit": 0.51,
        "force_constant_fraction": 1.0,

    }

    restraint_5 = {
        "atom_id_1": 226,
        "atom_id_2": 1758,
        "lower_distance_limit": 0.0,
        "atoms_distance": 0.27,
        "upper_distance_limit": 0.47,
        "force_constant_fraction": 1.0,

    }

    restraint_6 = {
        "atom_id_1": 259,
        "atom_id_2": 1725,
        "lower_distance_limit": 0.0,
        "atoms_distance": 0.31,
        "upper_distance_limit": 0.51,
        "force_constant_fraction": 1.0,

    }

    hairpin_labeled.add_restraint(restraint_1)
    hairpin_labeled.add_restraint(restraint_2)
    hairpin_labeled.add_restraint(restraint_3)
    hairpin_labeled.add_restraint(restraint_4)
    hairpin_labeled.add_restraint(restraint_5)
    hairpin_labeled.add_restraint(restraint_6)
    #hairpin_labeled.apply_restraints()
    hairpin_labeled.run_simulation_steps()
