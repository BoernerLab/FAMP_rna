# Simulation

Movement of the RNA and dyes is simulated in the pipeline using MD simulations. 
The GROMACS program is addressed here. MD simulations can be executed with a simple parameter definition. 
The implementation of restraints between atoms is also implemented. 

The functions are implemented in the `MDSimulation` class. An object is defined by calling the class and defining the 
following class attributes.
- `working_dir`: Folder in which the results are to be saved
- `file_path_input`: Path to the input structure (*.pdb*) for the MD simulation
- `md_parameter`: Python dictionary with the defined parameters

Example Code for initialising a MD simulation class: 

```python
simulation = famp.MDSimulation(working_dir=f"{os.getcwd()}/data/simulation",
                               file_path_input=f"{os.getcwd()}/data/input.pdb",
                               md_parameter=simulation_parameter)
```

First, the MD simulation must be prepared in the specified working directory. The following command is used for this: 

```
simulation.prepare_new_md_run()
```
Here the input structure and required files such as force field and parameter files are copied into the working directory. Parameters are also updated. 

The update_parameter() function can be used to update parameters. This reads the parameters from the python dictionary and adjusts them in the GROMACS parameters. 

## Solvation 

In the first step of the MD simulation with GROMACS, the molecule is virtually dissolved in water. The individual steps for the solvation of the molecule can be read here. Command for the molecule solvation with energy minimization. 

```
simulation.solvate_molecule()
```
The results are stored in the working directory under the folder *em*. 

## MD run

With the command: 
```
simulation.run_simulation_steps()
```
The simulation is started. The pressure and temperature are each set to 300 ps and then the full MD run is run with the set simulation time. A detailed description can be found here. 

The results of the MD run can be found in the *md0* folder.

## Restraints

Restraints can be used to force distances between atoms in the simulation. This is useful, for example, to obtain bonds or to induce bonds. 

Restraints are defined in the pipeline using a python dictionary. 
```
restraint_1 = {
        "atom_id_1": 250,
        "atom_id_2": 1831,
        "lower_distance_limit": 0.0,
        "atoms_distance": 0.32,
        "upper_distance_limit": 0.5,
        "force_constant_fraction": 0.5,

    }
```
- `atom_id_1`: Atom number of the first atom for the restraint. The number should be taken from the *.gro* file. 
- `atom_id_2`: Analog to atom 1
- `lower_distance_limit`: Distance for upper limit of free restraint
- `atoms_distance`: Distance to be kept between the atoms
- `upper_distance_limit`: Distance from which the harmonic potential is linearly defined
- `force_constant_fraction`: Divisor by which the force constant is to be multiplied. 

Restraints are defined as GROMACS distance limits using a harmonic potential. Further information can be found <a href="https://manual.gromacs.org/2023.3/reference-manual/functions/restraints.html#distance-restraints" target="_balnk">here</a>

Restraints are collected in a list. Each defined restraint must be added with: 
```
simulation.add_restraint(restraint_1)
```

To activate the restraints, the following command must be executed before the MD run: 
```
simulation.apply_restraints()
```

## Parameter
```python
simulation_parameter = {
    "simulation_name": "simulation1",
    "c_magnesium_ions[mol/l]": 0.02,
    "simulation_time[ns]": 1,
    "temperature[°C]": 25,
    "dist_to_box[nm]": "1",
    "water_model": "tip3p",
    "distance_restraints": True
}
```

The parameters for the simulation are a reduced selection of Rosetta flag parameters, which are required for the tertiary structure modeling.

- *simulation_name*: Name of the current simulation (Name for directory)
- *c_magnesium_ions[mol/l]*: Concentration of divalent magnesium iones in $\frac{mol}{l}$ used in the simulaiton. Divalent ions are limeted to Mg(II)
- *simulation_time[ns]*: The total time the simulation should run in ns.
- *temperature[°C]*: Temperature used in the simulation. 
- *dist_to_box[nm]*: Distance of the molecule to the edges of the box. Takes influcence of the box size.
- *water_model*: *tip3p* or *tip4p*. It is possible to use the 3-Point or the 4-point water model in the simulation.  
- *distance_restraints*: Parameter for activating distance restraints in the simulation. 