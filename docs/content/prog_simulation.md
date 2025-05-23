# Simulation

This section describes how to run molecular dynamics (MD) simulations with the FAMP pipeline using **GROMACS**. The goal is to simulate the motion of RNA molecules and attached dyes under realistic conditions.

All simulation functions are encapsulated in the `MDSimulation` class.

---

## Initialization
To begin an MD simulation, you need to initialize the `MDSimulation` class with the following attributes:

- `working_dir`: Directory where the simulation results will be saved
- `file_path_input`: Path to the input structure file in PDB format
- `md_parameter`: Dictionary containing MD-specific parameters (see [parameter section](parameter.md))

```python
import famp
import os

simulation = famp.simulation.MDSimulation(
    working_dir=f"{os.getcwd()}/data/simulation",
    file_path_input=f"{os.getcwd()}/data/input.pdb",
    md_parameter=simulation_parameter
)
```

---

## Step 1: Preparation
To prepare the MD run, call:
```python
simulation.prepare_new_md_run()
```
This function:
- sets up the directory structure,
- copies the force field and parameter files,
- imports the RNA 3D structure,
- updates simulation parameters based on the provided dictionary.

The function `update_parameter()` can be used to update specific values from the parameter set programmatically.

---

## Step 2: Solvation
Sets up the virtual water box and solvates the RNA 3D structure. This process includes:
- generating a simulation box,
- solvating the RNA with water molecules,
- adding counterions to neutralize the system,
- performing energy minimization.


```python
simulation.solvate_molecule()
```

All intermediate files and results are saved in the `em/` subdirectory of the working directory.

---

## Step 3: Equilibration and Production Run
Start the full MD simulation using:
```python
simulation.run_simulation_steps()
```
This step performs:
1. **Temperature equilibration** (NVT ensemble, 300 ps)
2. **Pressure equilibration** (NPT ensemble, 300 ps)
3. **Full MD production run** according to `simulation_time[ns]` defined in the parameter file

All result files (trajectories, logs, structures) are saved in dedicated subdirectories (e.g., `nvt/`, `npt/`, `md0/`).

---

## Notes
- Short test runs (1–2 ns) are recommended before long production runs.
- All necessary files are automatically generated, but advanced users may customize the `.mdp` files manually if needed.
