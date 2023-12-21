# Data Analysis

Here, on the one hand, the Ergennis files from the MD simulation are further processed and, on the other hand, photophysical parameters of the dyes or the MACV are calculated.

The functions are made usable with the initialization of the class. 

The class `DataAnalysis` is initialized with the following attributes:

- `working_dir`: Folder in which the results are to be saved
- `path_sim_results`: Path to the folder in which the MD simulation results are saved. The following files must exist: *.gro*, *.tpr*, *.xtc*. The names of the files must be the same as the *simulation_name* parameter in the analysis parameters. 
- `md_parameter`: Python dictionary with the defined parameters

Example Code for initialising a DataAnalysis class: 

```python
md_analysis = DataAnalysis(working_dir="/home/felix/Documents/md_pipeline_testfolder",
                           path_sim_results="/home/felix/Documents/md_pipeline_testfolder/m_tlr_ub",
                           analysis_parameter=analysis_paras, macv_label_pars=dye_acv_parameter)
```

## The Parameter

The Parameter here are defined to point out the MD results and to define the type of dyes and their attachment points

```python
analysis_paras = {
    "simulation_name": "simulation1",
    "input_structure_name": "input",
    "Donor_residue_name_number": ("C3W", 10),
    "Acceptor_residue_name_number": ("C5W", 45),
}
```
- *simulation_name*: Name of the simulation. Specifies the directory of the results. Should be the same name as in the simulation parameters.
- *input_structure_name*: File name of the input structure of the MD simulation
- *Donor_residue_name_number*: Tuple that defines the force field name of the donor dye and the residue number where the dye is attached.
- *Acceptor_residue_name_number*: Tuple that defines the force field name of the acceptor dye and the residue number where the dye is attached.

## Usage

The first function build the setup for the analysis. Here, directories are created in an *analysis* folder. The xtc file 
is reduced to the RNA and dyes and the periodic boundary conditions are solved. Then all the needed files for the analysis 
where copied from the MD results. 

```python
md_analysis.make_data_analysis_results_dirs()
```

Now the calculation of distances and $\kappa^2$ values from the trajectory of explicit dyes can be performed with: 

```python
md_analysis.generate_r_kappa_from_dyes()
```
