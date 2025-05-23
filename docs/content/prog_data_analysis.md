# Data Analysis

This section describes how to analyze MD simulation results and calculate photophysical properties related to FRET using the FAMP pipeline.

All analysis functions are encapsulated in the `DataAnalysis` class.

---

## Initialization
To begin analysis, initialize the `DataAnalysis` class with the following attributes:

- `working_dir`: Path to the output directory for analysis results
- `path_sim_results`: Path to the directory containing MD output files (`.gro`, `.tpr`, `.xtc`). These file names must match the `simulation_name` used in the parameters.
- `analysis_parameter`: Dictionary with simulation and dye information (see [parameter section](parameter.md))
- `macv_label_pars`: Dictionary with dye parameters for mACV calculation (see [parameter section](parameter.md))

```python
md_analysis = famp.data_analysis.DataAnalysis(
    working_dir="/home/user/famp_project",
    path_sim_results="/home/user/famp_project/simulation1",
    analysis_parameter=analysis_paras,
    macv_label_pars=dye_acv_parameter
)
```


## Step 1: Set Up the Analysis Environment
Creates the required directory structure and initializes input files.

```python
md_analysis.make_data_analysis_results_dirs()
```

---

## Step 2: Calculate FRET Parameters from Explicit Dye Trajectories
This step extracts the inter-dye distances and orientation factors (κ²) from MD trajectories from explicit dyes.

```python
md_analysis.generate_r_kappa_from_dyes()
```

Output: `rkappa.dat` file in the directory `analysis/explicit_dyes/`

---

## Step 3: Calculate FRET Parameters from mACV Model
If no explicit dye simulation is available, this step calculates FRET observables using the multiple accessible contact volume (mACV) model.

```python
md_analysis.generate_r_kappa_from_macv()
```

Output: `rkappa.dat` file in `analysis/MACV/` and additional summary data in `macv.pkl`

---

## Step 4: Simulate FRET Distributions
Using the `rkappa.dat` files, FRET distributions and anisotropy decay curves can be simulated with FRETraj.

### Example
```python
import fretraj as ft

macv_bursts = ft.burst.Experiment(
    path="analysis/MACV/",
    parameter=burst_parameter,
    compute_anisotropy=False,
    units="nm"
)
```

For explicit dye data:
```python
dye_bursts = ft.burst.Experiment(
    path="analysis/explicit_dyes/",
    parameter=burst_parameter,
    compute_anisotropy=True,
    units="nm"
)
```

Description of the `burst_paramter` can be found <a href="https://rna-fretools.github.io/fretraj/background/parameter_file.html#burst-simulation" target="_blank">here</a>

---

## Output Overview
The results of the analysis include:
- `rkappa.dat`: Inter-dye distance and κ² over time
- `macv.pkl`: Serialized FRET trajectory from mACV
- `dye_bursts.FRETefficiencies`: Simulated FRET efficiencies
- `dye_bursts.anisotropy`: Anisotropy decay (donor and acceptor channels)

These results can be visualized or directly compared to experimental FRET data.

---