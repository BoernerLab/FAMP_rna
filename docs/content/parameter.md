# Parameter 

The parameters listed here are used in the corresponding modules. The parameters are reduced parameters for programs 
such as GROMACS or Rosetta. They have been reduced to essential aspects for modeling and simulation.

(content:parameter:rosetta)=
## Rosetta Parameter

```
rosetta_parameter = {
    "path_to_rosetta": "rna_denovo.default.macosclangrelease",
    "nstruct": 50,
    "fasta": "rna_tlr_sequence.fasta",
    "minimize_rna": True,
    "cycles": 200
}
```

### Description
- `path_to_rosetta`: Path to the Rosetta executable for RNA 3D structure prediction (FARFAR2).
- `nstruct`: Number of structures to generate in a single modeling run.
- `fasta`: Path to the input RNA sequence in FASTA format.
- `minimize_rna`: Whether to perform structure minimization after generation.
- `cycles`: Number of Monte Carlo sampling cycles used in FARFAR2 modeling.

---


(content:parameter:simulation)=
## MD Simulaiton Parameter

```python
simulation_parameter = {
    "simulation_name": "KLTL",
    "c_magnesium_ions[mol/l]": 0.02,
    "simulation_time[ns]": 1000,
    "temperature[C]": 25,
    "dist_to_box[nm]": 1.25,
    "water_model": "tip3p"
}
```

### Description
- `simulation_name`: Identifier for this MD run (used to name output folders).
- `c_magnesium_ions[mol/l]`: Magnesium ion concentration (Mg²⁺) in mol/l to match experimental conditions.
- `simulation_time[ns]`: Total length of the MD simulation in nanoseconds.
- `temperature[C]`: Temperature at which the simulation is run (in °C).
- `dist_to_box[nm]`: Distance from the solute to the edge of the periodic box (defines box size).
- `water_model`: Specifies water model used in simulation (e.g., `tip3p` or `tip4p`).

---

(content:parameter:analysis)=
## Analysis Prameter
```
analysis_paras = {
        "simulation_name": "KL_TL",
        "input_structure_name": "KL_TL_input",
        "Donor_residue_name_number": ("C3W", 10),
        "Acceptor_residue_name_number": ("C5W", 45),
    }
```

### Description
- `simulation_name`: Must match the MD run identifier.
- `input_structure_name`: File name of the RNA structure used for the MD run.
- `Donor_residue_name_number`: Tuple (residue name, number) of donor dye attachment site.
- `Acceptor_residue_name_number`: Tuple (residue name, number) of acceptor dye attachment site.

---

(content:parameter:dyes)=
## Reduced dye Parameter
```
dye_acv_parameter = {
        "Acceptor": {
            "name": "sCy5",
            "linker_length": 20,
            "linker_width": 3.5,
            "dye_radius1": 9.5,
            "dye_radius2": 3,
            "dye_radius3": 1.5,
            "cv_fraction": 0.99,
            "cv_thickness": 3
        },
        "Donor": {
            "name": "sCy3",
            "linker_length": 20,
            "linker_width": 3.5,
            "dye_radius1": 8.0,
            "dye_radius2": 3,
            "dye_radius3": 1.5,
            "cv_fraction": 0.99,
            "cv_thickness": 3,
        },
        "Distance": {"sCy3-sCy5":
                         {"R0": 54,
                          "n_dist": 10 ** 6}
                     }

    }
```

### Description
- `name`: Identifier for the fluorophore (e.g., `sCy3`, `sCy5`).
- `linker_length`: Length of the dye linker (Å).
- `linker_width`: Width parameter for linker flexibility (Å).
- `dye_radius1/2/3`: Shape parameters for the accessible contact volume.
- `cv_fraction`: Contact volume fraction for stacking interaction modeling.
- `cv_thickness`: Thickness of the contact shell in ACV model.
- `R0`: Förster radius between donor and acceptor (in Å).
- `n_dist`: Number of sampled distances used in burst simulations.

---
