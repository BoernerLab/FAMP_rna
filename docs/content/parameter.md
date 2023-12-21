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


(content:parameter:simulation)=
## MD Simulaiton Parameter

```
simulation_parameter = {
    "simulation_name": "KLTL ",
    "c_magnesium_ions[mol/l]": 0.02,
    "simulation_time[ns]": 1000,
    "temperature[°C]": 25,
    "dist_to_box[nm]": "1.25",
}
```


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

(content:parameter:burst)=
## Burst Parameter
```
burst_parameter = {
    "dyes": {
        "tauD": 1.14,
        "tauA": 1.5,
        "QD": 0.2,
        "QA": 0.35,
        "etaD": 0.37,
        "etaA": 1,
        "dipole_angle_abs_em": 10.5
    },
    "sampling": {
        "nbursts": 10000,
        "skipframesatstart": 0,
        "skipframesatend": 0,
        "multiprocessing": true
    },
    "fret": {
        "R0": 54,
        "kappasquare": 0.6666,
        "gamma": true,
        "quenching_radius": 10
    },
    "species": {
        "name": ["all"],
        "unix_pattern_rkappa": ["*.dat"],
        "unix_pattern_don_coords": ["Acceptor*.txt"],
        "unix_pattern_acc_coords": ["Donor*.txt"],
        "probability": [1],
        "n_trajectory_splits": null
    },
    "bursts": {
        "lower_limit": 10,
        "upper_limit": 50,
        "lambda": -2.3,
        "QY_correction": false,
        "averaging": "all",
        "burst_size_file": None
    }
}
```