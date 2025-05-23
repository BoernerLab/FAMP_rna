# Modeling

The modeling module within the FAMP pipeline is designed for *de novo* RNA structure prediction. Two external tools are employed:
- **RNAfold** for secondary structure (2D) prediction
- **Rosetta (FARFAR2)** for tertiary structure (3D) prediction

All modeling-related functions are encapsulated in the `Modeling` class. This module forms the first step in the FAMP pipeline and enables generation of structural input models for FRET-based analysis.

## Initialization
To start modeling, initialize the class by providing:
- `working_dir`: Output directory for generated files
- `file_path_sequence`: Path to the input RNA sequence in FASTA format
- `modeling_parameter`: Python dictionary defining the modeling parameters (see the [parameter section](parameter.md))

```python
import famp
import os

current_dir = os.getcwd()
modeling = famp.modeling.Modeling(
    working_dir=f"{current_dir}",
    file_path_sequence=f"{current_dir}/input_data/RNA_Hairpin.fasta",
    modeling_parameter=rosetta_parameter
)
```


The input FASTA file should contain a single RNA sequence in the following format:
```
>Example_Sequence
caauauuuauuaauaucuuccggauauuaauaaauauug
```

---

## Secondary Structure Prediction (2D)

The function `predict_2d_structure()` performs secondary structure prediction using **RNAfold**. The RNA sequence is read from the input FASTA file, and the result is stored as a dot-bracket notation file.

- Output directory: `secondary_prediction/`
- Output file: `dot_bracket.secstruct`

This file can be edited manually after prediction to define or modify secondary structures.

```python
modeling.predict_2d_structure()
```

---

## Tertiary Structure Prediction (3D)

3D structure prediction is performed using Rosetta’s **FARFAR2** module via the `predict_3d_structure()` function. It requires the RNA secondary structure file (from `dot_bracket.secstruct`)

- Output directory: `rosetta_results`
- Output file: `silent_out.out`

```python
modeling.predict_3d_structure("dot_bracket.secstruct")
```

---

## Extracting PDB Structures

Once Rosetta modeling is complete, the `extract_pdb()` function is used to convert the generated structure files into PDB format.

- The models are ranked by Rosetta’s "res4" score.
- The function returns the top `n` structures, where `n` is defined by `number_of_structures`.

```python
modeling.extract_pdb(number_of_structures=5)
```

The resulting `.pdb` files are stored in the output directory and can be used for subsequent MD simulations or in silico dye labeling.

---


