# Modeling

Modeling within the pipeline program is limited to *de novo* modeling. Two external programs are used for this. The RNA 
Fold program is used for secondary structure prediction. Rosetta is used for tertiary structure prediction. 

All functions are implemented in the `Modeling` class. The class is initialized by calling the class and defining the following attributes:
- `working_dir` = directory in which the results are to be saved
- `file_path_sequence` = path to the input fasta sequence
- `modeling_parameter` = Python dictiononary with the defined modeling parameters (see parameters)

Example Code for initialising a modeling class:
```python
import famp
current_dir = os.getcwd()
modeling = famp.modeling.Modeling(working_dir=f"{current_dir}",
                    file_path_sequence=f"{current_dir}/input_data/RNA_Hairpin.fasta",
                    modeling_parameter=rosetta_parameter)
```

The input RNA sequence should be in **fasta** format:

```
> Example_Sequence
caauauuuauuaauaucuuccggauauuaauaaauauug
```

## 2D Prediction

The `predict_2d_structure()` function is used to perform a secondary structure prediction. This uses RNAfold with standard parameters and 
predicts the secondary structure. The sequence is automatically read from the fasta file.
The result is stored in the folder *secondary_prediction* as *dot_bracket.secstruct* 
file. In this file, the secondary structure can be edited manually afterward. 

```python
modeling.predict_2d_structure()
```
## 3D Prediction

The tertiary structure modeling is done with the Rosetta program and is executed with the `predict_3d_structure`
function. The sequence and a secondary structure are required as input. These are located in the *dot_bracket.secstruct* 
file. Parameters for the modeling are defined in a Python dictionary (here). To extract the structure ensemble from the 
modeling result, the pdb files can be extracted with `extract_pdb()` sorted by the Rosetta score. 


```python
modeling.predict_3d_structure(f"{os.getcwd()}/secondary_prediction/dot_bracket.secstruct")
modeling.extract_pdb(5)
```

(content:modeling:parameter)=
## Modeling Parameter

```
rosetta_parameter = {
    "path_to_rosetta": "rna_denovo.default.macosclangrelease",
    "nstruct": 5,
    "minimize_rna": True,
    "cycles": 200
}
```

The parameters for the modeling are a reduced selection of Rosetta flag parameters, which are required for the tertiary structure modeling. 

- *path_to_rosetta* : This parameter refers to the Rosetta modeling module used, which must be specified depending on the operating system A containerized Docker environment should make this parameter obsolete in the future. 
- *nstruct*: Defines the number of structures to be modeled
- *minimize_rna*: Optimization of the RNA after assembling the fragments
- *cycles*: Number of Monte carlo sampling cycles 
