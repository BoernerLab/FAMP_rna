# Installation

The code can be downloaded under the green button "Code" or cloned via the terminal:
```
git clone https://github.com/felixErichson/FAMP_rna.git
```
<br>

## Requirements
As a dependency of this notebook the following programms should be preinstalled. Please follow the installation instructions of the tools.

[**ViennaRNAPackge**](https://www.tbi.univie.ac.at/RNA/documentation.html#install)

[**Rosetta**](https://new.rosettacommons.org/docs/latest/build_documentation/Build-Documentation)

[**GROMACS**](https://manual.gromacs.org/documentation/2021.2/install-guide/index.html)

To start MD simulations using the pipeline, the force field must be added to the GROMACS installation directory.
Add the directory `amber14sb_OL15.ff` from `FAM_Piline/src/scripts/gromacs` (this repository) to `<GROMACS_installtion>/share/gromacs/top` 

To avoid GROMACS errors, when working with dyes, you should also 

Make sure that the programms are added to the bashrc or zshrc file.

Please create an anaconda environment for this Jupyter Notebook by importing the environment.yml file. 

Before you start open Pymol in your environment and install the plugins FRETraj and FRETlabel. Locate the 
package with `fretlabel --path` and remember this path. Then open PymMol and go to Plugin -> Plugin manager -> Install 
New Plugin -> Choose file ... -> and select fretlabel_gui.py, wich can be found at the path from the step bevore. The 
same works for the FRETraj GUI.