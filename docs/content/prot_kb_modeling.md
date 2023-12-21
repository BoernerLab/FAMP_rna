(content:protocol:modeling)=
# Protocol template modeling

In this protocol for knowledge-based modeling, manual structure assembly from an existing PDB structure is described 
using an example RNA sequence. 

The sequence11 of the construct consists of three structural elements. These are three hairpin loops located in the 25S 
subunit of the rRNA (PDB-ID: 3JCT) which together form a tertiary contact. H68 is the tetraloop containing a GAAA motif. 
H88 and H22 together form a kissing loop, which acts as a receptor for the tetraloop. The three hairpin loops are 
connected by poly adenine linkers. This results in the following primary structure:

**<font color="#f0b97a">H68</font>** – **poly-A-linker(8xA)** - **<font color="#d28282">H88</font>** – **poly-A-linker(11xA)** - **<font color="#a0c27c">H22</font>**

The sequence is saved in a FASTA file and is used as input to the modeling pipeline. 

```{admonition} Fasta File
:class: note 
<p style="font-family:monospace">>KLTL sequence construct <br>
<font color="#f0b97a">ugaagaaauuca</font>aaaaaaaa<font color="#d28282">gcucggaauuugagc</font>aaaaaaaaaaaa<font color="#a0c27c">cggugguaaauuccaucg</font></p> 
```

The secondary structure was predicted with RNAFold by using the `predict_2D_structure()` function of the pipeline. 
The results are stored in a dot-bracket-format. This shows canonical base pairs (round brackets) and free bases (dots). 
The kissing loop can’t be predicted due to computational limits of predicting pseudoknots. Here we introduce the 
kissing loop as pseudoknot knowledge based onto the secondary structure with square brackets.


```{admonition} Secondary Structure
:class: note 
<p style="font-family:monospace">
Result RNAFold: <br>
((((....))))........(((((.....)))))............((((((......)))))) <br>
Introducing the KL as pseudoknot: <br>
((((....))))........((((.[[[[[[))))............(((((..]]]]]]))))) <br>

</p>
```

Secondary structure visualization was performed with Forna: <a href="http://rna.tbi.univie.ac.at/forna/" target="_blank">http://rna.tbi.univie.ac.at/forna/ </a>

Tertiary structures were predicted with the FARFAR2 extension of the Rosetta modelling suite. The module for RNA modeling, based on FARFAR2, is called rna_denovo.
The following objects are necessary for the modeling process:
-	the sequence in a fasta file
-	Rosetta parameter as a flag file or parameter dictionary in FAMP
The flag file contains all additional flags of the rna_denovo command

```{admonition} Rosetta Flag File
:class: topic 
<p style="font-family:monospace; text-align:left">
-fasta BTL.fasta <br>
-minimize_rna true <br>
-cycles 20000 <br>
-nstruct 20 <br>
-save_times <br>
-secstruct "((((....))))........((((.[[[[[[))))............(((((..]]]]]])))))" <br>
-superimpose_over_all <br>
-out:file:silent BTL_modeling.out <br>
</p>
```

The command:
```bash
~/Documents/Rosetta/rosetta_bin_mac_2021.16.61629_bundle/main/source/bin/rna_denovo.default.macosclangrelease @flags
```

runs the modeling process. The path to the Rosetta source folder must be specified depending on the installation 
location and the program rna_denovo.default must be specified. The @flags calls the flag file and complements the 
command.

When Rosetta finished the modeling, the results are written to BTL_modeling.out. To export PDB structures, the out 
file must be converted using the extract_pdbs function. 

Rossetta calculates a SCORE for every modeled structure based on its own energy functions. So the structure with the 
best SCORE maybe contains the best structure. The models were sorted depending on its score and then exported as pdb 
structures. 
Since the structures with the best SCORE are desired, these values are transferred first into a file and sorted 
descending:

```bash
grep '^SCORE' BTL_modeling.out > test.sc
sort -k1,1 -k2n test.sc
```
Now the best structures can be chosen, the names of the structures are noted and collected in a .tag file. With the 
help of the .tag file the structures which are deposited in the .tag file are then converted into PDB structures with:

```bash
~/Documents/Rosetta/rosetta_bin_mac_2021.16.61629_bundle/main/source/bin/extract_pdbs.default.macosclangrelease 
-in:file:silent BTL_modeling.out -in:file:tagfile tags.tag
```
Afterwards, the structures can be selected with which further work is to be done. In case of the KL-TLGAAA, structures 
were predicted where the poly-A-linker blocking the kissing loop binding site. Others did not but had lower SCORE 
values. Since no model predicted a binding of the KL to TL a remodeling for the bound state is necessary. 

With the FAMP pipeline, 3D structure prediction is simplified by two functions. With predict_3D_structure() the 
modeling process is started. The flags are automatically generated using the defined parameters. 
As further input the predicted secondary structure is needed. This should be specified as a .sectstruct file. 
With extract_pdb() the output can be converted into a PDB file after the modeling process. 
The results of the structure modeling are stored in rosetta_results/out/1/ folder and can be further processed from 
there.
