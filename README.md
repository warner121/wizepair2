# Matched Molecular Pairs Analysis (MMPA) based upon Maximum Common Substructure (MCS)
Python (re-)implementation of the algorithm described in: **WizePairZ: A Novel Algorithm to Identify, Encode, and Exploit Matched Molecular Pairs with Unspecified Cores in Medicinal Chemistry** by Daniel J. Warner, Edward J. Griffen, and Stephen A. St-Gallay: *J. Chem. Inf. Model. 2010, 50, 8, 1350–1357, August 6, 2010.*

### Installation
This implementation builds on the open source cheminformatics tooklik RDKit and the graph toolkit networkx. You must have a working copy of RDKit with python bindings installed. If you have issues following the steps below refer to the RDKit installation guide here: https://www.rdkit.org/docs/Install.html.  

1. Clone the wisepair2 repository e.g.  
```
$ git clone https://github.com/warner121/wizepair2.git
```
2. Create your anaconda environment from the provided yaml.  
```
$ conda env create -f environment.yaml
```
3. Activate you new anaconda session.  
```
$ conda activate wizepair-env
```

### Execution
Access the help by executing the program with the --help argument:
```
(my-rdkit-env) ~/wizepair2$ python mmpa -h

Usage: mmpa [options]  
  
Copyright 2013 Daniel Warner  
Licensed under the Apache License 2.0  
http://www.apache.org/licenses/LICENSE-2.0  

Options:  
  --version            show program's version number and exit  
  -h, --help           show this help message and exit  
  -i FILE, --in=FILE   set input path [default: ./mmpa.in]  
  -o FILE, --out=FILE  set output path [default: ./mmpa.out]  
  -v, --verbose        set verbosity level [default: none]
```
Any .csv file featuring the columns "Molecule_L" and "Molecule_R" will serve as a suitable input file. An example input file is available here: https://github.com/warner121/wizepair2/blob/master/mmpa.in. 
Simply specify input and output files and execute.
```
(my-rdkit-env) ~/wizepair2$ python mmpa -i mmpa.in -o mmpa.out
infile = mmpa.in
outfile = mmpa.out
```
The results may then be found in the designated output file.
