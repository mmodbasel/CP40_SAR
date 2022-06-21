# CP40 SAR

This repository contains the source code of the program interaction-analyzer.x to evaluate the per-residue intramolecular and intermolecular hydrogen bonding and hydrophobic interactions of C3b with Cp40 and analogues thereof.
The calculation of the hydrophobic interaction component is based on the hydrophobicity score by Li et al. (DOI: [10.1002/prot.23106](https://doi.org/10.1002/prot.23106)). In our implementation, the backbone carbonyl carbon atoms are not  considered hydrophobic. The directional H-bonding interaction components are calculated using parameters by Vedani et al. (DOI: [10.1021/ja00168a021](https://doi.org/10.1021/ja00168a021)). 
The program source code consist of several files that need to reside in one directory:
- main program source code:  interaction-analyzer.cpp
- three header files: definitions.h, points_vectors.h, structures.h

## Requirements
A recent version of the g++ compiler is needed to compile the code.

## Compilation
In order to compile the source code into executable program, following command shall be used (on a Linux operating system):
```
g++ -Wall -pedantic -O2 interaction-analyzer.cpp -o interaction-analyzer.x
```

## Usage

The compiled program *interaction-analyzer.x* can be executed from the command line providing an MD-trajectory snapshot in the PDB formatted file (*C3b---CP40_MD_frame0.pdb* file is provided as an example) as a command line parameter:

```
./interaction-analyzer.x  C3b---CP40_MD_frame0.pdb
```

**Important note:** The ligand part (Cp40 or alike cyclic peptide) in the complex is automatically recognized by the program by finding of a tryptophan (TRP) residue with number 7 or 8. This is hard-coded on the line 797 of the *interaction-analyzer.cpp* source file and can be simply edited if needed.

## Output

The program prints all individual interaction energies between ligand and protein residues in human readable format within a fraction of a second. Energy unit used is kcal/mol.

## Reference

When using this source code, please cite the following paper: *to be filled*

