# MolSanitizer - A package to prepare SMILES databases

MolSannitizer is a package for preparation (remove salts, stereoisomers enumeration, protonation) and filtering of undesireable substructures (PAINS, reactive functional groups) for drug discovery projects.

# Setup Environment

We will set up the environment using [Anaconda](https://docs.anaconda.com/anaconda/install/index.html). Clone the
current repository:

    git clone https://github.com/Isra3l/MolSanitizer.git
    

This is an example for how to set up a working conda environment to run the code:

    conda env create -f MolSanitizer/environment.yml
    conda activate msani
    pip install -e MolSanitizer

# Input

The program requires a white-space or tab delimited file containing two columns (SMILES, moleculeID) without headers.

```
COCCC(=O)Nc1ncc(s1)Br  CP000000418470  -22.99
C1CC(C(=O)NC1)SCCC=CBr  CP000000432409  -19.54
CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597  -21.29
c1c(coc1Br)C(=O)NC2CCSC2  CP000001645677  -19.28
c1c(c([nH]n1)C(=O)NCC2(CC2)N)Br  CP000001647414  -12.96
```

# Usage

## **Overview**
The pipeline contains five preparation and/or filtering steps, which could be used simultaneously to prepare the database:

This is an example for a lazy pipeline which use all the preparation steps:

```bash
msani -i example.smi --removesalts --pains --unwanted all --stereoisomers --protonation
```

