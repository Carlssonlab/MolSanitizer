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

# Usage

## Input

The program requires a white-space or tab-delimited file containing two columns (SMILES, moleculeID) without headers.

```
COCCC(=O)Nc1ncc(s1)Br  CP000000418470
C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597
c1c(coc1Br)C(=O)NC2CCSC2  CP000001645677
c1c(c([nH]n1)C(=O)NCC2(CC2)N)Br  CP000001647414
```

## **Overview**
The pipeline contains five preparation and/or filtering steps, which could be used simultaneously to prepare the database:

This is an example for a lazy pipeline which use all the preparation steps:

```bash
msani -i example.smi --removesalts --tautomers --pains --unwanted all --stereoisomers --protonation
```

Use the `--help (-h) flag` for more information.

```bash
msani -h
```

By default, the program produces a new file with a **_clean** suffix. If the PAINS or unwanted filters are applied, the rejected molecules with reason of rejection will be output by **_rejected** suffix.

The program by default will conduct the preparation and filtering in the order as below:

```
    Put an image of the workflow here.
```

## **Step 1: Remove salts**
 To use the remove salts function, simply use `--removesalts flag`. The program uses a predefined salt list in **Data/salt_stripping.txt** to remove the salts, which contains both organic and inorganic salts that are commonly used in medicinal chemistry. 

*Caution:* if the entry is an organic salt (eg. sodium acetate CH_{3}COO^{-}Na{+}), all the entry will be removed.
```bash
msani -i example.smi --removesalts
```

## **Step 2: Tautomer enumeration**
The tautomers could be generated using a `--tautomers flag`. The program uses a predefined SMARTS rules to generate the possible tautomers of more specifically the conjugated ring systems containing Nitrogen.
```bash
msani -i example.smi --tautomers
```

## **Step 3: PAINS filtering**
Molecules that contain PAINS substructures could be efficiently eliminated using the `--pains flag`. The violated structures would be stored in the **_rejected** file.

```bash
msani -i example.smi --pains
```

## **Step 4: Unwanted substructures filtering**

## **Step 5: Stereoisomers enumeration**

## **Step 6: Protonation**
