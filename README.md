# MolSanitizer - A package to prepare SMILES databases
[![python](https://img.shields.io/badge/python-v3.10--3.12-blue)]()
[![anaconda](https://img.shields.io/badge/Anaconda.org-2.1.0-green.svg?style=flat-square)](https://docs.anaconda.com/anaconda/install/index.html)
[![Documentation](https://img.shields.io/badge/docs-0.2.2-orange)](https://msani.readthedocs.io/)
![GitHub forks](https://img.shields.io/github/forks/:user/:repo)
[![license](https://img.shields.io/badge/license-GPLv2-yellow)](LICENSE)


MolSanitizer is a package for preparation (remove salts, stereoisomers enumeration, protonation, ...) and filtering undesirable substructures (PAINS, reactive functional groups, ...) for drug discovery projects.

## Installation and getting started

We will set up the environment using [Anaconda](https://docs.anaconda.com/anaconda/install/index.html). Clone the current repository:

    git clone https://github.com/phonglam3103/MolSanitizer.git
    OR
    git clone https://ghp_token@github.com/phonglam3103/MolSanitizer.git  #Put your personal token so that you don't have to sign in every time.
    
Example of how to set up a working conda environment to run the code on Mac OS and Linux:
    
    conda env create -f MolSanitizer/environment.yml # Trick: use mamba (if you have) for much faster installation
    conda activate msani
    pip install -e MolSanitizer

For Windows users, we recommend them to follow the instruction [here](https://msani.readthedocs.io/en/latest/installation.html). More information on the installation and dependencies could be found in the same page

## New C++ implementation since October 18, 2025

Since version 0.5.0, MolSanitizer has transferred the heavily demanding parts to C++, hence more dependencies are required to build and compile the program. If the users happen to have created the conda environment using the above instruction, they need to install the additional dependencies for C++:

    mamba env update --name msani --file environment.yml

After that, re-install the package:

    pip install -e MolSanitizer

## Documentation

To start to use msani, use the `-h` or `--help` flag for available options:

    msani -h

Documentation on the theory behind MolSanitizer and how to use it can be found [here](https://msani.readthedocs.io)

<img src="./plots/Workflow.png" width="1000">

## Notes about AMSOL

By default, all the dependencies are automatically installed by conda and pip, except for AMSOL. The user is asked to place the compiled version of (named `amsol7.1`) to [msani`/amsol](msani/amsol). In that folder, there will be a README on how to compile it on the modern Linux systems.

*As for the current evaluation version, the precompiled AMSOL version is provided. It will be removed once the repository is publicly available.*


## Contribution

We warmly welcome contributions of all kinds, whether it's reporting a bug, suggesting a feature, or developing new functionality. To get started, please refer to our [CONTRIBUTING.md](CONTRIBUTING.md) guide, which details the steps for contributing, from opening an issue to submitting a pull request.


## Feedback
MolSanitizer is a rule-based program that relies on our experience from previous drug discovery projects. We are committed to continuously improving the program's performance by adding more rules to the filters and tautomers/protonation. If you have any ideas or suggestions, please don't hesitate to open an issue or contact us.

## Contact
1. Thua-Phong Lam, phong.lam@icm.uu.se
2. Szymon Pach, szymon.pach@icm.uu.se
3. Israel Cabeza de Vaca Lopez, israel.cabezadevaca@icm.uu.se
