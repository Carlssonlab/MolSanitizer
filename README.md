# MolSanitizer - A package to prepare SMILES databases
[![python](https://img.shields.io/badge/python-v3.11-blue)]()
[![anaconda](https://img.shields.io/badge/Anaconda.org-2.1.0-green.svg?style=flat-square)]()
[![Documentation](https://img.shields.io/badge/docs-0.2.2-orange)](https://msani.readthedocs.io/)
![GitHub forks](https://img.shields.io/github/forks/:user/:repo)
[![license](https://img.shields.io/badge/license-GPLv2-yellow)]()


MolSanitizer is a package for preparation (remove salts, stereoisomers enumeration, protonation, ...) and filtering undesirable substructures (PAINS, reactive functional groups, ...) for drug discovery projects.

## Installation and getting started

We will set up the environment using [Anaconda](https://docs.anaconda.com/anaconda/install/index.html). Clone the
current repository:

    git clone https://github.com/phonglam3103/MolSanitizer.git
    
Example of how to set up a working conda environment to run the code:

    conda env create -f MolSanitizer/environment.yml
    conda activate msani
    pip install -e MolSanitizer

More information on the installation and dependencies could be found [here](https://msani.readthedocs.io/en/latest/installation.html).


## Documentation

To start to use MolSanitizer, use the `-h` or `--help` flag for available options:

    msani -h

Documentation on the theory behind MolSanitizer and how to use it can be found [here](https://msani.readthedocs.io/)

## Notes about AMSOL

By default, all the dependencies are automatically installed by conda and pip, except for AMSOL. The user is asked to place the compiled version of (named `amsol7.1`) to [MolSanitizer/amsol](MolSanitizer/amsol). In that folder, there will be a README on how to compile it on the modern Linux systems.

*As for the current evaluation version, the precompiled AMSOL version is provided. It will be removed once the repository is publicly available.*


## Contribution

We warmly welcome contributions of all kinds, whether it's reporting a bug, suggesting a feature, or developing new functionality. To get started, please refer to our [CONTRIBUTING.md](CONTRIBUTING.md) guide, which details the steps for contributing, from opening an issue to submitting a pull request.


## Feedback
MolSanitizer is a rule-based program that relies on our experience from previous drug discovery projects. We are committed to continuously improving the program's performance by adding more rules to the filters and tautomers/protonation. If you have any ideas or suggestions, please don't hesitate to open an issue or contact us.

## Contact
1. Thua-Phong Lam, phong.lam@icm.uu.se
2. Szymon Pach, szymon.pach@icm.uu.se
3. Israel Cabeza de Vaca Lopez, israel.cabezadevaca@icm.uu.se
