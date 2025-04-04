# EirVS - A package to prepare SMILES databases
[![python](https://img.shields.io/badge/python-v3.9--3.12-blue)]()
[![anaconda](https://img.shields.io/badge/Anaconda.org-2.1.0-green.svg?style=flat-square)](https://docs.anaconda.com/anaconda/install/index.html)
[![Documentation](https://img.shields.io/badge/docs-0.2.2-orange)](https://EirVS.readthedocs.io/)
![GitHub forks](https://img.shields.io/github/forks/:user/:repo)
[![license](https://img.shields.io/badge/license-GPLv2-yellow)](LICENSE)


EirVS is a package for preparation (remove salts, stereoisomers enumeration, protonation, ...) and filtering undesirable substructures (PAINS, reactive functional groups, ...) for drug discovery projects.

## Installation and getting started

We will set up the environment using [Anaconda](https://docs.anaconda.com/anaconda/install/index.html). Clone the current repository:

    git clone https://github.com/phonglam3103/EirVS.git
    OR
    git clone https://ghp_token@github.com/phonglam3103/EirVS.git  #Put your personal token so that you don't have to sign in every time.
    
Example of how to set up a working conda environment to run the code:
    
    conda env create -f EirVS/environment.yml # Use mamba instead of conda for much faster installation
    conda activate eirvs
    pip install -e EirVS

More information on the installation and dependencies could be found [here](https://eirvs.readthedocs.io/en/latest/installation.html).


## Documentation

To start to use EirVS, use the `-h` or `--help` flag for available options:

    eirvs -h

Documentation on the theory behind EirVS and how to use it can be found [here](https://eirvs.readthedocs.io)

<img src="./plots/Workflow.png" width="1000">

## Notes about AMSOL

By default, all the dependencies are automatically installed by conda and pip, except for AMSOL. The user is asked to place the compiled version of (named `amsol7.1`) to [EirVS`/amsol](EirVS/amsol). In that folder, there will be a README on how to compile it on the modern Linux systems.

*As for the current evaluation version, the precompiled AMSOL version is provided. It will be removed once the repository is publicly available.*


## Contribution

We warmly welcome contributions of all kinds, whether it's reporting a bug, suggesting a feature, or developing new functionality. To get started, please refer to our [CONTRIBUTING.md](CONTRIBUTING.md) guide, which details the steps for contributing, from opening an issue to submitting a pull request.


## Feedback
EirVS is a rule-based program that relies on our experience from previous drug discovery projects. We are committed to continuously improving the program's performance by adding more rules to the filters and tautomers/protonation. If you have any ideas or suggestions, please don't hesitate to open an issue or contact us.

## Contact
1. Thua-Phong Lam, phong.lam@icm.uu.se
2. Szymon Pach, szymon.pach@icm.uu.se
3. Israel Cabeza de Vaca Lopez, israel.cabezadevaca@icm.uu.se
