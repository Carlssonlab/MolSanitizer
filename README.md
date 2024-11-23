# MolSanitizer - A package to prepare SMILES databases
[![python](https://img.shields.io/badge/python-v3.11-blue)]()
[![anaconda](https://img.shields.io/badge/Anaconda.org-2.1.0-green.svg?style=flat-square)]()
[![Documentation](https://img.shields.io/badge/docs-0.2.1-orange)](https://msani.readthedocs.io/)
![GitHub forks](https://img.shields.io/github/forks/:user/:repo)


MolSanitizer is a package for preparation (remove salts, stereoisomers enumeration, protonation, ...) and filtering undesirable substructures (PAINS, reactive functional groups, ...) for drug discovery projects.

Documentation on how to install, usage, and theory behind MolSanitizer can be found [here](https://msani.readthedocs.io).

By default, all the dependencies are automatically installed by conda and pip, except for AMSOL. The user is asked to place the compiled version of (named `amsol7.1`) to [MolSanitizer/amsol](MolSanitizer/amsol). In that folder, there will be a README on how to compile it on the modern Linux systems.

*As for the current evaluation version, the precompiled AMSOL version is provided. It will be removed once the repository is publicly available.*

# Feedback
MolSanitizer is a rule-based program that relies on our experience from previous drug discovery projects. We are committed to continuously improving the program's performance by adding more rules to the filters and tautomers/protonation. If you have any ideas or suggestions, please don't hesitate to open an issue or contact us.

# Contact
1. Thua-Phong Lam, phong.lam@icm.uu.se
2. Szymon Pach, szymon.pach@icm.uu.se
3. Israel Cabeza de Vaca Lopez, israel.cabezadevaca@icm.uu.se
