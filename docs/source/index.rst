MolSanitizer's documentation
#########################################

About MolSanitizer
******************

**MolSanitizer (msani)** is a Python package for preparation (remove salts, tautomer standardization, stereoisomers enumeration, protonation) and filtering undesirable substructures (PAINS, reactive functional groups) for drug discovery projects. MolSanitizer is also capable of building DB2 files (for DOCK3.8) and PDBQT files (for AutoDock Vina) programs. OpenEye OEB files for FRED/HYBRID and general conformer formats (SDF, MOL2) can also be generated. OEB output requires the optional OpenEye toolkits and a valid licence.

Check out the :doc:`usage` section for further information, including
how to :doc:`installation` the project.

.. note::

   This project is under active development.
   For any feedback, please refer to the :doc:`feedback` section.

Contents
************************

.. toctree::
   :maxdepth: 2
   :caption: General Documentation

   introduction
   validation
   installation
   quickstart
   usage
   outputs
   config
   citation
   feedback

.. toctree::
   :maxdepth: 1
   :caption: Advanced Topics
   
   changes
   python
   contributors

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   api