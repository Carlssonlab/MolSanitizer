.. _config:
Config file
============

The config file is designed to be reused across different projects. It is a YAML file that contains the configuration for the MolSanitizer. The file can be customized to suit specific needs, such as defining the chemical space of interest, setting up the protonation, tautomerization, 3D generation options, and specifying the output format.


Nearly all the options are available, except ones that are designed to be triggered once only like ``--version``, ``--help``, or ``--smiles``.


The keywords in the config file need to match that with the long_name of the command line options. For example, if you want to set the ``--protonation`` option, you would use `protonation: true` in the config file. 

Template
~~~~~~~~~~

To create a template for the config file, you can run the following command:

.. code-block:: console

    $ msani --create_config

The file will be created in the current working directory with the name ``config.yaml``. You can then edit this file to customize the settings for your project. Lines start with '#' are commented out. Here is the template of the config file:

.. code-block:: yaml

    # Input-output options
    input_files: [file1.smi, file2.smi]
    extended: false

    # Filtering options (other options include hba, hbd, mw, chiral)
    removesalts: false
    unwanted: [regular, optional]
    pains: false
    # ha: 17-25
    # logp: 350

    # SMILES processing options
    tautomers: true
    protonation: true
    pH: 7
    pH_range: 0
    # prot_lib: None #In case of modified protonation library
    # tau_lib: None  # In case of modified tautomer library

    # Conformer generation options
    # gen3d: false
    # format: [db2.tgz]
    # method: rdkit
    # numconfs: 2000
    # energywindow: 25
    # rigid: None
    # torsion: None # In case of modified torsion definitions

Multiple-value arguments
~~~~~~~~~~~~~~~

The following arguments accept multiple values and can be specified as a list in the config file:

- input_files:  list of input files
- unwanted: choice of category for unwanted substructures filtering (options: all, regular, special, optional)
- format: 3D output formats (options: db2.tgz, db2, mol2, pdbqt, sdf)

These arguments can be specified as a list in the config file in the Python format:

.. code-block:: yaml

    input_files: [file1.smi, file2.smi, file3.smi]

or as multiple lines:

.. code-block:: yaml

    input_files:
    - file1.smi
    - file2.smi

Usage
~~~~~~~~~~~

To use the config file, you can use the ``-c`` or ``--config`` option followed by the path to your config file when running the MolSanitizer command. For example, if you have a config file named `config.yaml`, you can run:

.. code-block:: console

    $ msani -c config.yaml

This will read the configuration from the ``config.yaml`` file and apply the settings specified in it. You can also specify additional command line options if needed, which will override the settings in the config file.
