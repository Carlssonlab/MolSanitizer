.. _config:

Config file
============

Personal defaults
~~~~~~~~~~~~~~~~~

MolSanitizer loads personal defaults automatically for ``msani`` and
``msani_batch``. You do not need to locate or edit files inside the installed
package. If no personal configuration exists, bundled defaults are used silently.
This is normal behavior and does not produce a warning in CLI runs, library
calls, or unit tests. Use ``msani --help`` for configuration guidance or
``msani config show`` to inspect defaults and their sources. Ordinary runs do
not create configuration files. Missing or invalid explicitly selected
configurations still produce errors.

Create your personal configuration with:

.. code-block:: console

    $ msani config init

The command prints the location of ``msani_configurations.yaml``. Edit that
file to change defaults for subsequent runs. Standard locations are:

.. list-table:: User configuration directory
   :header-rows: 1

   * - Operating system
     - Directory
   * - Linux
     - ``~/.config/msani/`` (or ``$XDG_CONFIG_HOME/msani/`` when set)
   * - macOS
     - ``~/Library/Application Support/msani/``
   * - Windows
     - ``%LOCALAPPDATA%\msani\``

For example, your file may contain only the settings you want to override:

.. code-block:: yaml

    CORINA: /opt/corina/corina
    SLURM_ACCOUNT: my-project
    NUMCONFS: 1000
    MAX_JOBS: 100

Missing keys retain bundled defaults. Setting names in this personal file
use the uppercase names from the generated template. Unknown keys, invalid
values, and unreadable files produce errors identifying the file or setting.
The bundled CORINA default is ``corina`` on PATH, and the default SLURM
account is unset.

Custom locations
^^^^^^^^^^^^^^^^

To create a configuration elsewhere and remember it, or select an existing
configuration, run:

.. code-block:: console

    $ msani config init --path /shared/mygroup/msani.yaml
    $ msani config use /shared/mygroup/existing.yaml

``use`` validates the existing file before remembering it. Paths containing
spaces should be quoted. Paths are stored as absolute paths after expanding
``~``. A small ``settings.json`` file in the standard user configuration
directory remembers the selected location. It is outside the installation
and survives wheel upgrades. Environments using the same user configuration
directory share this preference.

Additional commands:

.. code-block:: console

    $ msani config path
    $ msani config show
    $ msani config reset-path

``path`` prints the selected YAML location, even if the default file has not
yet been created. ``show`` displays the merged personal and bundled defaults
and the source of each setting; it does not include a particular job's CLI
arguments. ``reset-path`` removes the remembered location without deleting
any YAML files. ``init`` without ``--path`` creates and selects the standard
file. Existing files are preserved unless ``init --force`` is explicitly used.
These configuration commands are also available through ``msani_batch``.

For a temporary override, set the ``MSANI_CONFIG`` environment variable to a
YAML filename. This takes priority over the remembered location without
changing it. For example, on Linux or macOS:

.. code-block:: console

    $ MSANI_CONFIG=/shared/project/msani.yaml msani -i input.smi

In Windows PowerShell:

.. code-block:: powershell

    $env:MSANI_CONFIG = 'C:\project\msani.yaml'
    msani -i input.smi

Unset the variable to return to the remembered or standard location. If a
remembered or environment-selected file is missing or invalid, MolSanitizer
reports an error instead of silently using different defaults.

Precedence and batch jobs
^^^^^^^^^^^^^^^^^^^^^^^^^

Settings are applied in this order, from highest to lowest priority:

1. Explicit command-line arguments.
2. Job settings passed with ``--config`` (described below).
3. Personal defaults from ``MSANI_CONFIG``, the remembered location, or the
   standard location, in that order. Only one personal file is selected.
4. Bundled defaults.

Batch submission writes ``msani_defaults.yaml`` beside ``submit_msani.sh``.
The generated script selects this snapshot through ``MSANI_CONFIG`` so
changes to personal defaults do not affect already submitted jobs. Keep the
snapshot with the job files. Explicit job options still take precedence.

Job configuration
~~~~~~~~~~~~~~~~~

The job configuration described below is separate from personal defaults:
``--create_config`` creates a job template using CLI option names, whereas
``config init`` creates persistent defaults using uppercase setting names.

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
~~~~~~~~~~~~~~~~~~~~~~~~

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
