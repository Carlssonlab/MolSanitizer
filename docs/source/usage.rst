Usage
=================================

Preparing the input file
-----------------------------------

The program requires a white-space or tab-delimited file containing two columns (SMILES, moleculeID) without headers. 

Input files may be uncompressed, gzip-compressed, bzip2-compressed, or
xz-compressed. MolSanitizer detects gzip, bzip2, and xz from the file contents
(so the extension is optional) and streams decompressed rows to pandas in
chunks rather than loading the whole file into memory.

.. code-block:: console
   
   COCCC(=O)Nc1ncc(s1)Br  CP000000418470
   C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
   CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597
   c1c(coc1Br)C(=O)NC2CCSC2  CP000001645677
   c1c(c([nH]n1)C(=O)NCC2(CC2)N)Br  CP000001647414

Extended SMILES is supported (by the ``-e`` or ``--extended`` flag, but the database needs to be tab-separated.

.. code-block:: console

    CC[C@H]1[C@H](C(=O)N[C@H](C)CCCC(=O)NOCC(F)(F)F)CCN1C |&1:2,3|	Cmp0001
    CCC(CC(=O)N(CC)CCC(=O)N1CCO[C@H]2COC[C@H]21)C(F)F |&1:17,21|	Cmp0002
    CC(C)CC(CNC(=O)C1CSC1)C(=O)N[C@H]1C[C@@H](O)[C@H](F)C1 |&1:16,18,20|	Cmp0003

Overview
-----------------------------------

The pipeline contains six preparation and/or filtering steps, which can be used simultaneously to prepare the database:

This is an example of a lazy pipeline that uses all the preparation and processing steps:

.. code-block:: console

    $ msani -i example.smi --extended --lazy  # For extended SMILES format
    $ msani -i example.smi --lazy

    # This is equivalent to:
    # msani -i example.smi --removesalts --tautomers --pains --unwanted all --stereoisomers --protonation


By default, the program produces a new file with a **_clean** suffix. If the PAINS or unwanted filters are applied, the rejected molecules with the reason for rejection will be output with a **_rejected** suffix.

The program by default will conduct the preparation and filtering in the order below:

.. image:: _static/Workflow.png
   :width: 800px



Default settings
-----------------------------------

To customize defaults for future runs, use ``msani config init`` and edit the
personal ``msani_configurations.yaml`` at the location printed by the command.
Both single and batch mode load this file automatically. Use ``msani config use
PATH`` to remember a different configuration location, and ``msani config show``
to inspect effective defaults. Explicit command-line options override personal
defaults. See :doc:`config` for platform locations, precedence, and examples.


The users are asked to provide CORINA path if the he/she wants to use it for the generation of 3D coordinates. The path should be provided in the `CORINA` field.

Below is the default configuration file:

.. code-block:: yaml
    
    #===============SINGLE MODE================
    EMBED_METHOD: 'rdkit' # choose from 'rdkit', 'obabel', 'corina'
    CORINA: '/soft/corina/corina-4.2/corina'
    ENERGY_WINDOW: 25
    NUMCONFS: 2000
    MAX_STEREOISOMERS: 8
    TIMEOUT: 2
    PH: 7
    PH_RANGE: 0 # 0 means choose specific pH=7 (default), 2 means will sample pH 5 and 9

    #================BATCH MODE=================
    SLURM_ACCOUNT: 'PROJECT_NAME' # The account that will be charged by the SLURM cluster for running tasks
    LINES_PER_JOB: 200
    TIME_LIMIT: 96
    MAX_ARRAY_SIZE: 2000
    MAX_JOBS: 1000
    MAX_LIMIT_PROJECT: 5000


Reproducible setup by the configuration file
-----------------------------------------------

The program can be run with the configuration file by using the ``--config`` flag. The path to the configuration file should be provided. Refer to the :doc:`config` section for more details on the configuration file.

Help message
-----------------------------------

**Use the** ``--help (-h)`` **flag for more information.**

.. code-block:: console

    usage: msani [--input_files INPUT_FILES [INPUT_FILES ...]] [--input_list INPUT_LIST] [--smiles SMILES [SMILES ...]] [--extended] [--prefix PREFIX] [--synthon] [--removesalts] [--create_custom] [--custom CUSTOM]
                [--unwanted [{all,regular,special,optional} ...]] [--pains] [--metal {strict,soft,off}] [--metal_error_file METAL_ERROR_FILE] [--neutralize | --no-neutralize | -neu] [--stereoisomers | --no-stereoisomers | -st]
                [--max_isomers MAX_ISOMERS] [--stereo_timeout STEREO_TIMEOUT] [--tautomers] [--extended-tautomers] [--protonation] [--pH PH] [--pH_range PH_RANGE] [--standardize] [--gen3d]
                [--format [{db2,db2.tgz,pdbqt,sdf,mol2,oeb,oeb.lib} ...]] [--method {rdkit,obabel,corina}] [--corinaPath CORINAPATH] [--numconfs NUMCONFS] [--timeout TIMEOUT] [--energywindow ENERGYWINDOW] [--nringconfs NRINGCONFS]
                [--mode {fixed,random,ignoretorlib}] [--allowNonring] [--rmsd RMSD] [--config CONFIG] [--create_config] [--lazy] [--numcores NUMCORES] [--help] [--help_advanced] [--version]

    MolSanitizer - A package to prepare SMILES databases

    Input and output options:
    --input_files, -i     Input files containing chemical structures (plain, gzip, bzip2, or xz)
    --input_list, -il     Path to a text file containing one or more input file paths (one per line).
    --smiles, -s          Input SMILES strings
    --extended, -e        Extended SMILES reading (tab-separated files supported only).
    --prefix, -pre        Prefix for the output files. (default: input file name).
    --synthon, -stn       Synthon mode (Additional metadata about the capping groups required)

    Filtering options:
    Supported formats for descriptor-based filters (ha, logp, hba, hbd, mw, tpsa, fsp3, chiral):
        Range: Specify a range using two values (e.g., "17-25").
        Greater / Less than or equal to: Use >= or <= (e.g., ">=17", "<=25").
        Greater than / Less than: Use > or < (e.g., ">17", "<25").
        Exact match: Match a specific value (e.g., 17).
        For logP, the exact match format applies as 'less than or equal to'.

        Use --ha, --logp, --hba, --hbd, --mw, --tpsa, --fsp3, --chiral to apply these filters.

    --removesalts         Remove salts from the structures.
                            Small fragments within the same molecule are also removed.
    --create_custom       Generate a template for customized substructure filtering.
    --custom              Filter out unwanted substructures using a customized list.
                            To generate an example list, use --create_custom.
    --unwanted            Filter out unwanted substructures using the default list
                            (Options: all, regular, special, optional).
    --pains               Remove PAINS violations from the structures.

    SMILES processing options:
    --metal, -mtl         Connect dissociated metal complexes (metal and ligands given as separate
                            fragments) with dative bonds, before any other step (default: strict).
                            - strict: one structure per record, filling the coordination sphere by donor priority.
                            - soft: every plausible binding mode, as separate records (id_1, id_2, ...).
                            - off: leave metal records unchanged.
    --metal_error_file    File for metal records that cannot be prepared (SMILES ID reason).
                            Default: <prefix or input name>_metal_error.smi; batch jobs merge theirs
                            into in/removed/metal_error.smi.
    --neutralize, -neu    Neutralize molecules.
                            Will be applied after removesalts and before tautomerization/protonation (use --no-neutralize to disable).
    --stereoisomers, -st  Stereoisomers enumeration.
                            Will be applied by default when gen3d is on (use --no-stereoisomers to disable)
    --max_isomers, -ms    Maximum number of stereoisomers to consider (default: 8)
    --stereo_timeout, -sto
                            Per-molecule timeout in seconds for stereoisomer enumeration (default: 60).
                            If the timeout is reached, any stereoisomers found up to that point are kept; if none are found, the input SMILES is kept unchanged.
    --tautomers, -tau     Tautomers enumeration.
    --extended-tautomers, -et
                            Extended tautomers enumeration using additional tautomerization rules (will also activate '--tautomers')
                            This option will enumerate to less probable tautomeric forms.
    --protonation, -prot  (De)protonate the structures
    --pH, -p              pH for the protonation (default: 7)
    --pH_range, -r        pH range for the protonation (default: 0)
    --standardize, -std   Standardize structures for machine learning using RDKit

    Generate 3D conformers options:
    --gen3d, -3d          Generate 3D conformers
    --format, -f          Output file format. Multiple formats simultaneously supported.
                            (Default: db2.tgz - Options: sdf, db2, db2.tgz, mol2, pdbqt, oeb, oeb.lib.)
                            oeb: one .oeb.gz per molecule; oeb.lib: one library per input file.
                            Both need the optional OpenEye toolkits and OE_LICENSE.
    --method, -m          Embedding method (default: rdkit - options: rdkit, obabel, corina)
    --corinaPath          Path to the CORINA executable 
    --numconfs, -nconfs   Maximum number of conformers to generate (default: 2000)
    --timeout, -to        Timeout for the initial embedding for each entry before using OpenBabel
                            Default: 2 minutes
    --energywindow, -w    Energy window for sampling the conformations (default: 25 kcal/mol)
    --nringconfs, -nr     Maximum number of ring conformers to generate (default: 1)
    --mode, -mode         Mode for generating conformers
                            Default: fixed - Options: fixed, random, ignoretorlib
    --allowNonring        Allow the full sampling of non-ring compounds (default undersample to 30 confs).
    --rmsd, -rmsd         Minimum RMSD between two conformers (default: 0.5 Å).

    Miscellaneous:
    --config, -c          Path to the YAML configuration file
    --create_config       Create a template for the configuration file
    --lazy                Implement all the processing and preparation steps
    --numcores, -j        Number of cores to use for parallel processing (default: 4)
    --help, -h            Show this help message and exit
    --help_advanced, -xh  Show advanced help message with additional options
    --version, -v         Show the current version of MolSanitizer

    Example input file (space or tab-separated file):
            COCCC(=O)Nc1ncc(s1)Br  CP000000418470
            C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
            CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597

    Extended SMILES is supported with the -e flag (SMILES and IDs need to be tab-separated):
            CC[C@H]1[C@H](C(=O)N[C@H](C)CCCC(=O)NOCC(F)(F)F)CCN1C |&1:2,3|  Cmp0001
            CCC(CC(=O)N(CC)CCC(=O)N1CCO[C@H]2COC[C@H]21)C(F)F |&1:17,21|    Cmp0002
            CC(C)CC(CNC(=O)C1CSC1)C(=O)N[C@H]1C[C@@H](O)[C@H](F)C1 |&1:16,18,20|    Cmp0003

    Example usage:
        msani -i example.smi --removesalts --pains --unwanted all --stereoisomers --protonation
        msani -i example.smi --logp "<=500" --hba "<=10" --hbd "<=5" --mw "<=500" -3d -f pdbqt
        msani -i example.smi --pains --unwanted regular optional --stereoisomers --protonation
        msani -i example.smi --pains --unwanted all -prot -p 7 -tau -st -3d -f db2.tgz


Available filters and preparation steps
-------------------------------------------

1. Metal complexes
~~~~~~~~~~~~~~~~~~

Metal complexes are often stored with the metal and its ligands as separate fragments, for example ``C1CC(C1)(C(=O)O)C(=O)O.[NH2-].[NH2-].[Pt+2]`` for carboplatin. Before any other step, MolSanitizer connects them into one structure with dative bonds, so that salt removal and the filters see a complex rather than loose fragments. Only records containing a metal (Sc-As, Y-Sb, La-Bi) are touched, and complexes already drawn as one fragment are left as they are.

- ``--metal strict`` (default) writes one structure per record, filling the coordination sphere by donor priority. When several binding options tie, the lexicographically first canonical SMILES is kept and a warning is written to the log.
- ``--metal soft`` writes every plausible binding mode, for example both the C- and the N-bound cyanide, as separate records numbered ``id_1``, ``id_2``, ... (at most ``--metal_max_variants``, default 16).
- ``--metal off`` leaves metal records unchanged.

Records that cannot be prepared, such as polynuclear clusters or ligands without a donor atom, are removed and written with their reason to ``<prefix or input name>_metal_error.smi``; use ``--metal_error_file`` to choose another file. Batch jobs merge theirs into ``in/removed/metal_error.smi``.

*Caution:* the regular unwanted-substructure list (``--unwanted``) removes compounds of the common organometallic elements (``most_popular_organometallic_compounds``), connected or not.

.. code-block:: console

    $ msani -i example.smi --metal soft

2. Remove salts
~~~~~~~~~~~~~~~


To use the remove salts function, simply use the ``--removesalts`` flag. The program uses a predefined salt list in `msani/Data/salt_stripping.txt <https://github.com/phonglam3103/msani/blob/main/msani/Data/salt_stripping.txt>`_ to remove the salts, which contain both organic and inorganic salts commonly used in medicinal chemistry. Since the 0.2.3 version, Salt Remover will also remove the smaller fragments in the same molecule entry and only retain the largest one.

*Caution:* If the entry is an organic salt (e.g., sodium acetate CH\ :sub:`3` COO\ :sup:`-` Na\ :sup:`+`), the whole entry will be removed.

.. code-block:: console

    $ msani -i example.smi --removesalts

3. Neutralization
~~~~~~~~~~~~~~~~~~~~~~

Molecules can be neutralized using the ``--neutralize`` or ``-neu`` flag. The neutralization will be applied after the salt removal and before the tautomerization/protonation steps. If the user does not want to neutralize the molecules, he/she can use the ``--no-neutralize`` flag.

.. code-block:: console

    $ msani -i example.smi --neutralize
    $ msani -i example.smi -neu  # Short version
    $ msani -i example.smi --removesalts --no-neutralize  # To disable neutralization when using removesalts

4. Tautomers enumeration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The tautomers could be generated using the ``--tautomers`` flag. msani uses a two-step approach for the enumeration of tautomers. First, the canonical tautomer from the scoring function of ``rdMolStandardize.TautomerEnumerator`` is used. Then, the exceptions are corrected using the expert-curated SMARTS rules. The SMARTS rules are readily accessible at `msani/Data/tautomers.txt <https://github.com/phonglam3103/msani/blob/main/msani/Data/tautomers.txt>`_.

.. code-block:: console

    $ msani -i example.smi --tautomers

5. Descriptor-based filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following descriptors are supported for filtering: heavy atoms (HA), logP, hydrogen bond acceptors (HBA), hydrogen bond donors (HBD), molecular weight (MW), topological polar surface area (TPSA), fraction of sp3 carbon atoms (FSP3), and number of unspecified chiral centers. The descriptors can be filtered using the following flags:
``--ha``, ``--logp``, ``--hba``, ``--hbd``, ``--mw``, ``--tpsa``, ``--fsp3``, and ``--chiral``. FSP3 is expressed as a fraction from 0 to 1. The filtering can be done using the following formats:

* Range: Specify a range using two values (e.g., "17-25").
* Greater / Less than or equal to: Use >= or <= (e.g., ">=17", "<=25").
* Greater than / Less than: Use > or < (e.g., ">17", "<25").
* Exact match: Match a specific value (e.g., 17).
* For logP, the exact match format applies as 'less than or equal to'.

For example, to filter the logP values less than or equal to 3.5, use the following command:

.. code-block:: console

    $ msani -i example.smi --logp "<=3.5"

For example, to retain molecules with TPSA up to 100 and FSP3 of at least 0.25:

.. code-block:: console

    $ msani -i example.smi --tpsa "<=100" --fsp3 ">=0.25"


6. PAINS filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Molecules that contain PAINS substructures can be efficiently eliminated using the ``--pains`` flag. The violated structures will be stored in the **_rejected** file.

.. code-block:: console

    $ msani -i example.smi --pains

Example of the **_rejected** output is as below:

.. code-block:: text

    CCOc1cccc(C=C2C(=O)N(Cc3ccccc3)C(C)=C2C(=O)OC)c1O Z57339064     "PAINS violation: Ene_five_het_c(85)"
    N#Cc1ccccc1COC(=O)c1cccc2c1C(=O)c1ccccc1C2=O      Z18301252     "PAINS violation: Quinone_a(370)"
    Nc1sc2c(c1C(=O)c1ccccc1)CCC2                      Z1259205366   "PAINS violation: Thiophene_amino_aa(45)"
    COCC1(CC(=O)NCc2cc(O)ccc2O)CC1                    Z2832180283   "PAINS violation: Mannich_a(296)"
    CCCCN(Cc1ccc(OS(=O)(=O)F)cc1)Cc1ccccc1O           Z4607533150   "PAINS violation: Mannich_a(296)"

7. Unwanted substructures filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Molecules that contain unwanted substructures can be efficiently eliminated using the ``--unwanted`` flag. msani uses an expert-curated list that contains undesirable substructures, accompanied by the reasons and references for filtering. The list can be obtained from `msani/Data/filter_out.csv <https://github.com/phonglam3103/msani/blob/main/msani/Data/filter_out.csv>`_.

There are four options accompanied by the ``--unwanted`` flag, which are *['all', 'regular', 'special', 'optional']*. If no option is specified, the *regular* filters will be applied. The choice of the options depends on the user and can vary between targets.

.. code-block:: console

    $ msani -i example.smi --unwanted
    $ msani -i example.smi --unwanted regular  # By default
    $ msani -i example.smi --unwanted regular special
    $ msani -i example.smi --unwanted all

It is also possible to filter out customized unwanted substructures, depending on the user's preference, using a customized SMARTS list. To generate a template for this list, use the ``--create_custom`` flag. This will result in the **templates.txt** file.

.. code-block:: console

    $ msani --create_custom

The first two columns (SMARTS and LABEL) are required for the program to parse, while the remaining columns will be omitted by the program. To filter using the customized list, use the ``--custom`` flag with the path to the customized list file. It is also possible to apply both the available filters and the customized filters.

.. code-block:: console

    $ msani -i example.smi --custom templates.txt
    $ msani -i example.smi --unwanted all --custom templates.tsv

8. Protonation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

msani supports the assignment of protonation states at various pH values using the ``--protonation`` flag. By default, the pH is set to 7 (configurable via ``-p`` or ``--pH``), and the pH range is set to 0 (specified using ``-r`` or ``--range``). This configuration protonates molecules at a specific pH of 7. However, it is also possible to enumerate potential protonation states across a pH range. For instance, setting ``--range 2`` explores pH values within 7 ± 2. The program evaluates each pH value in the specified range and assigns the possible protonation states of the molecule at those pH levels. Only unique products are output to a file. Functional groups with multiple protonation possibilities (e.g., piperazine, amidine) are expanded, with an underscore (`_`) appended to their names to indicate variations.

The program employs SMARTS-based reactions to iteratively assign protonation states to atoms, considering the pKa of functional groups and the queried pH. Detailed SMARTS reaction definitions are available in the following resource: `msani/Data/ionizations.txt <https://github.com/phonglam3103/msani/blob/main/msani/Data/ionizations_v2.txt>`_.

.. code-block:: console

    $ msani -i example.smi --protonation # Default pH 7 +- 0
    $ msani -i example.smi --protonation --pH 7 --range 2 # Enumerate protonation states at pH 7 +- 2
    $ msani -i example.smi --protonation -p 7 -r 2 # Short version


.. code-block:: text

   Input:
   O=C(N1C(C2C(C1)C2O)C(O)=O)CN3CCNCC3 mol4

   Output:
   O=C([O-])C1C2C(O)C2CN1C(=O)C[NH+]1CCNCC1 mol4_1
   O=C([O-])C1C2C(O)C2CN1C(=O)CN1CC[NH2+]CC1 mol4_2


9. Stereoisomers enumeration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Stereoisomers enumeration will be considered for unspecified chiral centers using the ``--stereoisomers`` flag. For an entry that contains multiple stereoisomers, its ID will be expanded (e.g., mol8 -> mol8.1, mol8.2).

.. code-block:: console

    $ msani -i example.smi --stereoisomers

.. code-block:: text

   Input:
   C1C2CC3CC1CC(C2)(C3O)N                            mol8

   Output:
   N[C@@]12C[C@@H]3C[C@@H](C[C@@H](C3)[C@H]1O)C2     mol8.1
   N[C@@]12C[C@@H]3C[C@@H](C[C@@H](C3)[C@@H]1O)C2    mol8.2

It is possible to define the maximum number of stereoisomers generated for each molecule by adding the ``--max_stereoisomers`` flag.

.. code-block:: console

    $ msani -i example.smi --stereoisomers --max_stereoisomers 32

10. Conformer generator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Basic usage
^^^^^^^^^^^^^^


The following supported flags:

.. code-block:: console

    Generate 3D conformers options:
    --gen3d, -3d          Generate 3D conformers
    --format, -f          Output file format. Multiple formats simultaneously supported.
                            (Default: db2.tgz - Options: sdf, db2, db2.tgz, mol2, pdbqt, oeb, oeb.lib.)
                            oeb: one .oeb.gz per molecule; oeb.lib: one library per input file.
                            Both need the optional OpenEye toolkits and OE_LICENSE.
    --method, -m          Embedding method (default: rdkit - options: rdkit, obabel, corina)
    --numconfs, -nconfs   Maximum number of conformers to generate (default: 2000)
    --timeout, -to        Timeout for the initial embedding for each entry before using OpenBabel
                            Default: 2 minutes
    --energywindow, -w    Energy window for sampling the conformations (default: 25 kcal/mol)
    --nringconfs, -nr     Maximum number of ring conformers to generate (default: 1)
    --mode, -mode         Mode for generating conformers
                            Default: fixed - Options: fixed, random, ignoretorlib
    --allowNonring        Allow the full sampling of non-ring comdpounds (default undersample to 30 confs).

The conformer generator platform can be triggered using the ``--gen3d`` or ``-3d`` flag. Three initial embeeder are supported (``-m`` or ``--method`` flag):

* RDKit srETKDG-v3 (default ``-m rdkit``): `Ref <https://pubs.acs.org/doi/10.1021/acs.jcim.0c00025>`__ 
* CORINA (``-m corina``): `Ref <https://doi.org/10.1016/0898-5529(90)90156-3>`__
* Open Babel (``-m obabel``): `Ref <https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0372-5>`__ 

Multiple aliphatic ring conformations are supported for RDKit and CORINA with the ``-nr`` flag.

As RDKit ETKDGv3 is based on distance geometry method, it may takes a long time to generate the initial conformer for some large molecules. In this case, the program will opt for the OpenBabel method, after the timeout (default: 2 minutes) is reached. The timeout can be modified using the ``--timeout`` flag. 

A modified version of `TorsionLibrary v3 <https://pubs.acs.org/doi/10.1021/acs.jcim.2c00043>`_ is used to drive the generation of conformations. The modifications made and the full library can be obtained `here <https://github.com/phonglam3103/msani/blob/main/msani/Data/modified_tor_lib_2020.xml>`_. The number of conformers are controlled by the ``--numconfs`` or ``-nconfs`` flag. The default value is 2000, but it can be modified to any number. The program will sample the conformers based on the energy window (default: 25 kcal/mol) using the ``--energywindow`` or ``-w`` flag.

Three sampling modes are supported (``--mode`` or ``-mode`` flag):

* `fixed (default)`: the peaks (with tolerance) in the TorsionLibrary is discretinized into central angles with +-30 degrees offset (as long as they are within the tolernace 2), then combinatorially sampled to generate the conformers until the number of conformers is reached. Symmetric substructures (such as phenyl, carboxylates, etc.) are removed by SMARTS matching in advance.
* `random`: each peak combination is sampled multiple times. Two conformers are regarded distinct if they differ by at least 30 degrees in any dihedral angle.
* `ignoretorlib`: the program will ignore the TorsionLibrary and sample every 30 degrees.


Supported output formats are ``db2``, ``db2.tgz``, ``pdbqt``, ``sdf``, ``mol2``,
``oeb``, and ``oeb.lib``. The default is the compressed DB2 archive
``db2.tgz``. Select one or more formats with ``--format`` or ``-f``.
PDBQT requires the optional Meeko dependencies; both OEB formats require
the optional OpenEye toolkits and a valid licence (see :doc:`installation`).
See :doc:`outputs` for OEB record grouping, file naming, and restart behavior.

For DB2 generation, AMSOLcpp assigns desolvation penalties and partial charges
directly from the in-memory RDKit molecule. These values are passed directly to
the DB2 layout writer without intermediate AMSOL or ``.solv`` files.

.. code-block:: console

    $ msani -i example.smi --protonation --stereoisomers -3d -f db2  # Generate DB2 files
    $ msani -i example.smi --tautomers --protonation --stereoisomers -3d -f sdf pdbqt #Generate SDF and PDBQT files

Customization of the torsion definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to add and/or modify the torsion definitions in the TorsionLibrary without modifying the Torlib file by the ``--torsion`` flag. The flag accepts a file which defines multiple SMARTS patterns definining rotatable bonds with the expected dihedral angles. To create a template for the torsion definition file, use the ``--create_torsion`` flag. The template will be saved in the **custom_torsion_templates.txt** file:

.. code-block:: console

    $ msani --create_torsion

The template file will look like this:

.. code-block:: text

    # Template for customized torsion driving
    # Format: Comma separated. Weight is optional, not set = equal weights
    # In case using commas in SMARTS, put it in quotation marks ("text") 
    # SMARTS, possible angle values (sep by space), weight (optional)
    # Example:
    # [O:1]=[C:2]-[N:3][H:4], 180
    # "[O:1]=[C:2]-[N,O,X:3][!H:4]", 0 180
    # [*:1]~[CX4:2]!@[OX2:3]~[*:4], 0 120 240, 1 1 1

To apply the torsion definition file, use the ``--torsion`` flag with the path to the file:

.. code-block:: console

    $ msani -i example.smi --torsion custom_torsion_templates.txt -3d -f db2

Advanced options
-----------------------------------

The advanced options can be accessed using the ``--help_advanced`` or ``-xh`` flag. The advanced options are not recommended for general users, but they can be useful for advanced users who want to customize the program's behavior. The following advanced options are available:

.. code-block:: console

    --noneutralize          Do not neutralize the molecule before tautomerization and protonation
    --notaurdkit            Do not use RDKit to canonicalize the tautomeric form of the input SMILES
    --rigid, -r             Only align the DB2 on this rigid scaffold in SMARTS format. All rings if not provided.
    --tolerance, -tol       Minimum angle for differentiating two conformers (default: 30)
    --nringconfs, -nr       Maximum number of ring conformers to generate (default: 1)
    --allowNonring          Allow the full sampling of non-ring compounds (default undersample to 30 confs).
    --eps                   The dielectric constant for electrostatic calculations (default: 4).
    --debug                 Enable debug mode
    --timing                Enable timing information
    --create_protlib      Create a template for customized protonation scheme
    --create_taulib       Create a template for customized tautomerization scheme
    --create_torsion      Create a template for customized torsion definition
    --protlib             Path to the protonation library file (default: msani/Data/ionizations_v3.txt).
    --taulib              Path to the tautomer library file (default:  msani/Data/tautomers_v3.txt).
    --torsion, -tor       Path to the customized torsion definitions.


Running in batch mode
-----------------------------------


msani now supports the batch mode ``msani_batch``, which allows handling bigger SMILES databases on a SLURM-based cluster. Nearly all flags supported by standalone ``msani`` are available in batch mode. The ``-l`` or ``--lines_per_job`` option (default: 200) defines how many input records each array job processes. By default, a maximum of 500 jobs run simultaneously, but this can be changed with ``--max_jobs``.

Batch mode accepts plain, gzip-compressed, bzip2-compressed, and xz-compressed
inputs. Each input gets a corresponding batch directory containing
``input.data``, ``input.offsets``, and the submission script. For plain input,
``input.data`` is a symbolic link to the original file. Compressed input is
decompressed once into the single ``input.data`` staging file. During that same
pass, MolSanitizer records one byte boundary per job in ``input.offsets``.
Array tasks use these offsets to copy only their assigned range to node-local
scratch space, so the shared filesystem does not contain thousands of input
chunks. The temporary task input is removed when the task exits.

Each input is validated and submitted independently. An input that exceeds the
array-size or currently available project capacity is skipped without
preventing later inputs from being considered. With cleanup enabled, shared
staging data and its offset index are removed after every task succeeds; they
remain available when tasks must be retried. Successful output shards are
sorted by their numeric job identifier and atomically concatenated into
``in/processed/processed.smi``. Rejected-output shards are similarly merged into
``in/removed/removed.smi``. Individual shards are deleted only after their merge
succeeds.

The additional flags supported by ``msani_batch`` so far:

.. code-block:: console

    --projectName, -A           The account that will be charged by the SLURM cluster for running tasks (default: PROJECT_NAME)
    --lines_per_job, -l         Number of lines to process per job (default: 200)
    --timelimit, -tl            Time limit in hours for each SLURM job (default: 96)
    --max_jobs, -mj             Maximum number of jobs to run simultaneously (default: 500)

Batch examples
~~~~~~~~~~~~~~

.. code-block:: console

    $ msani_batch -i example.smi -l 50 -3d -f db2
    $ msani_batch -i example.smi -l 50 --stereosiomers --protonation -3d -f db2 --nocleanup
    $ msani_batch -i example.smi -l 50 -A snic2021-3-32 -tl 2 -3d -f db2

It is also possible to submit the batch jobs for multiple input files. The program will automatically detect the input files and submit the jobs accordingly.

.. code-block:: console

    $ msani_batch -i example.smi example2.smi -3d -f db2 --protonation --stereoisomers


Structure standardization
-------------------------

Use ``--standardize`` to obtain a parent representation for applications such
as descriptor calculation or machine learning:

.. code-block:: console

    $ msani -i example.smi --standardize

After invalid and unsupported structures are removed, the standardization
branch applies RDKit ``Cleanup``, selects the ``FragmentParent``, and applies
``Uncharger`` before writing canonical SMILES. It also excludes entries
containing the metal elements listed by the standardization filter.
Neutralization is attempted where possible; permanent charges can remain.
This does not assign protonation states at a specified pH.

In the current implementation, explicitly requested salt removal and descriptor
filters still run before standardization. The standardization branch then
returns before the separate neutralizer, tautomerizer, PAINS/custom/unwanted
filters, ionizer, and stereoisomer enumerator. Consequently, those preparation
options do not run in combination with ``--standardize``. Requested ``--gen3d``
output is still generated afterwards from the standardized structures.
For standardization alone, use the command above without other processing flags.

The output follows the usual ``_clean`` naming convention. See :doc:`outputs`
for output files, failures, and reruns.
