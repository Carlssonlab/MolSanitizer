Quickstart
==========

Install MolSanitizer following :doc:`installation`. Work in a new directory
so the example outputs are easy to identify. Create a file named ``example.smi``
with these two lines, without a header (both tab and space are accepted as separators):

.. code-block:: text

    Nc1nc(C(=CCC(=O)O)C(=O)NC2C(=O)N3C(C(=O)O)=CCSC23)cs1 Ceftibuten
    CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O Ibuprofen
    C[C@]1(c2cccc(c2C(=O)C3=C([C@]4([C@@H](C[C@@H]31)[C@@H](C(=C(C4=O)C(=O)N)O)N(C)C)O)O)O)O Tetracycline
    O=C1NC(=O)NC(=O)C1(c2ccccc2)CC Phenobarbital

CXSMILES (extended SMILES) is recognized automatically. Unquoted CXSMILES needs
tab-separated structures and IDs throughout the file. Alternatively, enclose the
complete CXSMILES field in double quotes to use space-separated input:

.. code-block:: text

    "C[C@H](O)F |&1:1|" compound_001
    CCO ethanol


Prepare SMILES
--------------

.. code-block:: console

    $ msani -i example.smi --removesalts --tautomers --protonation --pH 7

With bundled defaults, ``example_clean.smi`` contains the prepared SMILES and
original identifiers. Acetic acid is deprotonated and benzene remains neutral:

.. code-block:: text

    Nc1nc(C(=CCC(=O)[O-])C(=O)NC2C(=O)N3C(C(=O)[O-])=CCSC23)cs1 Ceftibuten
    CC(C)Cc1ccc([C@@H](C)C(=O)[O-])cc1 Ibuprofen
    C[NH+](C)[C@@H]1C([O-])=C(C(N)=O)C(=O)[C@@]2(O)C(O)=C3C(=O)c4c(O)cccc4[C@@](C)(O)[C@H]3C[C@@H]12 Tetracycline_1
    C[NH+](C)[C@@H]1C([O-])=C(C(N)=O)C(=O)[C@@]2(O)C([O-])=C3C(=O)c4c(O)cccc4[C@@](C)(O)[C@H]3C[C@@H]12 Tetracycline_2
    CCC1(c2ccccc2)C(=O)NC(=O)NC1=O Phenobarbital_1
    CCC1(c2ccccc2)C(=O)[N-]C(=O)NC1=O Phenobarbital_2

Personal defaults can affect a run; inspect them with ``msani config show``.
To modify the default configurations or save options for a project, see
:doc:`config`.

Generate SDF conformers
-----------------------

.. code-block:: console

    $ msani -i example_clean.smi --gen3d --format sdf --numconfs 20

Successful generation creates ``sdf/Ibuprofen.sdf``, ``sdf/Ceftibuten.sdf``, etc.
Each SDF can contain multiple conformer records, up to the requested limit;
these rigid examples may produce fewer than 20. The command also writes
``example_clean_clean.smi``. Unspecified stereochemistry is enumerated by
default for 3D generation, which can expand identifiers on other inputs.

Generate PDBQT for docking
--------------------------

Install the optional dependencies, then prepare PDBQT files:

.. code-block:: console

    $ pip install "molsanitizer[pdbqt]"
    $ msani -i example_clean.smi --gen3d --format pdbqt

Successful generation creates pdbqt files within the ``pdbqt/`` directory.
PDBQT is prepared with Meeko from the initial embedded
structure, with rotatable bonds for docking; it does not export the sampled
SDF ensemble. To request both formats in one run, use ``--format sdf pdbqt``.

Generate DB2 for docking
--------------------------

No additional dependencies are required to prepare DB2 files for docking with DOCK3 and DOCK6. Run:

.. code-block:: console

    $ msani -i example_clean.smi --gen3d --format db2

Successful generation creates db2 files in the ``db2/`` directory, such as
``db2/Ceftibuten.db2`` and ``db2/Ibuprofen.db2``.

.. caution::

    DOCK reads the DB2 molecule name as a fixed 16-character field, so
    MolSanitizer, like the ZINC pipeline, keeps only the last 16 characters of
    longer IDs: `ZINC000000198632_1` becomes `NC000000198632_1`.
    `OUTDOCK`, the poses, and `top_poses.scores` show this short name. The
    full ID is kept only in the DB2 filename, which DOCK does not report, so map
    results back from there (`full_name[-16:] -> full_name`).

    If two IDs end with the same 16 characters, DOCK cannot tell them apart and
    `top_poses.py` keeps only the better-scoring one. IDs should therefore be
    unique in their last 16 characters; ZINC and ChEMBL IDs normally are.

Generate an OpenEye docking library
----------------------------------

Install the ``oe`` extra and configure ``OE_LICENSE`` as described in
:doc:`installation`, then run:

.. code-block:: console

    $ msani -i example_clean.smi --gen3d --format oeb.lib --numconfs 20

This writes ``oeb/example_clean.oeb.gz`` for use as a FRED/HYBRID docking
library. Each prepared molecule normally has one multi-conformer record,
combining its sampled ring-conformer groups. To write individual files such as
``oeb/acetate.oeb.gz`` instead, select ``--format oeb``. Both formats trigger
conformer sampling. You can export an SDF ensemble alongside the library with
``--format sdf oeb.lib``.

Submit SLURM batch jobs
----------------------

On a SLURM cluster, activate the environment containing MolSanitizer and use
``msani_batch`` to submit processing as a job array:

.. code-block:: console

    $ msani_batch -i example_clean.smi -l 200 -A my-project -tl 2 \
        --gen3d --format sdf --numconfs 20

Each task processes up to 200 input records (``-l``) with a two-hour time limit
(``-tl``). Replace ``my-project`` with your SLURM account; omit ``-A`` if your
cluster does not require one. Add ``--partition`` if you need a specific
partition. The command creates an ``example_clean/`` batch directory containing
``submit_msani.sh`` and submits it automatically with ``sbatch``. See
:doc:`usage` for more batch options and :doc:`outputs` for recovery of failed
tasks. The default batch configurations could be changed according to :doc:`config`.

See :doc:`usage` for filtering and other processing options, and :doc:`outputs`
for logs, output naming, and recovery when a molecule fails.
