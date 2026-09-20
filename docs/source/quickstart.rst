Quickstart
==========

Install MolSanitizer following :doc:`installation`. Work in a new directory
so the example outputs are easy to identify. Create a file named ``example.smi``
with these two lines, without a header:

.. code-block:: text

    CC(=O)O acetate
    c1ccccc1 phenyl

Prepare SMILES
--------------

.. code-block:: console

    $ msani -i example.smi --removesalts --tautomers --protonation --pH 7

With bundled defaults, ``example_clean.smi`` contains the prepared SMILES and
original identifiers. Acetic acid is deprotonated and benzene remains neutral:

.. code-block:: text

    CC(=O)[O-] acetate
    c1ccccc1 phenyl

Personal defaults can affect a run; inspect them with ``msani config show``.
For a project, save your selected options in a job configuration as described
in :doc:`config`.

Generate SDF conformers
-----------------------

.. code-block:: console

    $ msani -i example_clean.smi --gen3d --format sdf --numconfs 20

Successful generation creates ``sdf/acetate.sdf`` and ``sdf/phenyl.sdf``.
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

Successful generation creates ``pdbqt/acetate.pdbqt`` and
``pdbqt/phenyl.pdbqt``. PDBQT is prepared with Meeko from the initial embedded
structure, with rotatable bonds for docking; it does not export the sampled
SDF ensemble. To request both formats in one run, use ``--format sdf pdbqt``.

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

See :doc:`usage` for filtering and other processing options, and :doc:`outputs`
for logs, output naming, and recovery when a molecule fails.
