Outputs
========

Files and identifiers
---------------------

Outputs are written relative to the working directory. For a plain input named
``example.smi``, the main outputs are:

.. list-table:: Output files
   :header-rows: 1
   :widths: 35 65

   * - File or directory
     - Contents
   * - ``example_clean.smi``
     - Prepared SMILES and identifiers, without a header.
   * - ``example_rejected.smi``
     - Rejections recorded by enabled filters, including reasons. This file
       need not exist when no rejections are recorded.
   * - ``sdf/``, ``mol2/``, ``pdbqt/``, ``db2/``
     - Requested 3D outputs, normally named by the prepared molecule ID.
       Multiple ring conformations may add an ``.nr`` suffix.
   * - ``db2/example.db2.tgz``
     - DB2 archive when ``--format db2.tgz`` is selected (the default 3D format).
   * - ``oeb/<molecule ID>.oeb.gz``
     - Per-molecule multi-conformer output with ``--format oeb``.
   * - ``oeb/example.oeb.gz``
     - Multi-molecule library with ``--format oeb.lib``.
   * - ``msani.log``
     - Execution log; ``--prefix run1`` changes the log name to ``run1.log``.
   * - ``msani_error.err``
     - Molecule-level conformer failures, with the molecule ID and reason.

``--prefix run1`` names the prepared and rejected files ``run1_clean.txt`` and
``run1_rejected.txt``. It does not rename molecule-based 3D output files.
Without a prefix, naming uses the input's final suffix: for example,
``example.smi.gz`` produces ``example.smi_clean.gz``. The prepared text writer
uses pandas' filename-based compression inference. An explicit prefix avoids
surprising output names for compressed inputs.

Tautomer/protomer expansion appends underscore indices, for example ``mol_1``
and ``mol_2``. Stereoisomer enumeration appends dot indices, such as ``mol_1.1``
and ``mol_1.2``. Larger expansions use zero-padded indices. Use unique input IDs
so per-molecule files do not collide, and keep the input to trace derived IDs.

OpenEye records and naming
--------------------------

Both OEB formats export sampled conformers. Ring-conformer groups are normally
merged into one multi-conformer ``OEMol`` per prepared molecule, titled with
its identifier, rather than separate records for each conformer. Groups that
fail the writer's atom-count or atomic-number-order checks are retained as
separate records. Distinct prepared tautomers, protomers, and stereoisomers
remain separate molecules.

The library filename uses the input basename with its final extension removed:
``example.smi`` produces ``oeb/example.oeb.gz`` and ``example.smi.gz`` produces
``oeb/example.smi.oeb.gz``. ``--prefix`` does not rename this library.
When requesting both ``oeb`` and ``oeb.lib``, ensure no prepared molecule ID
matches that input basename: both formats use the same directory and extension,
so their files would collide. Use separate working directories if necessary.

Rejected and failed molecules
-----------------------------

A filter rejection means a molecule failed an enabled selection criterion.
A conformer failure means preparation could not complete a requested 3D step.
A molecule can appear in the clean SMILES output while lacking a 3D file;
check ``msani_error.err`` and the execution log. The rejected file is not a
complete ledger of every removed input: invalid SMILES and some other early
removal paths are reported separately in the log.

See :doc:`troubleshooting` for DB2 identifier limits, conformer clashes,
missing outputs, and recovery of interrupted runs or failed batch tasks.
