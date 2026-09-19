Outputs and troubleshooting
===========================

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

Rejected and failed molecules
-----------------------------

A filter rejection means a molecule failed an enabled selection criterion.
A conformer failure means preparation could not complete a requested 3D step.
A molecule can appear in the clean SMILES output while lacking a 3D file;
check ``msani_error.err`` and the execution log. The rejected file is not a
complete ledger of every removed input: invalid SMILES and some other early
removal paths are reported separately in the log.

Common problems
---------------

* **Invalid SMILES:** check the input syntax, the two-column layout, and the
  absence of a header. Extended SMILES require tab-separated input and
  ``--extended``.
* **No molecules remain:** inspect the rejection reasons and enabled filters.
  ``--lazy`` enables all unwanted-filter categories; choose narrower options
  when appropriate for your project.
* **Missing PDBQT output:** install ``molsanitizer[pdbqt]`` in the environment
  running the command and inspect the molecule's error message.
* **Embedding timeout:** ``--timeout`` is measured in minutes (default 2).
  When RDKit embedding times out, MolSanitizer tries Open Babel. If Open Babel
  is unavailable or the fallback fails, the reason is recorded and processing
  continues with other molecules. ``--timeout_conf`` separately limits
  torsional sampling; consult ``--help_advanced`` for its current setting.
* **Fewer conformers than requested:** ``--numconfs`` is a maximum. Symmetry,
  molecular rigidity, RMSD pruning, energy/clash filters, and time limits can
  all reduce the ensemble size.

Reruns and recovery
-------------------

A file-processing rerun replaces the matching clean and rejected text outputs.
Per-molecule files may be overwritten, while files for molecules absent from
the new run can remain. Use a fresh working directory when changing chemistry
settings or inputs; changing only ``--prefix`` does not isolate the 3D outputs.

For an existing ``db2.tgz`` archive, MolSanitizer attempts to restore completed
DB2 members and skip their molecule IDs. This supports restarting interrupted
work with the same inputs and settings. An unreadable archive is logged and
processing starts fresh. Use a separate directory for a new preparation rather
than reusing an old archive with changed settings.

SLURM batch recovery
--------------------

Each batch directory contains ``submit_msani.sh``, the defaults snapshot
``msani_defaults.yaml``, and staged input files ``input.data`` and
``input.offsets``. Keep these together while jobs are unfinished. For plain
input, ``input.data`` links to the original file, so keep that file available
and unchanged as well.

Inspect the SLURM task logs and scheduler state to identify failed array task
IDs. After correcting the cause, resubmit only those tasks from the batch
directory, for example:

.. code-block:: console

    $ sbatch --array=3,7 submit_msani.sh

Here ``3,7`` are examples; replace them with your failed task IDs. Check that
the staging files still exist before resubmitting. Do not assume a scheduler
success means every molecule produced every requested format: molecule-level
failures are recorded separately.

With cleanup enabled, successful shards are merged into
``in/processed/processed.smi`` and ``in/removed/removed.smi``. Staging data is
retained while tasks remain unsuccessful, and shards are deleted only after
their merge succeeds. ``--nocleanup`` retains the intermediate files for
inspection. If staging data has already been removed, prepare a new batch
from the original input instead of resubmitting the old script.
