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

Common problems
---------------

* **Invalid SMILES:** check the input syntax, the two-column layout, and the
  absence of a header. Use tab-separated input for unquoted CXSMILES, or quote
  the entire CXSMILES field in double quotes for space-separated input.
  ``--extended`` forces tab-separated parsing; omit it for space-separated input.
* **No molecules remain:** inspect the rejection reasons and enabled filters.
  ``--lazy`` enables all unwanted-filter categories; choose narrower options
  when appropriate for your project.
* **Missing PDBQT output:** install ``molsanitizer[pdbqt]`` in the environment
  running the command and inspect the molecule's error message.
* **Missing OEB output:** install ``molsanitizer[oe]`` in the environment
  running the command and set ``OE_LICENSE`` to a valid licence file. Missing
  toolkits and an unlicensed OEChem installation produce distinct errors.
  Inspect the execution log and ``msani_error.err`` for details.
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

For file input with ``oeb.lib``, an existing library is moved aside and its
records are copied back before new molecules are added. This retains records
across input chunks and supports reruns with the same input and settings.
Restored molecule titles are used to avoid adding those molecules again; use
unique IDs. With only ``oeb.lib`` selected, restored molecules skip conformer
generation. With other formats selected, they may still be processed for those
outputs. This is not a guarantee of recovery from a corrupt or partially
written library; keep a backup of valuable output before retrying.

When ``db2.tgz`` and ``oeb.lib`` are requested together, IDs already restored
from the DB2 archive skip all processing, including OEB output. A rerun therefore
does not repair an OEB library missing those IDs. Regenerate in a fresh directory
if the two outputs are inconsistent.

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
