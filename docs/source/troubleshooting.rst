Troubleshooting
===============

See :doc:`outputs` for output filenames, identifiers, logs, and rejection files.

.. _db2-name-limit:

DB2 molecule names are limited to 16 characters
-----------------------------------------------

DOCK reads the DB2 molecule name from a fixed 16-character field. MolSanitizer
therefore writes only the **last 16 characters** of a longer prepared molecule
ID into that field. The full prepared ID remains in the DB2 filename, including
archive member names, but docking results identify the molecule by its short
internal name.

For example, ``ZINC000000198632_1`` becomes ``NC000000198632_1`` inside the
DB2 file. ``OUTDOCK``, poses, and ``top_poses.scores`` use that short name.
Keep a mapping from ``full_id[-16:]`` to the full prepared ID, using DB2
filenames to recover it. Account for suffixes added by tautomer, protomer, and
stereoisomer enumeration when building this mapping.

Two distinct IDs with the same last 16 characters become indistinguishable
in docking results. Give compounds short, unique IDs before preparation, keep
a mapping to their original names, and check that the **prepared IDs** remain
unique after truncation. Do not rely on differences at the beginning of a long
name. Changing ``--prefix`` does not change molecule IDs or this format limit.

.. _conformer-clashes:

MOE or Maestro reports clashes in generated conformers
------------------------------------------------------

Different programs use different clash definitions. Passing MolSanitizer's
filter does not guarantee that another program will accept every contact.

MolSanitizer's torsional sampler rejects a conformer if a checked non-bonded
atom pair satisfies:

.. math::

    d_{ij} < s \left(r_i^{\mathrm{vdW}} + r_j^{\mathrm{vdW}}\right)

Here ``s`` is ``--clash_scale`` (``-clash``), which defaults to **0.6**.
The radii come from RDKit's periodic table. Directly bonded atom pairs and
pairs sharing a bonded neighbor (1-2 and 1-3 interactions) are excluded.

For comparison, Maestro's contact display classifies a distance below
**0.75** times the sum of van der Waals radii as an **ugly** contact, and uses
**0.89** as the upper cutoff for **bad** contacts; see the
`Schrödinger clash-detection documentation
<https://learn.schrodinger.com/public/python_api/2024-1/api/schrodinger.comparison.struc.html#schrodinger.comparison.struc.has_clash>`_.
A distance ratio between 0.6 and 0.75 can therefore pass MolSanitizer's distance
filter and still be flagged as ugly in Maestro. Radii and atom-pair exclusions
can differ too, so matching the scale alone does not ensure identical results.

MOE can also flag contacts in MolSanitizer output. Check the criterion and
settings of the MOE tool you are using. If it uses a force-field energy cutoff
for clashes, that energy threshold has no direct conversion to MolSanitizer's
distance scale; do not assume one universal MOE cutoff.

To reject more close contacts during sampling, **increase** the MolSanitizer
clash scale, for example from 0.6 to 0.75:

.. code-block:: console

    $ msani -i example.smi --gen3d --format sdf --clash_scale 0.75
    $ msani_batch -i example.smi --gen3d --format db2.tgz --clash_scale 0.75

Use a fresh working directory when comparing settings, especially for
``db2.tgz`` and ``oeb.lib``, whose restart logic can reuse completed molecules.
Inspect the regenerated ensemble in the downstream program to select a value
appropriate for your workflow. See :ref:`conformer-clash-cutoff` for job YAML
and related sampling controls.

**Tradeoff:** a larger scale rejects more conformers and can leave fewer
conformers per compound, or no accepted conformer for some compounds. It can
also reduce the fraction of compounds with a conformer close to an experimental
or other reference structure. A stricter cutoff may discard a useful strained
or compact reference-like conformation. Compare both clash counts and reference
coverage before adopting it for a large database.

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
