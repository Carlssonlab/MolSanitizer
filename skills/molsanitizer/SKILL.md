---
name: molsanitizer
description: Prepare, filter, protonate, and generate small-molecule conformers or docking libraries with MolSanitizer (msani), using its CLI or Python API. Use for MolSanitizer workflows, output inspection, and failed molecule recovery.
---

# MolSanitizer

Use the installed MolSanitizer CLI for file workflows and its Python API when
integrating preparation into a script. Preserve the requested chemistry and
record the settings used; preparation can expand or remove input records.

## Establish the run

- Inspect `msani --version`, `msani --help`, and `msani config show` in the
  environment that will run the job. Use `msani --help_advanced` for sampling
  controls. Check flags against that version before using them.
- If installation is needed, prefer `python -m pip install molsanitizer` in
  the intended environment. Native source builds are a development workflow.
  PDBQT needs `molsanitizer[pdbqt]`; OEB needs `molsanitizer[oe]` and a valid
  OpenEye licence selected through `OE_LICENSE`.
- Determine the input, output format, preparation steps, pH/range when
  protonating, and available compute from the task. Ask about missing choices
  that change the intended chemistry; label assumptions in example workflows.
- Use a fresh working directory for new inputs or changed settings. Outputs
  are relative to the working directory. A new `--prefix` does not isolate
  molecule-based 3D files. Resume in an existing directory only with the same
  inputs and settings.

## Input and chemistry

File input has two columns, `SMILES ID`, without a header. Use unique IDs that
are safe as filenames and keep the original input for tracing derived IDs.
Gzip, bzip2, and xz input are recognized by content. Unquoted CXSMILES requires
tabs between structures and IDs throughout the file; alternatively quote the
entire CXSMILES field for space-separated input. `--extended` forces tabs.

Choose preparation steps explicitly:

- `--removesalts`, `--tautomers`, and `--protonation` enable the respective
  transformations. The CLI also enables neutralization with these steps
  unless `--no-neutralize` is supplied.
- Specify both `--pH` and `--pH_range` for reproducible protonation intent.
  Since 0.8.0 the bundled range is 1 (pH 7 spans 6–8), previously 0.
  Explicit pH options enable protonation; inherited pH defaults alone do not.
- `--gen3d` enables unspecified stereochemistry enumeration by default.
  Use `--no-stereoisomers` when the task requires disabling enumeration.
- Add PAINS, unwanted-substructure, or descriptor filters only when the task
  calls for those selection criteria. `--lazy` also enables PAINS and all
  unwanted-filter categories, so it can discard compounds unexpectedly.
  Quote descriptor comparisons, for example `--logp "<=3.5"`.
- `--standardize` is a separate predictor-preparation mode that skips the
  other processing steps; do not treat it as an additional docking step.

Configuration precedence is explicit CLI arguments, then job YAML (`--config`),
then personal defaults, then bundled defaults. `msani config show` describes
personal/bundled settings, not the complete effective settings of a job.
Use CLI arguments or job YAML for project choices without changing global
personal defaults. Job YAML uses CLI names such as `pH_range`; personal
configuration uses uppercase names such as `PH_RANGE`.

## CLI examples

These examples assume a headerless `input.smi`, salt removal and tautomer
preparation, and protonation around pH 7 ± 1. Adapt these choices to the task.
Run in a new output directory with the input copied or referenced by its path.

```bash
msani -i input.smi --removesalts --tautomers --protonation --pH 7 --pH_range 1
msani -i input_clean.smi --gen3d --format sdf --numconfs 20 --numcores 2
```

The first command writes `input_clean.smi`; the second normally writes
`sdf/<prepared-ID>.sdf` and `input_clean_clean.smi`. Twenty conformers is an
example limit for a small run, not a production sampling recommendation.
`--numconfs` is a maximum, not a guaranteed record count. For a large library,
check a representative sample with the intended settings before scaling up.

Select formats explicitly; the default 3D format is `db2.tgz`:

| Requested artifact | CLI format | Important behavior |
| --- | --- | --- |
| SDF ensemble | `--format sdf` | Sampled conformers; multiple records per molecule. |
| MOL2 ensemble | `--format mol2` | Sampled conformers. |
| AutoDock Vina ligand | `--format pdbqt` | Meeko uses the initial embedded structure with rotatable bonds, not the sampled SDF ensemble. |
| DOCK ligand/library | `--format db2` or `--format db2.tgz` | DB2 files or an archive; internal names keep only the last 16 ID characters. |
| OpenEye molecules | `--format oeb` | Per-molecule multi-conformer files in `oeb/`. |
| FRED/HYBRID library | `--format oeb.lib` | Multi-molecule library at `oeb/<input-basename>.oeb.gz`. |

Multiple formats can be requested together, for example `--format sdf pdbqt`.
For DB2, verify uniqueness of the last 16 characters of **prepared** IDs,
including enumeration suffixes, and retain a mapping to full IDs. For combined
`oeb` and `oeb.lib`, a molecule ID matching the input basename would cause a
filename collision; use separate directories or formats in that case.

## Python integration

Use `Msani.run()` for a DataFrame with `smiles` and `ids` columns. Its constructor
does not load personal CLI defaults. Set the chemistry explicitly, including
`neutralize` when matching the CLI example above. For example:

```python
import pandas as pd
from msani.api import Msani


def main():
    data = pd.DataFrame({"smiles": ["CC(=O)O"], "ids": ["acetic_acid"]})
    with Msani(removesalts=True, neutralize=True, tautomers=True,
               protonation=True, pH=7, pH_range=1, numcores=2) as processor:
        prepared = processor.run(data)
    prepared[["smiles", "ids"]].to_csv(
        "prepared.smi", sep=" ", header=False, index=False)


if __name__ == "__main__":
    main()
```

Save worker-starting code in an importable script with the main guard; worker
startup uses `forkserver` on Linux and `spawn` on macOS/Windows. Reuse a `Msani`
context manager across input chunks to share workers. Do not share a processor
between concurrent threads. `Msani.run()` prepares structures; for Python 3D
generation, consult the `ConformerGenerator` examples and API documentation.

## Verify and recover

- Inspect the execution log (normally `msani.log`), rejected output when
  present, `msani_error.err`, and the requested format files. A successful
  exit or a clean SMILES row does not guarantee successful 3D output.
- Rejected files cover enabled filters, not all losses: invalid input and
  other early removals can be logged separately. Inspect any additional
  failure files reported by the installed version.
- Tautomer/protomer expansion adds underscore indices; stereoisomers add dot
  indices. Compare input coverage through derived IDs instead of assuming
  equal input/output row counts or counting conformers as distinct compounds.
- Existing `db2.tgz` and `oeb.lib` outputs can reuse completed IDs on a rerun.
  With both formats requested, IDs restored from DB2 skip all processing and
  may leave missing OEB records unrepaired. Use a fresh directory to regenerate
  inconsistent outputs. See troubleshooting before retrying failed work.
- Report the command or script, versions, configuration and chemistry choices,
  output paths, prepared structures/coverage, rejections, and unresolved
  failures. Distinguish verified counts from estimates.

For SLURM tasks, read the batch usage and recovery documentation first.
`msani_batch` immediately submits via `sbatch`; invoke it only when cluster
submission is within the user's request. Derive account, partition and limits
from the task/environment. Keep `submit_msani.sh`, `msani_defaults.yaml`, staged
input files, and the original input together until jobs finish. Scheduler
success still requires molecule-level output checks.

## Detailed documentation

Prefer documentation matching the installed release. In a source checkout,
the corresponding sources are `docs/source/*.rst`; `latest` documentation can
describe unreleased changes.

- [Installation](https://msani.readthedocs.io/en/stable/installation.html)
- [CLI options and SLURM](https://msani.readthedocs.io/en/stable/usage.html)
- [Configuration](https://msani.readthedocs.io/en/stable/config.html)
- [Output naming](https://msani.readthedocs.io/en/stable/outputs.html)
- [Failures and recovery](https://msani.readthedocs.io/en/stable/troubleshooting.html)
- [Python examples](https://msani.readthedocs.io/en/stable/python.html)
- [API reference](https://msani.readthedocs.io/en/stable/api.html)
