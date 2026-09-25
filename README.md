# MolSanitizer

[![Python](https://img.shields.io/badge/python-3.10--3.14-blue)](pyproject.toml)
[![Documentation](https://img.shields.io/badge/docs-msani.readthedocs.io-orange)](https://msani.readthedocs.io/)
[![License](https://img.shields.io/badge/license-Apache%202.0-yellow)](LICENSE)

MolSanitizer (`msani`) prepares SMILES libraries for structure-based drug
discovery. It standardizes molecules (salt removal, tautomer and protonation
state enumeration, stereoisomer enumeration), filters undesirable substructures
(PAINS, reactive functional groups, and more), and generates 3D conformers in
formats ready for docking.

<img src="./docs/source/_static/Workflow.png" width="1000" alt="MolSanitizer workflow">

## Installation

```bash
pip install molsanitizer
```

Prebuilt wheels are available for Linux and Windows. For building from source
and details on dependencies, see the
[installation guide](https://msani.readthedocs.io/en/latest/installation.html).

### Optional dependencies

| Feature | Install | Notes |
| --- | --- | --- |
| PDBQT output | `pip install "molsanitizer[pdbqt]"` | Requires Meeko. |
| OEB output | `pip install "molsanitizer[oe]"` | Requires OpenEye toolkits and a licence (`export OE_LICENSE=/path/to/oe_license.txt`). |
| Open Babel embedding | `conda install openbabel` | Used when selected as the embedding method, or as a fallback after an RDKit embedding timeout. |

## Quick start

Create `example.smi` (one SMILES and ID per line, no header):

```text
CC(=O)O acetate
c1ccccc1 phenyl
```

Prepare the SMILES (strip salts, enumerate tautomers, protonate at pH 7):

```bash
msani -i example.smi --removesalts --tautomers --protonation --pH 7
```

This writes `example_clean.smi`, where acetic acid is deprotonated:

```text
CC(=O)[O-] acetate
c1ccccc1 phenyl
```

Generate 3D structures for docking from the prepared file:

```bash
msani -i example_clean.smi --gen3d --format sdf --numconfs 20   # SDF conformer ensembles
msani -i example_clean.smi --gen3d --format db2                 # DOCK3/DOCK6
msani -i example_clean.smi --gen3d --format pdbqt               # AutoDock Vina (needs [pdbqt])
msani -i example_clean.smi --gen3d --format oeb.lib             # OpenEye FRED/HYBRID (needs [oe])
```

Each format is written to its own directory (`sdf/`, `db2/`, `pdbqt/`, `oeb/`).
Several formats can be requested at once, e.g. `--format sdf db2`. Run
`msani -h` for all options, and see the
[quickstart](https://msani.readthedocs.io/en/latest/quickstart.html) and
[output documentation](https://msani.readthedocs.io/en/latest/outputs.html)
for details.

## Configuration

MolSanitizer runs with bundled defaults out of the box. To set personal
defaults (for example, tool paths or cluster settings), create a config file:

```bash
msani config init
```

The command prints the location of `msani_configurations.yaml`
(`~/.config/msani/` on Linux, `~/Library/Application Support/msani/` on macOS,
`%LOCALAPPDATA%\msani\` on Windows). Add only the settings you want to override:

```yaml
CORINA: /opt/corina/corina
SLURM_ACCOUNT: my-project
NUMCONFS: 1000
MAX_JOBS: 100
```

Useful commands:

```bash
msani config show                            # merged settings and where each comes from
msani config use /shared/mygroup/msani.yaml  # select an existing config file
msani config reset-path                      # return to the standard location
```

Command-line arguments always take priority over the config file. For a
one-off override, set `MSANI_CONFIG=/path/to/config.yaml`. See the
[configuration guide](https://msani.readthedocs.io/en/latest/config.html) for
precedence rules and batch jobs.

## Documentation

Full documentation, including the methodology behind MolSanitizer, is available
at [msani.readthedocs.io](https://msani.readthedocs.io).


## Contributing

Contributions are welcome — bug reports, feature suggestions, and pull
requests alike. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before getting
started. Python workflows live in `msani/`; native code and bindings are in
[`msani/cpp/`](msani/cpp/README.md).

MolSanitizer is rule-based and draws on experience from prior drug discovery
projects. Suggestions for new filter, tautomer, or protonation rules are
especially appreciated — please open an issue.

## Citation

If you use MolSanitizer in your research, please cite:

> Lam, T.-P.; Pach, S.; Ullmann, P.; et al. MolSanitizer: An open-source
> pipeline to prepare large-scale small-molecule databases for virtual
> screening. *ChemRxiv* **2026**. DOI: [TODO](TODO)

## License

Distributed under the Apache License 2.0. See [LICENSE](LICENSE).

## Contact

- Thua-Phong Lam — phong.lam@icm.uu.se
- Szymon Pach — szymon.pach@icm.uu.se
- Israel Cabeza de Vaca Lopez — israel.cabezadevaca@icm.uu.se
