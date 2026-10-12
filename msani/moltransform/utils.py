import io
import subprocess

from rdkit import Chem


def run_reaction_preserving_fragments(rxn, mol):
    """Run a single-reactant rule and reattach untouched input components.

    RDKit reaction products omit disconnected components that do not react.
    Use atom provenance to identify those spectators, preserving duplicates,
    stereochemistry, isotopes, and atom properties. Callers sanitize products
    after recombination, just as they sanitize ordinary reaction products.
    """
    outcomes = rxn.RunReactants((mol,))
    if not outcomes:
        return outcomes
    fragment_indices = Chem.GetMolFrags(mol)
    if len(fragment_indices) < 2:
        return outcomes

    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    preserved = []
    for outcome in outcomes:
        represented = {
            atom.GetIntProp('react_atom_idx')
            for product in outcome for atom in product.GetAtoms()
            if atom.HasProp('react_atom_idx')
        }
        product = outcome[0]
        for indices, fragment in zip(fragment_indices, fragments):
            if represented.isdisjoint(indices):
                product = Chem.CombineMols(product, fragment)
        preserved.append((product, *outcome[1:]))
    return tuple(preserved)


def enumerate_stereoisomers_corina(
    smiles: str,
    max_isomers: int,
    onlyUnassigned: bool,
    corina_path: str,
    timeout: float | None = None,
) -> list[str]:
    """Enumerate stereoisomers with CORINA and return isomeric SMILES."""
    if not corina_path:
        raise ValueError("A CORINA executable path is required for CORINA enumeration.")

    stergen_options = "stergen"
    if max_isomers > 0:
        stergen_options += f",msi={max_isomers}"

    if onlyUnassigned:
        stergen_options += ",preserve"

    command = [
        str(corina_path),
        "-i", "t=smiles,scn=1,ncn=2",
        "-o", "t=sdf",
        "-d", stergen_options,
    ]
    try:
        result = subprocess.run(
            command,
            input=f"{smiles} molecule\n".encode("utf-8"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except FileNotFoundError as error:
        raise RuntimeError(f"CORINA executable was not found: {corina_path}") from error
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"CORINA stereoisomer enumeration timed out after {timeout} seconds."
        ) from error

    if result.returncode != 0:
        error_message = result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(
            f"CORINA stereoisomer enumeration failed with exit code "
            f"{result.returncode}: {error_message or 'no error message'}"
        )

    # CORINA may include diagnostic comment lines in its SDF stream.
    output = result.stdout.decode("utf-8", errors="replace")
    cleaned_output = "\n".join(
        line for line in output.splitlines()
        if not line.lstrip().startswith("#")
    )
    supplier = Chem.ForwardSDMolSupplier(
        io.BytesIO(cleaned_output.encode("utf-8")),
        removeHs=True,
        sanitize=True,
    )

    stereoisomers = []
    seen = set()
    for mol in supplier:
        if mol is None:
            continue
        isomer = Chem.MolToSmiles(mol, isomericSmiles=True)
        if isomer not in seen:
            seen.add(isomer)
            stereoisomers.append(isomer)

    return stereoisomers
