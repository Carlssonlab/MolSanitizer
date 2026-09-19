import io
import subprocess

from rdkit import Chem


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
