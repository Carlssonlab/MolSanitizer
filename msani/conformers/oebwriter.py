"""OpenEye .oeb.gz output: RDKit conformers -> one multi-conformer OEMol.

Why this exists
---------------
MolSanitizer already covers db2 (UCSF DOCK), pdbqt (AutoDock/Vina) and the
generic sdf/mol2 formats. The one docking format still missing is OpenEye's
binary .oeb.gz, which is what FRED and HYBRID read fastest. Adding it means a
single msani pass can feed DOCK, AutoDock/Vina and OpenEye docking from the
same sanitized library, with no per-user conversion step in between.

What it does
------------
FRED and HYBRID expect ONE OEMol per ligand that carries all of its conformers;
if conformers arrive as separate records (as they do in an sdf), each is read as
a separate molecule and the ligand identity is lost unless the caller regroups
the records by title afterwards. So instead of serialising conformers to text,
this writer transfers the molecular graph once through a single Mol block and
then copies the coordinates of every RDKit conformer straight into the OEMol.
That keeps ligand identity intact, preserves the coordinates exactly rather
than rounding them through a text format, and is about twice as fast as an
SDF round-trip.

Output layout
-------------
``--format oeb``      one ``oeb/<molecule>.oeb.gz`` per molecule
``--format oeb.lib``  one ``oeb/<input>.oeb.gz`` library per input file, which is
                      what FRED and HYBRID take as their database (the same
                      relation db2.tgz has to db2)

Either way the ring-conformer groups of a molecule are merged into a single
OEMol, so every ligand is one record and docks to one result.

Licensing
---------
The OpenEye toolkits are commercial, so they are an optional dependency here
(``pip install MolSanitizer[oe]``) and need the user's own licence file:

    export OE_LICENSE=/path/to/oe_license.txt

The toolkits read that variable themselves; msani only checks that a licence is
present and reports it clearly when it is not.

Authors: Hocine El Khaoudi Enyoury, Phong Lam
"""
from rdkit import Chem

try:
    from openeye import oechem
    OE_AVAILABLE = True
except ImportError:
    OE_AVAILABLE = False
    oechem = None


def require_openeye():
    """Raise unless the OpenEye toolkits are both installed and licensed."""
    if not OE_AVAILABLE:
        raise ImportError(
            'The OpenEye toolkits are required for OEB output. Install the optional '
            'dependency with "pip install MolSanitizer[oe]".'
        )
    if not oechem.OEChemIsLicensed():
        raise RuntimeError(
            'The OpenEye toolkits are installed but not licensed. Point OE_LICENSE at '
            'your licence file, e.g. export OE_LICENSE=/path/to/oe_license.txt'
        )


def rdkit_to_oemol(rdkit_mol, name=None, sddata=None):
    """
    RDKit Mol with N conformers -> one multi-conformer OpenEye OEMol.

    No multi-conformer SDF is created, even in memory.
    A single Mol block is used only to transfer the molecular graph.
    """

    rd_confs = list(rdkit_mol.GetConformers())
    if not rd_confs:
        return None

    # ----------------------------------------------------------
    # 1. Transfer graph/topology using ONE conformer
    # ----------------------------------------------------------
    first_cid = rd_confs[0].GetId()

    molblock = Chem.MolToMolBlock(
        rdkit_mol,
        confId=first_cid
    )

    ims = oechem.oemolistream()
    ims.SetFormat(oechem.OEFormat_MDL)

    if not ims.openstring(molblock):
        raise RuntimeError("Could not open RDKit Mol block")

    graph = oechem.OEGraphMol()

    if not oechem.OEReadMolecule(ims, graph):
        raise RuntimeError("Could not parse RDKit Mol block with OEChem")

    ims.close()

    # ----------------------------------------------------------
    # 2. Promote graph to multiconformer OEMol
    # ----------------------------------------------------------
    oe = oechem.OEMol(graph)

    # OEMol initially contains the conformer copied from graph.
    # Remove it because we'll recreate all conformers from RDKit.
    oe.DeleteConfs()

    oe_atoms = list(oe.GetAtoms())

    if len(oe_atoms) != rdkit_mol.GetNumAtoms():
        raise RuntimeError(
            f"Atom count mismatch: "
            f"RDKit={rdkit_mol.GetNumAtoms()}, "
            f"OpenEye={len(oe_atoms)}"
        )

    # Optional but useful sanity check that atom order survived
    for i, oe_atom in enumerate(oe_atoms):
        rd_atom = rdkit_mol.GetAtomWithIdx(i)

        if rd_atom.GetAtomicNum() != oe_atom.GetAtomicNum():
            raise RuntimeError(
                f"Atom-order mismatch at index {i}: "
                f"RDKit Z={rd_atom.GetAtomicNum()}, "
                f"OpenEye Z={oe_atom.GetAtomicNum()}"
            )

    # ----------------------------------------------------------
    # 3. Copy every RDKit conformer directly
    # ----------------------------------------------------------
    # OEChem atom indices need not be contiguous, so size the array by
    # GetMaxAtomIdx() and address it by the atom's own GetIdx(); rd_idx indexes
    # the RDKit side, where indices are always 0..N-1. The check above ties the
    # two orders together.
    ncoords = 3 * oe.GetMaxAtomIdx()

    for rd_conf in rd_confs:

        coords = oechem.OEFloatArray(ncoords)

        for rd_idx, oe_atom in enumerate(oe_atoms):
            pos = rd_conf.GetAtomPosition(rd_idx)

            idx = 3 * oe_atom.GetIdx()

            coords[idx]     = pos.x
            coords[idx + 1] = pos.y
            coords[idx + 2] = pos.z

        oe.NewConf(coords)

    # ----------------------------------------------------------
    # 4. Metadata
    # ----------------------------------------------------------
    if name is None and rdkit_mol.HasProp("_Name"):
        name = rdkit_mol.GetProp("_Name")

    if name is not None:
        oe.SetTitle(name)

    for key, value in (sddata or {}).items():
        oechem.OESetSDData(oe, str(key), str(value))

    return oe


def append_conformers(oe, rdkit_mol):
    """Append the conformers of another RDKit mol of the SAME molecule to `oe`.

    Used to merge ring-conformer groups: they are copies of one molecule, so the
    atom order is identical and only the coordinates differ. Returns False, and
    leaves `oe` untouched, if the graphs do not line up.
    """
    oe_atoms = list(oe.GetAtoms())
    if len(oe_atoms) != rdkit_mol.GetNumAtoms():
        return False
    for i, oe_atom in enumerate(oe_atoms):
        if rdkit_mol.GetAtomWithIdx(i).GetAtomicNum() != oe_atom.GetAtomicNum():
            return False

    ncoords = 3 * oe.GetMaxAtomIdx()
    for rd_conf in rdkit_mol.GetConformers():
        coords = oechem.OEFloatArray(ncoords)
        for rd_idx, oe_atom in enumerate(oe_atoms):
            pos = rd_conf.GetAtomPosition(rd_idx)
            idx = 3 * oe_atom.GetIdx()
            coords[idx]     = pos.x
            coords[idx + 1] = pos.y
            coords[idx + 2] = pos.z
        oe.NewConf(coords)
    return True


def merge_to_oemols(rdkit_mols, name=None, sddata=None):
    """Ring-conformer groups of one molecule -> a single multi-conformer OEMol.

    FRED and HYBRID dock every record separately, so a ligand split over several
    records would come back as several results. Merging the groups gives one
    record, and one result, per ligand. A group whose graph does not match the
    first one is kept as its own molecule rather than dropped, hence the list.
    """
    merged = []
    for rdkit_mol in rdkit_mols:
        if rdkit_mol.GetNumConformers() == 0:
            continue
        if merged and append_conformers(merged[0], rdkit_mol):
            continue
        oe = rdkit_to_oemol(rdkit_mol, name=name, sddata=sddata)
        if oe is not None:
            merged.append(oe)
    return merged


def open_oeb(filename):
    """Open an .oeb/.oeb.gz output stream. Opening truncates: open a library ONCE."""
    require_openeye()
    ofs = oechem.oemolostream()
    if not ofs.open(filename):
        raise RuntimeError(f"Could not open {filename} for writing")
    return ofs


def write_to_stream(ofs, rdkit_mols, name=None, sddata=None):
    """Write one molecule (all its ring-conformer groups) to an open stream."""
    written = 0
    for oe in merge_to_oemols(rdkit_mols, name=name, sddata=sddata):
        if oechem.OEWriteMolecule(ofs, oe) != oechem.OEWriteMolReturnCode_Success:
            raise RuntimeError(f"Could not write {name} to the OEB stream")
        written += 1
    return written


def write_oeb(rdkit_mols, filename, name=None, sddata=None):
    """Write one molecule to its own .oeb.gz file."""
    ofs = open_oeb(filename)
    try:
        return write_to_stream(ofs, rdkit_mols, name=name, sddata=sddata)
    finally:
        ofs.close()
