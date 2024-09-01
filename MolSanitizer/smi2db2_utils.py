from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms, rdForceFieldHelpers, rdDistGeom, rdMolAlign
from rdkit.Chem.rdchem import Mol, Conformer
import os, glob, shutil

def chair_like(mol, conf, matches_6_member_aliphatic):
    """Check if the conformers contain a non-chair conformation."""
    for ring_atoms in matches_6_member_aliphatic:
        dihedrals = [
            (ring_atoms[i], ring_atoms[(i + 1) % 6], ring_atoms[(i + 2) % 6], ring_atoms[(i + 3) % 6])
            for i in range(6)
        ]
        
        for dihedral in dihedrals:
            angle = rdMolTransforms.GetDihedralDeg(mol.GetConformer(conf), *dihedral)
            # Chair conformation typically has dihedrals close to ±60°
            if not (45 <= abs(angle) <= 75):
                return False
    return True

def find_6_member_aliphatic(mol_H: Mol):
    """Find all 6-membered aliphatic rings in the molecule."""
    return mol_H.GetSubstructMatches(Chem.MolFromSmarts("A1AAAAA1"))

def find_sulfonamide_like_scaffolds(mol_H: Mol):
    """Find all Sulfonamide-like scaffolds (S(O2)-N(R1)R2 or (S(O)(N)-N(R1)(R2))."""
    return mol_H.GetSubstructMatches(Chem.MolFromSmarts("[*:1][S;$(S(=*)=*):2]-!@[N&+0;!$([NH2]):3](-[*,#1:4])-[*,#1:5]"))


def find_non_planar_rings(mol):
    """Find non-planar rings in a molecule."""
    ssr = [set(ring) for ring in Chem.GetSymmSSSR(mol)]
    fused_rings = []

    while ssr:
        fused_set = ssr.pop(0)
        fused = True

        while fused:
            fused = False
            for other_ring in ssr:
                if fused_set.intersection(other_ring):
                    fused_set.update(other_ring)
                    ssr.remove(other_ring)
                    fused = True

        fused_rings.append(fused_set)
    # Exclude aromatic rings
    non_aromatic_rings = [
        ring for ring in fused_rings if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
    ]
    return non_aromatic_rings

def calculate_dihedrals_for_rings(conf: Conformer, ring_atoms):
    """Calculate dihedral angles for all torsions involving four consecutive atoms in the ring."""
    ring_atoms = list(ring_atoms)
    dihedrals = [
        (idx1, idx2, idx3, idx4, rdMolTransforms.GetDihedralDeg(conf, idx1, idx2, idx3, idx4))
        for idx1, idx2, idx3, idx4 in [
            (ring_atoms[i], ring_atoms[(i + 1) % len(ring_atoms)], ring_atoms[(i + 2) % len(ring_atoms)], ring_atoms[(i + 3) % len(ring_atoms)])
            for i in range(len(ring_atoms))
        ]
    ]
    return dihedrals

def normalize_angle(angle):
    """Normalize the angle to the range -180 to 180 degrees."""
    return (angle + 180) % 360 - 180


def are_dihedrals_similar(dihedrals1, dihedrals2, tolerance=10):
    """Check if two sets of dihedrals are within a specified tolerance."""
    for d1, d2 in zip(dihedrals1, dihedrals2):
        # Normalize both angles to -180 to 180 range
        d1_norm = normalize_angle(d1[4])
        d2_norm = normalize_angle(d2[4])
        # Calculate the difference
        difference = abs(d1_norm - d2_norm)
        
        # Consider the circular nature
        if difference > 180:
            difference = 360 - difference
        
        if difference > tolerance:
            return False
            
    return True

def ring_conf_clusters(current_conf, remaining_confs, non_planar_rings, tolerance=10):
    """Cluster conformers based on the dihedral angles of the fused ring systems."""
    clusters = []
    ref_dihedrals = []
    for ring in non_planar_rings:
        dihedrals = calculate_dihedrals_for_rings(current_conf, ring)
        ref_dihedrals.extend(dihedrals)
    #print(f"Ref: {ref_dihedrals}")
    
    for energy, conf, _ in remaining_confs:
        conf_dihedrals = []
        for ring in non_planar_rings:
            dihedrals = calculate_dihedrals_for_rings(conf, ring)
            conf_dihedrals.extend(dihedrals)
        
        if not(are_dihedrals_similar(ref_dihedrals, conf_dihedrals, tolerance)):
            clusters.append((energy, conf, _))
            #print(f"Different cluster: {conf_dihedrals}")
    return clusters

def move_and_rename_mol2_files(name: str, num_rigid_scaffolds: int, VERBOSE = False):
    # Loop over the range of scaffold indices
    for i in range(num_rigid_scaffolds):
        # Define the source file path for the current scaffold index
        src_file = f"3d/{name}/{name}_mol{i}_rotated.mol2"
        
        # Define the new file path and destination
        dest_file = src_file.replace('_rotated.mol2', '.mol2').replace('3d', 'db2')
        
        # Ensure the destination directory exists
        os.makedirs(os.path.dirname(dest_file), exist_ok=True)
        
        # Move and rename the file
        shutil.move(src_file, dest_file)
        if VERBOSE: print(f"Moved and renamed: {src_file} -> {dest_file}")
