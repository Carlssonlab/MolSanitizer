import subprocess
import os
import shutil
import random
import logging

import numpy as np
import yaml

from pandas import DataFrame, concat, read_csv  # only what you use
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolAlign
from rdkit.Chem.rdchem import Mol, Conformer
from scipy.spatial.distance import pdist, squareform
from pathlib import Path

# from msani.filtering import strain_filter
from msani.db2 import mol2db2, mol2


logger = logging.getLogger('msani')

# Define SMARTS patterns for various functional groups
sulfonamide_like_substructure = Chem.MolFromSmarts("[*:1][S;$(S(=*)=*):2]-!@[N&+0;!$([NH2]):3](-[*,#1;!$(C=A):4])-[*,#1;!$(C=A):5]")
substituted_C_cyclohexane = Chem.MolFromSmarts('[!#1]-!@[CH]1-[*]~[*]~[*]~[*]-[A]-1')
flippable_Ns_1 = Chem.MolFromSmarts("[!#1:1]-!@[NH+;!$(N-*=*):2]1-[A:3]-[A:4]-[A]-[A:6]-[A:5]-1")
flippable_Ns_2 = Chem.MolFromSmarts("[*:1]-!@[N+0;!$(N-*=*):2]1-[A:3]-[A:4]-[A]-[A:6]-[A:5]-1")
# substituted_C_cyclohexane = Chem.MolFromSmarts('[!#1:1]-!@[CH:2]1-[A^3:3]-[A^3:4]-[A]-[A^3:6]-[A^3:5]-1')
# substituted_C_cyclohex_23_enyl = Chem.MolFromSmarts('[!#1:1]-!@[CH:2]1-[*^2]~[*^2]-[A^3]-[A^3]-[A^3]-1')
# substituted_C_cyclohex_34_enyl = Chem.MolFromSmarts('[!#1:1]-!@[CH:2]1-[A^3]-[*^2]~[*^2]-[A^3]-[A^3]-1')

conjugated_substituted_nitrogen_5aro = Chem.MolFromSmarts('*-[nX3&+0:1]1[a:2][a:3][a:4][a:5]1')
conjugated_substituted_nitrogen_6aro = Chem.MolFromSmarts('*-[nX3&+0:1]1[a:2][a:3][a:4][a:5][a:6]1')

aro_5_patt = Chem.MolFromSmarts('*-[a:1]1[a:2][a:3][a:4][a:5]1')  # 5 aromatic atoms
aro_6_patt = Chem.MolFromSmarts('*-[a:1]1[a:2][a:3][a:4][a:5][a:6]1')  # 6 aromatic atoms
cycloheptatriene_smarts = Chem.MolFromSmarts('[*^2]1~[*^2]-[*^2]~[*^2]-[*^2]~[*^2]-[A;$([A^3]),$([N^2]),$(A=!@[*!X1])]-1')
cyclohepta_1_4_diene_3_sp2_smarts = Chem.MolFromSmarts('[*^2]1~[*^2]-[A^3,O]-[A^3,O]-[*^2]~[*^2]-[C^2,O,NH0]-1')

barbiturate = Chem.MolFromSmarts('[C;$(C~[OX1,SX1]):1]1~[N:2]~[C;$(C~[OX1,SX1]):3]~[*^2:4]~[*^2:5]~[*:6]~1') # To 0 iteratively four consecutive atoms
hydantoin = Chem.MolFromSmarts('[C;$(C~[OX1,SX1]):1]1~[N:2]~[C;$(C~[OX1,SX1]):3]~[*^2:4]~[A:5]~1') # To 0 iteratively four consecutive atoms
substituted_N_barbi_hydan_like = Chem.MolFromSmarts('*~[C^2,N^2:1][C^2,N^2:2][C^2,N^2:3]')
amide_substructure = Chem.MolFromSmarts('[$(C=O):1]!@[NX3&+0:2]') #primary, secondary amide for constrained planarity only

const_rule = [(-120, 30, 30, 1), (-60, 30, 30, 1), (0, 30, 30, 1), (60, 30, 30, 1), (120, 30, 30, 1), (180, 30, 30, 1)]


rigid_rule_files = Path(__file__).parent.parent / 'Data' / 'rigid_part_rules.txt'
rigid_rules = read_csv(rigid_rule_files, header=None, sep =r'\s+', names=['SMARTS','label'])
rigid_rules['mol'] = rigid_rules['SMARTS'].apply(lambda x: Chem.MolFromSmarts(x))
rotatable_pattern=r'*~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~*'

with open(Path(__file__).parent.parent / 'msani_configurations.yaml') as confFile:
    msani_configurations = yaml.full_load(confFile)
CORINA_EXE = msani_configurations['CORINA']

# These below are for the new more deterministic method
symmetric_patterns_file = Path(__file__).parent.parent / 'Data' / 'symmetric_smarts.txt'
symmetric_patterns_df = read_csv(symmetric_patterns_file, sep=r'\s+', header=None, names=['pattern', 'name', 'num_scaled'])
symmetric_patterns_df['mol'] = symmetric_patterns_df['pattern'].apply(lambda x: Chem.MolFromSmarts(x))
prim_amidines_guanidines_pattern_mol = [Chem.MolFromSmarts('[#1:1][NH2,NX3H1:2]!@-[CX3+0;$(C(~[NH2])(~[NH2])~*):3]~[NH2:4]'), 
                                        Chem.MolFromSmarts('[#1:1][NX3H2:2]!@-[#6:2]~[#7&+1]')] #in ring

def embed_smiles_corina(smiles, name, numringconfs, VERBOSE):
    '''
    Embed the SMILES string using CORINA and return the mol, net_charge,
    rigid_scaffolds, and flexible_scaffolds
    '''
    # Input SMILES string
    smiles_input = f"{smiles} {name}"

    # Command to run CORINA
    # rc: multiple ring confs; flapn: Flap ring nitrogen atoms,
    # de=6: energy window, mc=nconfs: write nring conf, wh: write hydrogens, sanpyr: make sulfonamide pyramidal
    command = [
        CORINA_EXE,
        "-i", "t=smiles,scn=1,ncn=2",
        "-o", "t=mol2",
        "-d", f"rc,flapn,de=6,mc={numringconfs},wh,sanpyr"
    ]

    # Run the command
    result = subprocess.run(
        command,
        input=smiles_input.encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    # Check for errors
    if result.returncode != 0:
        print("Error:", result.stderr.decode())
    else:
        # Decode and process output
        output = result.stdout.decode()

        # Remove comment lines (starting with #)
        lines = [line for line in output.splitlines() if not line.strip().startswith("#")]

        # Rejoin cleaned lines
        cleaned_output = "\n".join(lines)

        # Split by the MOL2 section header (and keep it in each block)
        raw_blocks = cleaned_output.split("@<TRIPOS>MOLECULE")

        # Add back the header to each block (except the first if it's empty)
        mol2_blocks = [
            "@<TRIPOS>MOLECULE\n" + block.strip() + "\n" 
            for block in raw_blocks if block.strip()
        ]
        if VERBOSE:
            print(f"\tNumber of ring conformers: {len(mol2_blocks)}")
        ring_confs = []
        # Now you have a list of strings, each containing one MOL2 molecule
        for i, mol in enumerate(mol2_blocks, 1):
            rdkit_mol = Chem.MolFromMol2Block(mol, removeHs=False, sanitize=True)
            if rdkit_mol:
                ring_confs.append(rdkit_mol)
        mol2_string = mol2_blocks[0] if len(mol2_blocks) > 0 else None
        return mol2_string, ring_confs
    

def find_flipped_nitrogen(mol_H: Mol):
    '''
    Find the flippable nitrogen in the molecule. Mainly for reinforcing the equatorial of the *-N of piperidine and piperazine
    '''
    return mol_H.GetSubstructMatches(flippable_Ns_1) + mol_H.GetSubstructMatches(flippable_Ns_2)

def find_flipped_carbon(mol_H: Mol):
    '''
    Find the flippable carbon in the molecule. Mainly for substituted cyclohexane
    '''
    return mol_H.GetSubstructMatches(substituted_C_cyclohexane) #+ \
        # mol_H.GetSubstructMatches(substituted_C_cyclohex_23_enyl) +\
            # mol_H.GetSubstructMatches(substituted_C_cyclohex_34_enyl)

def find_conjugated_substituted_nitrogen_5aro(mol_H: Mol):
    '''
    Find the conjugated substituted nitrogen in the molecule with 5 aromatic atoms
    Format: *-[nX3&+0:1]1[a:2][a:3][a:4][a:5]1

    '''
    return mol_H.GetSubstructMatches(conjugated_substituted_nitrogen_5aro)


def find_conjugated_substituted_nitrogen_6aro(mol_H: Mol):
    '''
    Find the conjugated substituted nitrogen in the molecule with 6 aromatic atoms
    Format: *-[nX3&+0:1]1[a:2][a:3][a:4][a:5][a:6]1

    '''
    return mol_H.GetSubstructMatches(conjugated_substituted_nitrogen_6aro)

def find_barbiturates(mol_H: Mol):
    matches = mol_H.GetSubstructMatches(barbiturate)
    return matches


def find_hydantoins(mol_H: Mol):
    matches = mol_H.GetSubstructMatches(hydantoin)
    return matches

def find_substituted_N_barbi_hydan_like(mol_H: Mol, barbiturate: tuple, hydantoin: tuple):
    matches = mol_H.GetSubstructMatches(substituted_N_barbi_hydan_like)
    filtered_matches = []
    for match in matches:
        for barbi in barbiturate:
            if set(match[1:]).issubset(barbi) and match[0] not in barbi:   
                filtered_matches.append(match)
        for hydan in hydantoin:
            if set(match[1:]).issubset(hydan) and match[0] not in hydan:   
                filtered_matches.append(match)
    return filtered_matches
        

def find_amide(mol_H: Mol):
    '''
        Find the amide in the molecule. O=C-N(-R)-H
    '''
    return mol_H.GetSubstructMatches(amide_substructure)

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


def ring_conf_clusters(current_conf_ring_descriptors, remaining_confs):
    clusters = []
    for remaining_conf in remaining_confs:
        if not current_conf_ring_descriptors == remaining_conf[2:-1]:
            clusters.append(remaining_conf)
    return clusters

def move_and_rename_mol2_files(name: str, num_rigid_scaffolds: int, numPossibleRigidAlignments: int, VERBOSE = False):
    # Loop over the range of scaffold indices
    for i in range(num_rigid_scaffolds):
        for j in range(numPossibleRigidAlignments):
            # Define the source file path for the current scaffold index
            src_file = f"3d/{name}/{name}_mol{i}_align{j}.mol2"
            
            # Define the new file path and destination
            dest_file = src_file.replace('3d', 'db2')
            
            # Ensure the destination directory exists
            os.makedirs(os.path.dirname(dest_file), exist_ok=True)
            
            # Move and rename the file
            shutil.move(src_file, dest_file)
            if VERBOSE: print(f"\tMoved and renamed: {src_file} -> {dest_file}")



def get_flexible_ring(mol: Mol, srlib):
    """Extract the flexible rings from the molecule based on the given rules 
    srlib. Only consider the rings from SSR to exclude fused rings.

    Args:
        mol (Mol): Rdkit molecule object
        srlib (SmallRingLibrary): SmallRingLibrary object containing the rules for flexible rings

    Returns:
        planar_rings (set): Set containing the indices of the planar rings
        non_planar_rings (list): List containing the indices of the non-planar rings
    """
    ssr = [set(ring) for ring in Chem.GetSymmSSSR(mol)]
    planar_rings = set()
    sorted_planar_rings = set()
    for rule in srlib.planar:
        matches = mol.GetSubstructMatches(rule[2])
        if len(matches) > 0:
            for match in matches:
                sorted_match = tuple(sorted(match))
                if sorted_match not in sorted_planar_rings:
                    sorted_planar_rings.add(sorted_match)
                    planar_rings.add((match))
    # Pre-sort the planar_rings and ssr sets
    sorted_ssr = {tuple(sorted(t)) for t in ssr}         
    non_planar_rings = []
    set_non_planar_rings= set()

    for rule in srlib.non_planar:
        matches = mol.GetSubstructMatches(rule[2])
        if len(matches) > 0:
            for match in matches:
                sorted_match = tuple(sorted(match))
                if sorted_match in sorted_planar_rings: continue 
                if sorted_match not in sorted_ssr: continue
                if sorted_match in set_non_planar_rings: continue
                set_non_planar_rings.add(sorted_match)
                non_planar_rings.append((rule[0], match, [reference_set for reference_set in rule[3]]))
    return planar_rings, non_planar_rings

def calculate_dihedrals_for_rings(conf: Conformer, ring_atoms):
    """Calculate dihedral angles for all torsions involving four consecutive atoms in the ring."""
    ring_atoms = list(ring_atoms)
    dihedrals = []
    for i in range(len(ring_atoms)):
        dihedrals.append(rdMolTransforms.GetDihedralDeg(conf, ring_atoms[i], ring_atoms[(i + 1) % len(ring_atoms)], ring_atoms[(i + 2) % len(ring_atoms)], ring_atoms[(i + 3) % len(ring_atoms)]))
    return dihedrals

def dfs(mol, current_idx, visited, path, forbidden_idxs):
    visited.add(current_idx)
    path.append(current_idx)
    
    for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        # Check if the atom is not visited, not forbidden, and not a hydrogen atom
        if neighbor_idx not in visited and neighbor_idx not in forbidden_idxs:
            if mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() != 1:  # Avoid hydrogen atoms
                dfs(mol, neighbor_idx, visited, path, forbidden_idxs)

def get_dfs_path(mol, start_idx, forbidden_idxs):
    visited = set()
    path = []
    dfs(mol, start_idx, visited, path, forbidden_idxs)
    return path

def identical_substituents(mol, idx2, idx3, idx4, idx5):
    '''
    A symmetrical check using the fragment SMILES and depth-first search
    '''
    # Create forbidden indices set
    forbidden_idxs = {idx3, idx2}

    # Get paths using DFS
    path4 = get_dfs_path(mol, idx4, forbidden_idxs)
    path5 = get_dfs_path(mol, idx5, forbidden_idxs)
    # Generate SMILES for the paths
    smiles4 = Chem.MolFragmentToSmiles(mol, atomsToUse=path4, rootedAtAtom=idx4)
    smiles5 = Chem.MolFragmentToSmiles(mol, atomsToUse=path5, rootedAtAtom=idx5)

    # Compare the SMILES strings
    return smiles4 == smiles5

def is_equatorial(conf, atom_idx):
    """
    Determines if a substituent is in an equatorial position on a cyclohexane ring.

    This function calculates two dihedral angles around a specific atom in the molecule
    and checks if they both fall within the range that indicates an equatorial position 
    (approximately 150-180 degrees in absolute value).

    Parameters
    ----------
    conf : rdkit.Chem.rdchem.Conformer
    atom_idx : list or tuple
        List of atom indices defining the relevant atoms for dihedral calculations.
        Requires at least 7 indices:
        - atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3]: First dihedral angle
        - atom_idx[0], atom_idx[1], atom_idx[6], atom_idx[5]: Second dihedral angle

    Returns
    -------
    bool
        True if the substituent is equatorial (both dihedrals between 150-180 degrees),
        False otherwise
    """
    
    dihedral1 = rdMolTransforms.GetDihedralDeg(conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3])
    dihedral2 = rdMolTransforms.GetDihedralDeg(conf, atom_idx[0], atom_idx[1], atom_idx[6], atom_idx[5])
    if 150 <= abs(dihedral1) <= 180 and 150 <= abs(dihedral2) <= 180:
            return True
    return False
    

def find_sulfonamide_like_scaffolds(mol_H: Mol):
    """Find all Sulfonamide-like scaffolds (S(O2)-N(R1)R2 or (S(O)(N)-N(R1)(R2))."""
    preliminary_sulfonamide = mol_H.GetSubstructMatches(sulfonamide_like_substructure)
    matches_sulfonamide = []
    for (a, b, c, d, e) in preliminary_sulfonamide:
        if identical_substituents(mol_H, b, c, d, e): continue
        else: matches_sulfonamide.append((a, b, c, d, e))
    return matches_sulfonamide

def find_cycloheptatriene(mol_H: Mol):
    """Find cycloheptatriene or cyclohepta-1,4-diene-3-sp2 substructure in the molecule."""
    return mol_H.GetSubstructMatches(cycloheptatriene_smarts) + mol_H.GetSubstructMatches(cyclohepta_1_4_diene_3_sp2_smarts)

def classify_confs(conf, 
                   energy, 
                   non_planar_rings, 
                   flippable_Ns, 
                   flippable_Cs, 
                   sulfo_matches, 
                   conf_ring_descriptors_df, 
                   tolerance=25):
    
    temp_dict = {
        'Conformer': conf,
        'Energy': energy
    }
    
    # Process non-planar rings
    for ring_idx, (ring_name, ring_atoms, reference_sets) in enumerate(non_planar_rings):
        dihedral_one_ring = calculate_dihedrals_for_rings(conf, ring_atoms)
        sign_vector = np.sign(dihedral_one_ring)
        
        ring_conf_id = -1
        for idx, reference_set in enumerate(reference_sets):
            if np.any(np.sign(reference_set) * sign_vector < 0):
                continue
            abs_diff = np.abs(np.array(dihedral_one_ring) - np.array(reference_set))
            if np.any(abs_diff > tolerance):
                continue
            ring_conf_id = idx
            break
        
        temp_dict[f"{ring_name}{ring_idx}"] = ring_conf_id

    # Process flippable Nitrogens
    flippable_N_descriptors = sum([1 if is_equatorial(conf, atom_idx) else 0 for atom_idx in flippable_Ns])
    # print(f"Flippable Nitrogens: {flippable_N_descriptors}")
    temp_dict['equatorial_subs_Ns'] = flippable_N_descriptors if flippable_Ns else -1

    # Process substituted cyclohexane
    aliphatic_cyclohexane_descriptors = sum([1 if is_equatorial(conf, atom_idx) else 0 for atom_idx in flippable_Cs])
    # print(f"Aliphatic cyclohexane: {aliphatic_cyclohexane_descriptors}")
    temp_dict['equatorial_subs_Cs'] = aliphatic_cyclohexane_descriptors if flippable_Cs else -1
    
    # Process sulfo matches
    sulfo_descriptors = tuple([1 if rdMolTransforms.GetDihedralDeg(conf, d, b, c, e) > 0 else 0 for (a, b, c, d, e) in sulfo_matches])
    temp_dict['sulfo_descriptors'] = sulfo_descriptors if sulfo_descriptors else [-1]

    # Create dataframe and append to conf_ring_descriptors_df
    temp_df = DataFrame([temp_dict])
    conf_ring_descriptors_df = concat([conf_ring_descriptors_df, temp_df], ignore_index=True)

    return conf_ring_descriptors_df

def remove_unfavorable_confs(conf_ring_descriptors_df: DataFrame, name: str ='0')-> DataFrame:
    for column in conf_ring_descriptors_df.columns[2:-2]:
        if (conf_ring_descriptors_df[column] == -1).all():
            conf_ring_descriptors_df.drop(columns=[column], inplace=True)
            #print(f"\t While handling {name}, found no good conformation of ring: {column}")

    # Step 2: Drop rows where any value in remaining columns (2 onwards) is -1
    for column in conf_ring_descriptors_df.columns[2:-2]:
        conf_ring_descriptors_df = conf_ring_descriptors_df[conf_ring_descriptors_df[column] != -1]
    return conf_ring_descriptors_df


def is_slurm_job():
    '''Check if SLURM_JOB_ID is present in environment variables'''
    return 'SLURM_JOB_ID' in os.environ

def remove_folders(folders_to_remove: list):
    ''' Remove specified folders and all their contents if they exist'''
    for folder in folders_to_remove:
        if os.path.exists(folder) and os.path.isdir(folder):
            try:
                shutil.rmtree(folder)
            except Exception as e:
                print(f"Error removing folder {folder}: {e}")

def canonicalize_if_smiles(query: str):
    """Canonicalize only if it's a valid SMILES; return unchanged if SMARTS or invalid."""
    mol = Chem.MolFromSmiles(query)
    if mol:
        return Chem.CanonSmiles(query)
    else: 
        return query  # Return SMARTS or invalid input unchanged



def precompute_bonded_and_same_parent_pairs(mol):
    """
    Precompute bonded atom pairs and pairs of atoms that share the same parent (common neighbor).

    Parameters:
    mol (rdkit.Chem.Mol): The RDKit molecule object.

    Returns:
    bonded_pairs (set): Set of tuples representing bonded atom pairs.
    same_parent_pairs (set): Set of tuples representing atoms that share the same parent atom.
    """
    bonded_pairs = set()
    same_parent_pairs = set()

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        neighbors = atom.GetNeighbors()
        atom_idx = atom.GetIdx()

        # Get bonded pairs
        for neighbor in neighbors:
            neighbor_idx = neighbor.GetIdx()
            bonded_pairs.add((atom_idx, neighbor_idx))
            bonded_pairs.add((neighbor_idx, atom_idx))  # Symmetric bond

        # Get atoms that share the same parent (common neighbors)
        if len(neighbors) > 1:  # Only meaningful for atoms with more than 1 neighbor
            neighbor_indices = [n.GetIdx() for n in neighbors]
            for i in range(len(neighbor_indices)):
                for j in range(i + 1, len(neighbor_indices)):
                    same_parent_pairs.add((neighbor_indices[i], neighbor_indices[j]))
                    same_parent_pairs.add((neighbor_indices[j], neighbor_indices[i]))  # Symmetric relation

    return bonded_pairs, same_parent_pairs


def check_too_close_nonbonded_atoms(conformer, mol, bonded_pairs, same_parent_pairs, threshold=1.6):
    """
    Check if any non-bonded and non-same-parent atoms in a given conformer are too close to each other.

    Parameters:
    conformer (rdkit.Chem.rdchem.Conformer): The RDKit conformer object.
    mol (rdkit.Chem.Mol): The RDKit molecule object.
    bonded_pairs (set): Precomputed set of bonded atom pairs.
    same_parent_pairs (set): Precomputed set of atoms that share the same parent atom.
    threshold (float): The distance threshold below which atoms are considered too close.

    Returns:
    bool: True if any non-bonded atoms are too close, False otherwise.
    """
    positions = conformer.GetPositions()
    num_atoms = mol.GetNumAtoms()

    # Compute pairwise distances between all atoms
    pairwise_distances = pdist(positions)

    # Convert pairwise distances into a square matrix form
    distance_matrix = squareform(pairwise_distances)

    # Iterate over all atom pairs and check distances
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            # Skip if the atoms are bonded or share a common parent
            if (i, j) in bonded_pairs or (i, j) in same_parent_pairs:
                continue

            # Check if the distance is below the threshold
            if distance_matrix[i, j] < threshold:
                #print(i, j, distance_matrix[i, j])
                return True

    return False



def count_confs_by_rotbonds(mol,
                            rot_bonds,
                            amide_bonds,
                            ignoretorlib=False,
                            torlib=None,
                            VERBOSE=False):
    """
    Count the number of conformations based on rotatable bonds and reorder bonds with terminal
    atoms at the beginning.

    Args:
    mol (rdkit.Chem.Mol): The RDKit molecule object.
    rot_bonds (list): List of rotatable bonds, each represented as a tuple of atom indices.
    amide_bonds (list): List of tuples representing amide linkages.
    ignoretorlib (bool): If True, ignore the torsion library and use a constant rule for amide atoms.
    torlib (TorsionLibrary): Torsion library object to match dihedral angles.
    VERBOSE (bool): If True, print detailed information about the rotatable bonds and matched rules.

    Returns:
    bond_to_rule (dict): A dictionary mapping each rotatable bond to its corresponding rule.
    """

    matched_rules = torlib.get_match_dihedral(mol, mode = 'random')
    amide_atoms = set()
    if (amide_bonds):
        for b, c in amide_bonds:
            amide_atoms.add(b)
            amide_atoms.add(c)
    # Pre-compute a dictionary of rules by bond
    rule_by_bond = {}
    for rule in matched_rules:
        central = tuple(sorted((rule[1][1], rule[1][2])))
        rule_by_bond[central] = rule

    # Then do a single pass through the rotatable bonds
    bond_to_rule = {}
    for bond in rot_bonds:
        if bond in rule_by_bond:
            rule = rule_by_bond[bond]
            rule_copy = list(rule)
            # Only apply ignoretorlib if not amide bonds
            if ignoretorlib and len(set(bond) & amide_atoms) <= 1:
                rule_copy[2] = const_rule
            bond_to_rule[bond] = rule_copy
    for _, sym_row in symmetric_patterns_df.iterrows():
        matches = mol.GetSubstructMatches(sym_row['mol'])
        if not matches:
            continue
        if VERBOSE:
            print(f"\tFound symmetric pattern: {sym_row['name']}")
        period = 360 / sym_row['num_scaled']
        check_symmetric_using_dfs = True if sym_row['name'].startswith('general') else False
        for match in matches:
            if check_symmetric_using_dfs:
                # Use a dfs to traverse in the Molecule graph, then compare the fragment smiles
                idx2, idx3, idx4, idx5 = match[0], match[1], match[2], match[-1]
                if not(identical_substituents(mol, idx2, idx3, idx4, idx5)): continue
            
            bond_key = tuple(sorted(match[:2]))
            if bond_key not in bond_to_rule:
                continue
            if VERBOSE:
                print(f'\tRemove duplicated rotation for {bond_key}')
            rule = bond_to_rule[bond_key]
            temp = []
           
            for peak in rule[2]:
                if period == 1:
                    temp.append((peak[0], peak[1], 30, peak[3]))
                    break
                is_similar = any(period - 10 < abs(peak[0] - existing_peak[0]) < period + 10
                    for existing_peak in temp)
                if not is_similar:
                    temp.append(peak)
            rule[2] = temp
    
    total_confs = 1
    
    match_torlib_clean = []
    for rule in bond_to_rule.values():
        match_torlib_clean.append(rule)
        total_angles = 0
        for peak in rule[2]:
            total_angles += (peak[2]*2/30 + 1)
        total_confs *= total_angles
    if VERBOSE:
        print("\nFinal processed torsion rules:")
        for rule in match_torlib_clean:
            print(f"\t{rule}")

    return int(total_confs), match_torlib_clean

def get_random_angle(mean, tolerance, method='gauss', rounding = False):
    """
    Generate a random angle value from a Gaussian distribution given the expected mean and standard deviation,
    and limit it within the specified range. Normalize the result to the [-180, 180] degree range.

    Args:
    mean (float): The expected mean angle.
    tolerance (float): The tolerance .
    method (str): The method to use for generating the random angle ('gauss' or 'uniform').

    Returns:
    float: A random angle normalized to the [-180, 180] degree range.
    """

    # Generate a random angle within the specified Gaussian distribution and range limits
    if tolerance == 0: return mean
    while True:
        if method == 'gauss':
            random_angle = random.gauss(mean, tolerance)
        elif method == 'uniform':
            random_angle = random.uniform(mean - tolerance, mean + tolerance)
        if rounding: random_angle = round(random_angle/30) * 30
        if mean-tolerance < random_angle < mean + tolerance: 
            break

    # Normalize to the [-180, 180] range
    normalized_angle = (random_angle + 180) % 360 - 180
    return normalized_angle


def find_rigid_part(mol, request_alignment=None):
    '''Find fused ring system of the molecule as the rigid part, 
    if none, use the hierarchical rules in rigid_part_rules.txt'''

    rigid_part = []
    rule_label = None
    
    # If the user request for only a specific ring as rigid segment.
    if request_alignment:
        matches = mol.GetSubstructMatches((request_alignment))
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        if len(matches) != 0:
            if len(matches) > 1: logger.warning(f"Multiple matches found for the requested alignment: {smiles}")
            for match in matches:
                rigid_part.append(match)
            rule_label = None
        else:
            rigid_part = []
            rule_label = None
        return rigid_part, rule_label
    
    # Find fused ring systems
    ssr = [set(ring) for ring in Chem.GetSymmSSSR(mol)]

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
        rigid_part.append(fused_set)
    
    rigid_part = sorted(rigid_part, key = lambda x: len(x), reverse = True)
    if len(rigid_part) == 0: 
        for rule in rigid_rules.itertuples():
            matches = mol.GetSubstructMatches(rule.mol)
            if len(matches) > 0:
                rigid_part = [matches[0]]
                rule_label = rule.label
                break
    return rigid_part, rule_label
    

def log_error(smiles, name):
    with open('msani_error.log', 'a') as f:
        f.write(f"{smiles} \t {name}\n")

def Align_ConvertToDb2(ring_conf, rigid_scaffold, solv_obj, name, smiles, longname):
    """
    Align the all the conformers to the rigid scaffold (ring) and convert it to the DB2 string.
    """
    mol2_obj = mol2.Mol2(mol2fileName=f'{name}.mol2')
    mol2_obj.cleanConfs()
    mol2_obj.longname = longname if longname else "fake"
    mol2_obj.smiles = smiles
    for conf_id in range(ring_conf.GetNumConformers()):
        mol2_obj.atomXyz.append([])  # Initialize a list for atom coordinates
        rdMolAlign.AlignMol(ring_conf, ring_conf, conf_id, 0, atomMap=[(i, i) for i in rigid_scaffold])
        conf = ring_conf.GetConformer(conf_id)
        for atom_idx in range(ring_conf.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            mol2_obj.atomXyz[-1].append((float(pos.x), float(pos.y), float(pos.z)))
            # Set the number of conformations for this Mol2 object
        mol2_obj.xyzCount = ring_conf.GetNumConformers()
    return mol2db2.mol2db2_quick_ver2(mol2_obj, solv_obj)

def is_similar_conformer(new_dihedrals, exist, tol = 30.0):
    if exist.shape[0] == 0:
        return False

    diffs = np.abs(exist - new_dihedrals)
    diffs = np.minimum(diffs, 360 - diffs)
    is_similar = np.all(diffs <= tol, axis=1)
    return np.any(is_similar)


# All deterministic version of torsional sampling will be available here
def bond_centrality(bond, dists):
    # Centrality: shortest average distance from bond atoms to all others

    a1, a2 = bond
    centrality = sum(dists[a1]) + sum(dists[a2])
    return centrality

def softmax_weights(scores, T=1.0):
    """Lower scores → larger weight.  T↑  ⇒ flatter, T↓ ⇒ steeper."""
    scores = np.asarray(scores, dtype=float)
    w = np.exp(scores / T)          
    return w / w.sum()

def get_importance_order(mol, rot_bonds, debug = False):
    """
    Calculate the importance order of rotatable bonds based on their centrality in the molecule.
    """
    dists = Chem.GetDistanceMatrix(mol)
    bond_scores = [bond_centrality((b[0],b[1]), dists) for b in rot_bonds]

    # Get sorting indices - lower score is more central (smaller distance sum)
    sorting_indices = np.argsort(bond_scores)

    # Create weights based on position - most central bonds get higher weights
    position_weights = [0]*len(sorting_indices)
    for idx in range(len(sorting_indices)):
        # Assign weights based on position in sorted order
        # More central bonds (lower indices) get higher weights
        # This is a simple linear weighting scheme
        position_weights[sorting_indices[idx]] = (len(sorting_indices) - idx) / 10 # 10 is to reduce the difference
    
    position_weights = softmax_weights(position_weights)
    if debug:
        print("Bond scores (lower is more central):", bond_scores)
        print("Softmax weights:", position_weights)
    return position_weights

def getDihedralMatches(mol):
    '''return list of atom indices of dihedrals'''
    global rotatable_pattern
    qmol = Chem.MolFromSmarts(rotatable_pattern)
    matches = mol.GetSubstructMatches(qmol)
    #these are all sets of 4 atoms, uniquify by middle two
    uniqmatches = []
    seen = set()
    for (_,b,c,_) in matches:
        if ((b,c) not in seen) and ((c,b) not in seen):
            seen.add((b,c))
            uniqmatches.append(tuple(sorted((b,c))))
    return uniqmatches


def discretinize_dihedrals(typical, tolerance, step = 30):
    '''
    Discretinize a dihedral angle into a list of angles.
    typical: the typical dihedral angle (peak of the distribution)
    tolerance: the tolerance around the typical angle
    step: the step size for discretinization
    '''
    if tolerance < step: return [typical]
    # n_steps = int(tolerance / step)
    # angles = [typical + i * step for i in range(-n_steps, n_steps + 1)]
    angles = [typical, typical - step, typical + step] #Only sample 3 angles for each dihedral

    normalized_angles = [round((angle + 180) % 360 - 180, 1) for angle in angles]
    return normalized_angles

def angular_diff(a, b):
    diff = abs(a - b) % 360
    return min(diff, 360 - diff)

def filter_symmetric_angles(angles, scores, symmetry_angle=180, tolerance=10):
    """
    Filter out angles that are approximately symmetry_angle degrees apart from a previously kept angle.

    Parameters:
        angles (list of float): Angles in degrees.
        scores (list of float): Corresponding scores.
        symmetry_angle (float): The symmetry angle (e.g., 180, 120).
        tolerance (float): Allowed deviation for approximate symmetry.

    Returns:
        Tuple: (filtered_angles, filtered_scores)
    """
    kept_angles = []
    kept_scores = []

    for _, (a, s) in enumerate(zip(angles, scores)):
        # Check if a is similar to any already kept angle
        is_similar = any(
            symmetry_angle - tolerance <= angular_diff(a, kept) <= symmetry_angle + tolerance
            for kept in kept_angles
        )
        if not is_similar:
            kept_angles.append(a)
            kept_scores.append(s)

    return kept_angles, kept_scores

def count_confs_by_rotbonds_v2(mol,
                               rot_bonds,
                               ignoretorlib = False,
                               amide_bonds = None,
                               torlib = None,
                               VERBOSE=False):
    """
    Estimates the number of conformations by analyzing rotatable bonds and torsion rules.
    Adjusts torsions for amide bonds and symmetric patterns.

    Args:
        mol (rdkit.Chem.Mol): The input molecule.
        rot_bonds (list): List of rotatable bonds.
        ignoretorlib (bool): If True, ignore the torsion library and use all possible angles differ by 30 degrees.
        amide_bonds (list): List of amide bonds to zero out fluctuations.
        torlib (torsions.Torsional Library object)
        VERBOSE (bool): If True, print detailed steps.

    Returns:
        tuple: (total_confs, angle_map, score_map)
            - total_confs (int): Estimated number of conformations.
            - angle_map (dict): Mapping of bond indices to possible angles.
            - score_map (dict): Mapping of bond indices to scores for each angle.
    """
    # Step 1: Match torsion rules and rotatable bonds
    rot_bonds = rot_bonds.copy()
    matched_rules = torlib.get_match_dihedral(mol, mode = 'fixed')

    amide_atoms = set()
    if (amide_bonds):
        for b, c in amide_bonds:
            amide_atoms.add(b)
            amide_atoms.add(c)

    # Pre-compute a dictionary of rules by bond
    rule_by_bond = {}
    for rule in matched_rules:
        central = tuple(sorted((rule[1][1], rule[1][2])))
        rule_by_bond[central] = rule

    # Then do a single pass through the rotatable bonds
    bond_to_rule = {}
    for bond in rot_bonds:
        if bond in rule_by_bond:
            rule = rule_by_bond[bond]
            rule_copy = list(rule)
            
            # Only apply ignoretorlib if not amide bonds
            if ignoretorlib and len(set(bond) & amide_atoms) <= 1:
                rule_copy[2] = const_rule
                
            bond_to_rule[bond] = rule_copy

    # Step 2: Handle planar substructures (amides, amidines - zero out fluctuations)
    if (amide_bonds):
        for amide_match in amide_bonds:
            bond_key = tuple(sorted(amide_match[0:2]))
            if bond_key not in bond_to_rule:
                continue
            rule = bond_to_rule[bond_key]
            adjusted_peaks = []
            for peak in rule[2]:
                if isinstance(peak, tuple) and len(peak) >= 3:
                    adjusted_peaks.append((peak[0], 0, 0) + peak[3:])
                else:
                    adjusted_peaks.append(peak)
            rule[2] = adjusted_peaks

    # Zero out fluctuations for primary amidines, guanidines:
    uniq_matches_prim_amidines_guanidines = set()
    for patt in prim_amidines_guanidines_pattern_mol:
        for match in mol.GetSubstructMatches(patt):
            uniq_matches_prim_amidines_guanidines.add(tuple(sorted(match[1:3])))

    for bond_key in list(uniq_matches_prim_amidines_guanidines):
        if bond_key not in bond_to_rule: continue
        rule = bond_to_rule[bond_key]
        adjusted_peaks = []
        for peak in rule[2]:
            if isinstance(peak, tuple) and len(peak) >= 3:
                adjusted_peaks.append((peak[0], 0, 0) + peak[3:])
            else:
                adjusted_peaks.append(peak)
        rule[2] = adjusted_peaks

    # Step 3: Apply symmetry filtering to avoid redundant angles
    for _, sym_row in symmetric_patterns_df.iterrows():
        matches = mol.GetSubstructMatches(sym_row['mol'])
        if not matches:
            continue
        if VERBOSE:
            print(f"\tFound symmetric pattern: {sym_row['name']}")
        period = 360 / sym_row['num_scaled']
        check_symmetric_using_dfs = True if sym_row['name'].startswith('general') else False
        for match in matches:
            if check_symmetric_using_dfs:
                # Use a dfs to traverse in the Molecule graph, then compare the fragment smiles
                idx2, idx3, idx4, idx5 = match[0], match[1], match[2], match[-1]
                if not(identical_substituents(mol, idx2, idx3, idx4, idx5)): continue
            
            bond_key = tuple(sorted(match[:2]))
            if bond_key not in bond_to_rule:
                continue
            if period == 1: 
                bond_to_rule.pop(bond_key)
                rot_bonds.remove(bond_key)
                continue

            if sym_row['name'] == 'tri_large_halogeno_methyl':
                # Special case for tri_large_halogeno_methyl, where we want to permute + 30o.
                rule = bond_to_rule[bond_key]
                peaks = rule[2]
                angles = [normalize_angle(peaks[0][0] + x) for x in [0, -30, 30]]
                scores = [peaks[0][3], round(peaks[0][3]/2, 2), round(peaks[0][3]/2, 2)]
                rule[2] = list(zip(angles, scores))
            if VERBOSE:
                print(f'\tRemove duplicated rotation for {bond_key}')
            rule = bond_to_rule[bond_key]
            angles, scores = [], []
            if len(rule[2][0]) < 4:
                # Already reduced for this angle, skip further reduced
                continue
            for peak in rule[2]:
                angle_list = discretinize_dihedrals(peak[0], peak[2])
                angles.extend(angle_list)
                scores.extend([peak[3]] * len(angle_list))
            filtered_angles, filtered_scores = filter_symmetric_angles(angles, scores, period, 10)
            rule[2] = list(zip(filtered_angles, filtered_scores))


    # Step 4: Discretize angles and estimate total possible conformations
    angle_map = {}
    score_map = {}
    total_confs = 1
    for bond_idx, rule in enumerate(bond_to_rule.values()):
        name, atom_indices, peak_list = rule
        angle_list = []
        score_list = []
        total_angles = 0
        for peak in peak_list:
            if len(peak) < 4:
                # fallback case (e.g., symmetric angle filtering already applied)
                angle_list = [peak[0] for peak in peak_list]
                score_list = [peak[1] for peak in peak_list]
                total_angles = len(angle_list)
                break
            angle_vals = discretinize_dihedrals(peak[0], peak[2])
            angle_list.extend(angle_vals)
            score_list.extend([peak[3]] * len(angle_vals))
        

        # Filter out too similar angles (within <±30 degrees)
        deduplicated_angles = []
        deduplicated_scores = []
        for angle, score in zip(angle_list, score_list):
            # Only add if not too similar to any existing angle
            if not any(abs(angular_diff(angle, existing)) < 30 for existing in deduplicated_angles):
                deduplicated_angles.append(angle)
                deduplicated_scores.append(score)

        total_angles = len(deduplicated_angles)
        angle_map[bond_idx] = [name, atom_indices, deduplicated_angles]
        score_map[bond_idx] = deduplicated_scores
        total_confs *= total_angles


    return total_confs, angle_map, score_map, rot_bonds

def within_tolerance(angle, center, tolerance):
    """
    Check if an angle falls within a specified tolerance range around a central angle in the [-180, 180] range.

    Args:
    angle (float): The angle to check.
    center (float): The central angle.
    tolerance (float): The tolerance range.

    Returns:
    bool: True if the angle is within the tolerance range, False otherwise.
    """
    # Normalize angle and center to the [-180, 180] range
    angle = (angle + 180) % 360 - 180
    center = (center + 180) % 360 - 180
    
    # Calculate the difference considering wrap-around
    diff = (angle - center + 180) % 360 - 180
    return abs(diff) <= tolerance

