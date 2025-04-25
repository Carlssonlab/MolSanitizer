from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolAlign
from rdkit.Chem.rdchem import Mol, Conformer

import yaml, subprocess
import os, shutil
import itertools
import copy
import numpy as np
import pandas as pd
import random
import logging

from scipy.spatial.distance import pdist, squareform
from collections import defaultdict
from pathlib import Path

from EirVS.filtering import strain_filter
from EirVS.db2 import mol2db2, mol2


logger = logging.getLogger('eirvs')

# Define SMARTS patterns for various functional groups
sulfonamide_like_substructure = Chem.MolFromSmarts("[*:1][S;$(S(=*)=*):2]-!@[N&+0;!$([NH2]):3](-[*,#1:4])-[*,#1:5]")
aliphatic_nitrogen_substructure = Chem.MolFromSmarts("[A:1]@[N&+0;!$(N-*=*):2](@[A:3])!@[*,#1:4]")
conjugated_substituted_nitrogen = Chem.MolFromSmarts('[a:1]:[a:2]:[a:3]:[nX3&+0:4]-*')
additional_substituted_nitrogen = Chem.MolFromSmarts('*-[nX3&+0:1]:[a:2]:[a:3]')
amide_substructure = Chem.MolFromSmarts('[O:1]=[CX3:2]!@[N&+0:3](-[!#1:4])-[#1:5]')
symmetric_ring = Chem.MolFromSmarts('*!@-a1[cH][cH][a][cH][cH]1')

Torlib = strain_filter.parse_torlib()
rigid_rule_files = Path(__file__).parent.parent / 'Data' / 'rigid_part_rules.txt'
rigid_rules = pd.read_csv(rigid_rule_files, header=None, sep =r'\s+', names=['SMARTS','label'])
rigid_rules['mol'] = rigid_rules['SMARTS'].apply(lambda x: Chem.MolFromSmarts(x))
rotatable_pattern=r'[*]~[*;!$(*#*)]-&!@[*;!$(*#*)]~[*]'

with open(Path(__file__).parent.parent / 'eirvs_configurations.yaml') as confFile:
    eirvs_configurations = yaml.full_load(confFile)
CORINA_EXE = eirvs_configurations['CORINA']



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
            with open(f"mol2/{name}_corina.mol2", "w") as f:
                f.write(mol2_blocks[0])
        ring_confs = []
        # Now you have a list of strings, each containing one MOL2 molecule
        for i, mol in enumerate(mol2_blocks, 1):
            rdkit_mol = Chem.MolFromMol2Block(mol, removeHs=False, sanitize=True)
            if rdkit_mol:
                ring_confs.append(rdkit_mol)
        mol2_string = mol2_blocks[0]
        return mol2_string, ring_confs
    
def find_symmetric_rings(mol_H: Mol):
    '''
    Find the symmetric rings in the molecule. *!@-a1[cH][cH][a][cH][cH]1
    '''
    return mol_H.GetSubstructMatches(symmetric_ring)

def find_flipped_nitrogen(mol_H: Mol):
    '''
    Find the flippable nitrogen in the molecule. n(R):c:c:c:c.
    '''
    return mol_H.GetSubstructMatches(aliphatic_nitrogen_substructure)

def find_conjugated_substituted_nitrogen1(mol_H: Mol):
    '''
    Find the conjugated substituted nitrogen in the molecule. c:c:n(R):c:c. 

    Only match to the atoms within the same ring as n.
    '''
    ring_info = mol_H.GetRingInfo()
    matches = [tuple(match[:4]) for match in mol_H.GetSubstructMatches(conjugated_substituted_nitrogen)]
    filtered_matches = []
    for match in matches:
        for ring in ring_info.AtomRings():
            if set(match).issubset(ring):
                filtered_matches.append(match)
    return filtered_matches

def find_conjugated_substituted_nitrogen2(mol_H: Mol):
    '''
        Another function to find conjugated substituted nitrogen in the molecule fo fix the dihedral

        Find two *-n:a:a matches for each Ns, then fix them to 180 to make them planar
    '''
    matches = mol_H.GetSubstructMatches(additional_substituted_nitrogen)
    ring_info = mol_H.GetRingInfo()
    filtered_matches = []
    for match in matches:
        for ring in ring_info.AtomRings():
            if set(match[1:]).issubset(ring):
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

def normalize_scores(peaks):
    """Normalize the scores for each peak in a list of peaks."""
    total_score = sum(peak[-1] for peak in peaks)
    return [(angles[:-1], (peak[-1] / total_score) * 100) for angles, peak in zip(peaks, peaks)]

def generate_combinations(data):
    """Generate all combinations of peaks from different dihedrals and calculate the combined score."""
    normalized_data = []
    combination_indices = []
    
    for smarts, atoms, peaks in data:
        normalized_peaks = normalize_scores(peaks)
        normalized_data.append(normalized_peaks)
        combination_indices.append(range(len(peaks)))

    # Generate all possible combinations of normalized peaks
    all_combinations = itertools.product(*combination_indices)
    
    # Create a dictionary to store the results
    combinations_dict = defaultdict(float)
    
    for combination in all_combinations:
        combined_probability = 1
        for dihedral_idx, peak_idx in enumerate(combination):
            combined_probability *= normalized_data[dihedral_idx][peak_idx][-1]
        combinations_dict[combination] = combined_probability
    
    return combinations_dict

def get_flexible_ring(mol: Mol, planar_lib:list, non_planar_lib:list):
    """Extract the flexible rings from the molecule based on the given rules 
    (planar_lib and non_planar_lib). Only consider the rings from SSR to exclude fused rings.

    Args:
        mol (Mol): Rdkit molecule object
        planar_lib and non_planar_lib (list): Lists containing 
            (name, smarts, rdkit mol from smarts, reference set of dihedral)

    Returns:
        planar_rings (set): Set containing the indices of the planar rings
        non_planar_rings (list): List containing the indices of the non-planar rings
    """
    ssr = [set(ring) for ring in Chem.GetSymmSSSR(mol)]
    planar_rings = set()
    for rule in planar_lib:
        matches = mol.GetSubstructMatches(rule[2])
        if len(matches) > 0:
            for match in matches:
                planar_rings.add((match))
    # Pre-sort the planar_rings and ssr sets
    sorted_planar_rings = {tuple(sorted(t)) for t in planar_rings}
    sorted_ssr = {tuple(sorted(t)) for t in ssr}         
    non_planar_rings = []
    set_non_planar_rings= set()
    for rule in non_planar_lib:

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
    
def find_sulfonamide_like_scaffolds(mol_H: Mol):
    """Find all Sulfonamide-like scaffolds (S(O2)-N(R1)R2 or (S(O)(N)-N(R1)(R2))."""
    preliminary_sulfonamide = mol_H.GetSubstructMatches(sulfonamide_like_substructure)
    matches_sulfonamide = []
    for (a, b, c, d, e) in preliminary_sulfonamide:
        if identical_substituents(mol_H, b, c, d, e): continue
        else: matches_sulfonamide.append((a, b, c, d, e))
    return matches_sulfonamide

def classify_confs(conf, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df, tolerance=20):
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
    flippable_N_descriptors = [1 if rdMolTransforms.GetDihedralDeg(conf, *flippable_N) > 0 else 0 for flippable_N in flippable_Ns]
    temp_dict['flippable_N_descriptors'] = flippable_N_descriptors if flippable_N_descriptors else [-1]

    # Process sulfo matches
    sulfo_descriptors = tuple([1 if rdMolTransforms.GetDihedralDeg(conf, d, b, c, e) > 0 else 0 for (a, b, c, d, e) in sulfo_matches])
    temp_dict['sulfo_descriptors'] = sulfo_descriptors if sulfo_descriptors else [-1]

    # Create dataframe and append to conf_ring_descriptors_df
    temp_df = pd.DataFrame([temp_dict])
    conf_ring_descriptors_df = pd.concat([conf_ring_descriptors_df, temp_df], ignore_index=True)

    return conf_ring_descriptors_df

def remove_unfavorable_confs(conf_ring_descriptors_df: pd.DataFrame, name: str)-> pd.DataFrame:
    for column in conf_ring_descriptors_df.columns[2:-2]:
        if (conf_ring_descriptors_df[column] == -1).all():
            conf_ring_descriptors_df.drop(columns=[column], inplace=True)
            #print(f"\t While handling {name}, found no good conformation of ring: {column}")

    # Step 2: Drop rows where any value in remaining columns (2 onwards) is -1
    for column in conf_ring_descriptors_df.columns[2:-2]:
        conf_ring_descriptors_df = conf_ring_descriptors_df[conf_ring_descriptors_df[column] != -1]
    return conf_ring_descriptors_df

def get_sdf_mol2_filename(name: str, rigid_scaffold_idx: int, align_copy: int):
    return f"{name}_mol{rigid_scaffold_idx}_align{align_copy}.sdf", f"{name}_mol{rigid_scaffold_idx}_align{align_copy}.mol2"

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

def getDihedralMatches(mol, pattern):
    '''return list of atom indices of dihedrals'''
    #this is rdkit's "strict" pattern
    qmol = Chem.MolFromSmarts(pattern)
    matches = mol.GetSubstructMatches(qmol)
    #these are all sets of 4 atoms, uniquify by middle two
    uniqmatches = []
    seen = set()
    for (a,b,c,d) in matches:
        if ((b,c) not in seen) and ((c,b) not in seen):
            seen.add((b,c))
            uniqmatches.append((a,b,c,d))
    return uniqmatches

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


def get_neighboring_atoms(mol, atom_index, excluded_indices):
    """
    Get the neighboring atom symbols and bond counts of a given atom index in the molecule,
    excluding the specified indices.

    Args:
    mol (rdkit.Chem.Mol): The RDKit molecule object.
    atom_index (int): The index of the atom.
    excluded_indices (list): List of atom indices to exclude.

    Returns:
    list: A list of tuples, each containing the neighboring atom symbol and bond count.
    """
    atom = mol.GetAtomWithIdx(atom_index)
    neighbors = [(neighbor.GetSymbol(), neighbor.GetDegree()) for neighbor in atom.GetNeighbors() if neighbor.GetIdx() not in excluded_indices]
    return neighbors
    
def is_terminal(mol, atom_indices):
    """
    Check if either of the two middle atoms in the dihedral is terminal.

    Args:
    mol (rdkit.Chem.Mol): The RDKit molecule object.
    atom_indices (list): List of two middle atom indices defining the dihedral.

    Returns:
    bool: True if the bond is terminal, False otherwise.
    """
    second_atom_neighbors = get_neighboring_atoms(mol, atom_indices[0], atom_indices)
    third_atom_neighbors = get_neighboring_atoms(mol, atom_indices[1], atom_indices)

    # Get the atom type (symbol) of the neighboring atoms
    second_terminal = all(bond_count == 1 for _, bond_count in second_atom_neighbors)
    third_terminal = all(bond_count == 1 for _, bond_count in third_atom_neighbors)
    return(second_terminal) or (third_terminal)

def is_symmetric(mol, atom_indices):
    """
    Check if either of the two middle atoms in the dihedral is terminal."""
    second_atom_neighbors = get_neighboring_atoms(mol, atom_indices[0], atom_indices)
    third_atom_neighbors = get_neighboring_atoms(mol, atom_indices[1], atom_indices)
    
    # Get the atom type (symbol) of the neighboring atoms
    second_atom_symbols = [symbol for symbol, _ in second_atom_neighbors]
    third_atom_symbols = [symbol for symbol, _ in third_atom_neighbors]
    # *-C(C3) should not be considered as symmetric
    if (second_atom_symbols.count('C') == 3)  or (third_atom_symbols.count('C') == 3):
        return False
    return(len(set(second_atom_symbols)) == 1 and len(second_atom_neighbors)==3) or (len(set(third_atom_symbols))==1 and len(third_atom_neighbors)==3)
   


def count_confs_by_rotbonds(mol, ignoreTorlib=False, VERBOSE=False):
    """
    Count the number of conformations based on rotatable bonds and reorder bonds with terminal
    atoms at the beginning.

    Args:
    mol (rdkit.Chem.Mol): The RDKit molecule object.
    VERBOSE (bool): If True, prints detailed information.

    Returns:
    tuple: Number of conformations, reordered rotatable bonds, and matched torsion rules.
    """
    rotatable_bonds = getDihedralMatches(mol, rotatable_pattern)

    match_torlib = strain_filter.get_match_dihedral(mol, Torlib)

    # Remove some redundant rotation related to the symmetric rings 
    # such as p-substituted benzenes or p-pyridine
    # Should only remove from one side of the ring. 
    # If two sides of a symmetric ring has 2 x 2 possible peaks that are 180 degrees apart, 
    # we should produce 2 possible combinations, 0x180 and 0x0, rather than produce
    # 180x180 0x180 0x0 180x0 as they are symmetrically equivalent.

    # symmetric_rings = smi2db2_utils.find_symmetric_rings(mol)
    # #print(f"\tSymmetric rings: {symmetric_rings}")
    # reduced_rotated_rings = set()
    # for ring in symmetric_rings:
    #     ring_key = set(ring[0:2])
    #     sorted_ring_atoms = tuple(sorted(ring[1:7]))  # Create a sorted tuple for consistent representation
    #     # Already processed this symmetric ring; skip
    #     if sorted_ring_atoms in reduced_rotated_rings: continue
    #     reduced_rotated_rings.add(sorted_ring_atoms)
    #     for rule in match_torlib:
    #         if ring_key == set(rule[1][1:3]):
    #             # Filter out redundant rotations where angles differ by 180 degrees
    #             rule[2] = [
    #                 angle for i, angle in enumerate(rule[2])
    #                 if all((angle[0] - other_angle[0]) % 180 != 0 for other_angle in rule[2][i + 1:])
    #             ]
    #             break

    for rotatable_bond in rotatable_bonds:
        if is_terminal(mol, rotatable_bond[1:3]) and is_symmetric(mol, rotatable_bond[1:3]):
        # Exclude redundant angles for symmetric terminal rotatable bonds (e.g., CH3, CF3)
            bond_key = set(rotatable_bond[1:3])
            for rule in match_torlib:
                if bond_key == set(rule[1][1:3]):
                    # Retain only unique angles that differ by at least 120 degrees
                    rule[2] = [
                        angle for i, angle in enumerate(rule[2])
                        if all((angle[0] - other_angle[0]) % 120 != 0 for other_angle in rule[2][i + 1:])
                    ]
                    while len(rule[2])>2: rule[2].pop()
                    break
                            
    if ignoreTorlib:
        amide_linkages = find_amide(mol)
        amide_atoms=set()
        for a, b, c, d, e in amide_linkages: 
            amide_atoms.add(b)
            amide_atoms.add(c)
        const_rule = [(-120, 30, 30, 1), (-60, 30, 30, 1), (0, 30, 30, 1), (60, 30, 30, 1), (120, 30, 30, 1), (180, 30, 30, 1)]

    num_confs = 1
    match_torlib_clean = []
    for bond in rotatable_bonds:
        for rule in match_torlib:
            if set(bond[1:3]) == set(rule[1][1:3]):
                if ignoreTorlib:
                    new_rule = copy.deepcopy(rule[:2])
                    if set(bond[1:3]).issubset(amide_atoms):
                        match_torlib_clean.append(rule)
                        num_confs *= (len(rule[2]))
                    else:
                        new_rule.append(const_rule)
                        num_confs *= (len(new_rule[2]))
                        match_torlib_clean.append(new_rule)
                else:
                    num_confs *= (len(rule[2]))
                    match_torlib_clean.append(rule)
                break

    if VERBOSE: 
        print(f"\t{rotatable_bonds}")
        for i in match_torlib_clean: print(f"\t{i}")

    return num_confs, match_torlib_clean

def get_random_angle(mean, tolerance):
    """
    Generate a random angle value from a Gaussian distribution given the expected mean and standard deviation,
    and limit it within the specified range. Normalize the result to the [-180, 180] degree range.

    Args:
    mean (float): The expected mean angle.
    std_dev (float): The standard deviation.
    lower_limit (float): The lower limit for the angle.
    upper_limit (float): The upper limit for the angle.
    seed (int, optional): The seed for the random number generator.

    Returns:
    float: A random angle normalized to the [-180, 180] degree range.
    """

    # Generate a random angle within the specified Gaussian distribution and range limits
    if tolerance == 0: return mean
    while True:
        random_angle = random.gauss(mean, tolerance)
        if mean-tolerance < random_angle < mean+tolerance: 
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
    with open('eirvs_error.log', 'a') as f:
        f.write(f"{smiles} \t {name}\n")

def Align_ConvertToDb2(ring_conf, rigid_scaffold, solv_obj, name, smiles, longname):
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