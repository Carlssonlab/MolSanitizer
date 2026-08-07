import subprocess
import os
import shutil
import logging
import time

import numpy as np
import yaml

from pandas import DataFrame, concat, read_csv  # only what you use
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolAlign
from rdkit.Chem.rdchem import Mol, Conformer

from pathlib import Path

from msani.db2 import mol2db2, mol2


logger = logging.getLogger('msani')

# Define SMARTS patterns for various functional groups
sulfonamide_like_substructure = Chem.MolFromSmarts("[*:1][S;$(S(=*)=*):2]-!@[N&+0;!$([NH2])!r:3](-[*,#1;!$(C=A):4])-[*,#1;!$(C=A):5]")
sulfonamide_cycloheptane = Chem.MolFromSmarts("[S;$(S(=*)=*)]-!@[N&+0&r7;!$([NH2])]1-[*;!$(C=A)]-[*]~[*]~[*]~[*]~[*;!$(C=A)]1")
substituted_C_cyclohexane = Chem.MolFromSmarts('[!#1;!$(*-!@[CH]1-[*^2;!O]~[*^2;!O]~*~[*^2;!O]~[*^2;!O]-1)]-!@[CH]1-[*]~[*]~[*]~[*]-[A]-1') # Ignore check for theoretically planar cyclohexanes
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
aliphatic_hydroxyl_thiol = Chem.MolFromSmarts('C-[OX2H,SX2H]') #aliphatic hydroxyls and thiols
phenol_thiolphenol = Chem.MolFromSmarts('a-[OX2H,SX2H]') #phenol and thiophenol
hydroxamic_acid = Chem.MolFromSmarts('[NX3;$(N(-C=O))]-[OX2H1]') #hydroxamic acid and its tautomeric aci form
hydroxyl_amine = Chem.MolFromSmarts('[NX3;!$(N~*=[O,S])!$(N=,#*):1]!@[OX2H1:2]') #hydroxylamine and its tautomeric aci form

const_rule = [(-120, 30, 30, 1), (-60, 30, 30, 1), (0, 30, 30, 1), (60, 30, 30, 1), (120, 30, 30, 1), (180, 30, 30, 1)]


rigid_rule_files = Path(__file__).parent.parent / 'Data' / 'rigid_part_rules.txt'
rigid_rules = read_csv(rigid_rule_files, header=None, sep =r'\s+', names=['SMARTS','label'])
rigid_rules['mol'] = rigid_rules['SMARTS'].apply(lambda x: Chem.MolFromSmarts(x))
rotatable_pattern=r'*~[!$(*(#*)-!@*#*)&!D1]-!@[!$(*(#*)-!@*#*)&!D1]~*'
rot_bond_mol = Chem.MolFromSmarts(rotatable_pattern)

with open(Path(__file__).parent.parent / 'msani_configurations.yaml') as confFile:
    msani_configurations = yaml.safe_load(confFile)
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
        "-o", "t=sdf",
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
        raw_blocks = cleaned_output.split("$$$$")

        # Add back the header to each block (except the first if it's empty)
        sdf_blocks = [
            block.strip() + "\n$$$$\n" 
            for block in raw_blocks if block.strip()
        ]
        if VERBOSE:
            print(f"\tNumber of ring conformers: {len(sdf_blocks)}")
        ring_confs = []
        # Now you have a list of strings, each containing one MOL2 molecule
        for i, mol in enumerate(sdf_blocks, 1):
            rdkit_mol = Chem.MolFromMolBlock(mol, removeHs=False, sanitize=True)
            if rdkit_mol:
                ring_confs.append(rdkit_mol)
        return ring_confs
    

def find_flipped_nitrogen(mol_H: Mol):
    '''
    Find the flippable nitrogen in the molecule. Mainly for reinforcing the equatorial of the *-N of piperidine and piperazine
    '''
    return mol_H.GetSubstructMatches(flippable_Ns_1) + mol_H.GetSubstructMatches(flippable_Ns_2)

def find_flipped_carbon(mol_H: Mol):
    '''
    Find the flippable carbon in the molecule. Mainly for substituted cyclohexane
    '''
    return mol_H.GetSubstructMatches(substituted_C_cyclohexane) 
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
        
def find_sulfonamide_cycloheptane(mol_H: Mol):
    '''
    Find the sulfonamide in the molecule with 7-membered ring
    '''
    return mol_H.GetSubstructMatches(sulfonamide_cycloheptane)

def find_amide(mol_H: Mol):
    '''
        Find the amide in the molecule. O=C-N(-R)-H
    '''
    return mol_H.GetSubstructMatches(amide_substructure)


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
    Determines if a substituent is in an equatorial position on a ring.

    This function calculates two dihedral angles around a specific atom in the molecule
    and checks if they both fall within the range that indicates an equatorial position 
    (approximately 150-180 degrees in absolute value).

    Parameters
    ----------
    conf : rdkit.Chem.rdchem.Conformer
    atom_idx : list or tuple
        List of atom indices defining the relevant atoms for dihedral calculations.
        Requires at least 6 indices:
        - atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3]: First dihedral angle
        - atom_idx[0], atom_idx[1], atom_idx[-1], atom_idx[-2]: Second dihedral angle
        
    Returns
    -------
    bool
        True if the substituent is equatorial (both dihedrals between 150-180 degrees),
        False otherwise
    """
    
    dihedral1 = rdMolTransforms.GetDihedralDeg(conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3])
    dihedral2 = rdMolTransforms.GetDihedralDeg(conf, atom_idx[0], atom_idx[1], atom_idx[-1], atom_idx[-2])
    if 140 <= abs(dihedral1) <= 180 and 140 <= abs(dihedral2) <= 180:
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
                    sulfo_7_ring=(), 
                    tolerance=25):
    """
    Classify a conformer based on ring conformations, flippable nitrogens,
    substituted cyclohexanes, and sulfonamide-like scaffolds.
    """
    
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
    temp_dict['equatorial_subs_Ns'] = flippable_N_descriptors if flippable_Ns else -1

    # Process substituted cyclohexane
    aliphatic_cyclohexane_descriptors = sum([1 if is_equatorial(conf, atom_idx) else 0 for atom_idx in flippable_Cs])
    temp_dict['equatorial_subs_Cs'] = aliphatic_cyclohexane_descriptors if flippable_Cs else -1
    
    # Process sulfo matches
    sulfo_descriptors = tuple([1 if rdMolTransforms.GetDihedralDeg(conf, d, b, c, e) > 0 else 0 for (a, b, c, d, e) in sulfo_matches])
    temp_dict['sulfo_descriptors'] = sulfo_descriptors if sulfo_descriptors else [-1]

    # Process flippable Nitrogens in sulfonamides within 7 membered rings, 
    # we want both axial and equatorials
    sulfo_7_descriptors = sum([1 if is_equatorial(conf, atom_idx) else 0 for atom_idx in sulfo_7_ring])
    temp_dict['sulfo_7_descriptors'] = sulfo_7_descriptors if sulfo_7_ring else -1
    return temp_dict

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
    # Heuristic: if the query contains SMARTS-specific characters, do not treat it as a standard SMILES
    if any(char in query for char in ['~', '*', '$', '!', '&', ',']) or '[#' in query:
        return query
    mol = Chem.MolFromSmiles(query)
    if mol:
        return Chem.CanonSmiles(query)
    else: 
        return query  # Return SMARTS or invalid input unchanged




def count_confs_by_rotbonds(mol,
                            rot_bonds,
                            amide_bonds,
                            ignoretorlib=False,
                            torlib=None,
                            VERBOSE=False):
    """
    Count the number of conformations based on rotatable bonds and torsion rules.

    This function analyzes the molecule's rotatable bonds, applies torsion rules (optionally ignoring the torsion library),
    handles special cases for amides and symmetric patterns, and returns the total number of conformations along with
    the processed torsion rules and heteroatom hydrogen bonds.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object.
        rot_bonds (list): List of rotatable bonds, each represented as a tuple of atom indices.
        amide_bonds (list): List of tuples representing amide linkages.
        ignoretorlib (bool): If True, ignore the torsion library for non-amide bonds and use a constant rule.
        torlib: Torsion library object to match dihedral angles.
        VERBOSE (bool): If True, print detailed information about symmetry and filtering.

    Returns:
        total_confs (int): Estimated total number of conformations.
        match_torlib_clean (list): List of processed torsion rules for each bond.
        hetero_H_bonds (list): List of rotatable bonds involving heteroatom hydrogens.
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
    hetero_H_bonds = _process_hetero_hydrogen_bonds(mol, bond_to_rule)

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
    
    # Reorder bonds and create final mappings
    rot_bonds_reordered = _reorder_bonds(rot_bonds, hetero_H_bonds)
    bond_to_rule_reordered = {b: bond_to_rule[b] for b in rot_bonds_reordered}

    total_confs = 1
    
    match_torlib_clean = []
    for rule in bond_to_rule_reordered.values():
        match_torlib_clean.append(rule)
        total_angles = 0
        for peak in rule[2]:
            total_angles += (peak[2]*2/30 + 1)
        total_confs *= total_angles
    if VERBOSE:
        print("\nFinal processed torsion rules:")
        for rule in match_torlib_clean:
            print(f"\t{rule}")

    return int(total_confs), match_torlib_clean, hetero_H_bonds


def find_rigid_part(mol, request_alignment=None):
    '''Find fused ring system of the molecule as the rigid part, 
    if none, use the hierarchical rules in rigid_part_rules.txt'''

    rigid_part = []
    rule_label = None
    
    # If the user request for only a specific ring as rigid segment.
    if request_alignment:
        matches = mol.GetSubstructMatches((request_alignment))
        if len(matches) != 0:
            if len(matches) > 1: 
                logger.warning(f"Multiple matches found for the requested alignment: {Chem.MolToSmiles(Chem.RemoveHs(mol))}. Will align by both.")
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

    global rot_bond_mol
    matches = mol.GetSubstructMatches(rot_bond_mol)
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
    n_steps = int(tolerance / step)
    angles = [typical + i * step for i in range(-n_steps, n_steps + 1)]
    #angles = [typical, typical - step, typical + step] #Only sample 3 angles for each dihedral

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

def _zero_out_fluctuations(rule):
    """Helper function to zero out fluctuations in a rule."""
    adjusted_peaks = []
    for peak in rule[2]:
        if isinstance(peak, tuple) and len(peak) >= 3:
            adjusted_peaks.append((peak[0], 0, 0) + peak[3:])
        else:
            adjusted_peaks.append(peak)
    rule[2] = adjusted_peaks

def _process_hetero_hydrogen_bonds(mol, bond_to_rule):
    """Process hydroxyl and thiol bonds with special angle rules."""
    hetero_H_bonds = []
    
    # Define bond patterns and their corresponding angle rules
    bond_patterns = [
        (aliphatic_hydroxyl_thiol, [(0, 0, 0, 1), (120, 0, 0, 1), (-120, 0, 0, 1)]),
        (phenol_thiolphenol, [(0, 0, 0, 1), (180, 0, 0, 1)]),
        (hydroxamic_acid, [(0, 0, 0, 1), (180, 0, 0, 1)]),
        (hydroxyl_amine, [(-120, 0, 0, 1), (120, 0, 0, 1)])
    ]
    
    for pattern, angles in bond_patterns:
        for match in mol.GetSubstructMatches(pattern):
            bond_key = tuple(sorted(match))
            if bond_key in bond_to_rule:
                bond_to_rule[bond_key][2] = list(angles)
                hetero_H_bonds.append(bond_key)
    
    return hetero_H_bonds

def _reorder_bonds(rot_bonds, hetero_H_bonds):
    """Reorder bonds to put hetero-H bonds at the end."""
    return ([b for b in rot_bonds if b not in hetero_H_bonds] + 
            [b for b in rot_bonds if b in hetero_H_bonds])

def _extract_angles_from_peaks(peak_list):
    """Extract angles and scores from peak list."""
    angle_list = []
    score_list = []
    
    for peak in peak_list:
        if len(peak) < 4:
            # Fallback case (e.g., symmetric angle filtering already applied)
            angle_list = [p[0] for p in peak_list]
            score_list = [p[1] for p in peak_list]
            break
        
        angle_vals = discretinize_dihedrals(peak[0], peak[2])
        angle_list.extend(angle_vals)
        score_list.extend([peak[3]] * len(angle_vals))
    
    return angle_list, score_list


def _deduplicate_angles(angle_list, score_list):
    """Remove angles that are too similar (within ±30 degrees)."""
    deduplicated_angles = []
    deduplicated_scores = []
    
    for angle, score in zip(angle_list, score_list):
        # Only add if not too similar to any existing angle
        is_similar = any(
            abs(angular_diff(angle, existing)) < 30 
            for existing in deduplicated_angles
        )
        if not is_similar:
            deduplicated_angles.append(angle)
            deduplicated_scores.append(score)
    
    return deduplicated_angles, deduplicated_scores

def _generate_angle_mappings(bond_to_rule_reordered):
    """Generate final angle and score mappings."""
    angle_map = {}
    score_map = {}
    total_confs = 1
    
    for bond_idx, rule in enumerate(bond_to_rule_reordered.values()):
        name, atom_indices, peak_list = rule
        
        # Extract angles and scores from peaks
        angle_list, score_list = _extract_angles_from_peaks(peak_list)
        
        # Remove duplicate angles
        deduplicated_angles, deduplicated_scores = _deduplicate_angles(angle_list, score_list)
        
        # Update mappings
        angle_map[bond_idx] = [name, atom_indices, deduplicated_angles]
        score_map[bond_idx] = deduplicated_scores
        total_confs *= len(deduplicated_angles)
    
    return angle_map, score_map, total_confs

def count_confs_by_rotbonds_v2(mol,
                               rot_bonds,
                               ignoretorlib = False,
                               amide_bonds = None,
                               torlib = None,
                               VERBOSE=False):
    """
    Count the number of conformations based on rotatable bonds and torsion rules (deterministic version).

    This function analyzes the molecule's rotatable bonds, applies torsion rules (optionally ignoring the torsion library),
    handles special cases for amides, symmetric patterns, and heteroatom hydrogens, and returns the total number of
    conformations along with detailed angle and score mappings for each bond.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object.
        rot_bonds (list): List of rotatable bonds (tuples of atom indices).
        ignoretorlib (bool): If True, ignore the torsion library for non-amide bonds and use a constant rule.
        amide_bonds (list): List of tuples representing amide bonds (optional).
        torlib: Torsion library object for matching dihedral rules.
        VERBOSE (bool): If True, print detailed information about symmetry and filtering.

    Returns:
        total_confs (int): Estimated total number of conformations.
        angle_map (dict): Mapping of bond indices to [name, atom_indices, list of angles].
        score_map (dict): Mapping of bond indices to list of scores for each angle.
        rot_bonds_reordered (list): Rotatable bonds reordered (hetero-H bonds at the end).
        hetero_H_bonds (list): List of rotatable bonds involving heteroatom hydrogens.
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
            rule_copy = list(rule_by_bond[bond])
            
            # Only apply ignoretorlib if not amide bonds
            if ignoretorlib and len(set(bond) & amide_atoms) <= 1:
                rule_copy[2] = const_rule
                
            bond_to_rule[bond] = rule_copy

    # Step 2: Handle planar substructures (amides, amidines - zero out fluctuations)
    if (amide_bonds):
        for amide_match in amide_bonds:
            bond_key = tuple(sorted(amide_match[0:2]))
            if bond_key in bond_to_rule:
                _zero_out_fluctuations(bond_to_rule[bond_key])

    # Zero out fluctuations for primary amidines, guanidines:
    uniq_matches_prim_amidines_guanidines = set()
    for patt in prim_amidines_guanidines_pattern_mol:
        for match in mol.GetSubstructMatches(patt):
            uniq_matches_prim_amidines_guanidines.add(tuple(sorted(match[1:3])))

    for bond_key in list(uniq_matches_prim_amidines_guanidines):
        if bond_key in bond_to_rule: 
            _zero_out_fluctuations(bond_to_rule[bond_key])

    # Step 3: Special treatment for hydroxyls/thiols
    hetero_H_bonds = _process_hetero_hydrogen_bonds(mol, bond_to_rule)


    # Step 4: Apply symmetry filtering
    _apply_symmetry_filtering(mol, bond_to_rule, rot_bonds, VERBOSE)
    
    # Step 5: Reorder bonds and create final mappings
    rot_bonds_reordered = _reorder_bonds(rot_bonds, hetero_H_bonds)
    bond_to_rule_reordered = {b: bond_to_rule[b] for b in rot_bonds_reordered}
    
    # Step 6: Generate angle and score mappings
    angle_map, score_map, total_confs = _generate_angle_mappings(bond_to_rule_reordered)
    
    return total_confs, angle_map, score_map, rot_bonds_reordered, hetero_H_bonds

def _apply_symmetry_filtering(mol, bond_to_rule, rot_bonds, VERBOSE):
    """Apply symmetry filtering to avoid redundant angles."""
    for _, sym_row in symmetric_patterns_df.iterrows():
        matches = mol.GetSubstructMatches(sym_row['mol'])
        if not matches:
            continue
            
        if VERBOSE:
            print(f"\tFound symmetric pattern: {sym_row['name']}")
            
        period = 360 / sym_row['num_scaled']
        check_symmetric_using_dfs = sym_row['name'].startswith('general')
        
        for match in matches:
            if check_symmetric_using_dfs:
                idx2, idx3, idx4, idx5 = match[0], match[1], match[2], match[-1]
                if not identical_substituents(mol, idx2, idx3, idx4, idx5):
                    continue
            
            bond_key = tuple(sorted(match[:2]))
            if bond_key not in bond_to_rule:
                continue

            # Handle complete rotation removal (period == 1)
            if period == 1:
                bond_to_rule.pop(bond_key)
                rot_bonds.remove(bond_key)
                continue

            # Handle special case for tri_large_halogeno_methyl
            if sym_row['name'] == 'tri_large_halogeno_methyl':
                _handle_tri_halogeno_methyl(bond_to_rule[bond_key])
                continue
            
            # Apply symmetry filtering
            _filter_symmetric_rule(bond_to_rule[bond_key], period, VERBOSE, bond_key)

def _handle_tri_halogeno_methyl(rule):
    """Handle special case for tri_large_halogeno_methyl pattern."""
    peaks = rule[2]
    angles = [normalize_angle(peaks[0][0] + x) for x in [0, -30, 30]]
    scores = [peaks[0][3], round(peaks[0][3]/2, 2), round(peaks[0][3]/2, 2)]
    rule[2] = list(zip(angles, scores))


def _filter_symmetric_rule(rule, period, VERBOSE, bond_key):
    """Apply symmetric filtering to a torsion rule."""
    if VERBOSE:
        print(f'\tRemove duplicated rotation for {bond_key}')
        
    # Skip if already processed
    if len(rule[2][0]) < 4:
        return
    
    angles, scores = [], []
    for peak in rule[2]:
        angle_list = discretinize_dihedrals(peak[0], peak[2])
        angles.extend(angle_list)
        scores.extend([peak[3]] * len(angle_list))
    
    filtered_angles, filtered_scores = filter_symmetric_angles(angles, scores, period, 10)
    rule[2] = list(zip(filtered_angles, filtered_scores))


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


def check_timeout(start_time, max_duration):
    """Check if the elapsed time has exceeded the maximum duration."""
    elapsed_time = time.time() - start_time
    return elapsed_time > max_duration