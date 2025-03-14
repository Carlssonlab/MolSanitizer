"""Based partly on: 
    https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py
    http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
    http://pubs.acs.org/doi/abs/10.1021/ci2004658
    https://greglandrum.github.io/rdkit-blog/posts/2024-02-11-more-multithreading.html

    Should try to sample all possible conformations based on dihedral angles sampling based on: https://github.com/dkoes/rdkit-scripts/blob/master/rdallconf.py
"""
# Author: Thua-Phong Lam, Jens Carlsson lab, Uppsala University

from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolTransforms, rdDistGeom, rdMolAlign

from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from EirVS.amsol import run_amsol
from EirVS.db2 import mol2db2, mol2, solv
from . import smi2db2_utils
from EirVS.filtering import strain_filter

import random
import time
import multiprocessing
import tarfile, io
import platform
import logging
import os,  shutil
import copy
import pandas as pd
import itertools
logger = logging.getLogger('eirvs')

#rotatable_pattern=r'[*]~[*;!$(*#*)!$([!#6&X2H])!$([!#6&X3H2])]-&!@[*;!$(*#*)!$([!#6&X2H])!$([!#6&X3H2])]~[*]'
rotatable_pattern=r'[*]~[*;!$(*#*)]-&!@[*;!$(*#*)]~[*]'

Torlib = strain_filter.parse_torlib()
rigid_rule_files = Path(__file__).parent.parent / 'Data' / 'rigid_part_rules.txt'
rigid_rules = pd.read_csv(rigid_rule_files, header=None, sep =r'\s+', names=['SMARTS','label'])
rigid_rules['mol'] = rigid_rules['SMARTS'].apply(lambda x: Chem.MolFromSmarts(x))

planar_lib, non_planar_lib = strain_filter.parse_sr_confs_library()

system = platform.system()
if system == 'Windows':
    AMSOLEXE = Path(__file__).parent.parent / "amsol" / "amsol7.1.exe"
elif system == 'Linux':
    AMSOLEXE = Path(__file__).parent.parent / "amsol" / "amsol7.1"
if not AMSOLEXE.exists():
    raise FileNotFoundError(f"AMSOL executable not found at {AMSOLEXE}. Check the amsol directory for instruction to install amsol")

def setup_env():
    env = os.environ.copy()
    script_dir = Path(__file__).parent.parent
    extra_libs_path = script_dir / "libs" / "extralibs-2"

    if 'LD_LIBRARY_PATH' in env:
        env['LD_LIBRARY_PATH'] += f":{script_dir}:{extra_libs_path}"
    else:
        env['LD_LIBRARY_PATH'] = f"{script_dir}:{extra_libs_path}"
    return env

def write_to_file(content, file):
    with open(file, "w") as f:
        f.write(content)

def write_to_tarball(ball, data, name):
    tar = tarfile.TarInfo(name=name)
    tar.size = len(data)
    ball.addfile(tar, io.BytesIO(data))

def convert(data, inf, otf):
    obConversion = ob.OBConversion()
    obMol = ob.OBMol()
    obConversion.SetInAndOutFormats(inf, otf)
    obConversion.ReadString(obMol, data)

    return obConversion.WriteString(obMol)


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
    return(len(set(second_atom_symbols))==1 and len(second_atom_neighbors)==3) or (len(set(third_atom_symbols))==1 and len(third_atom_neighbors)==3)
   

def convert_sdf_mol2(sdf_file, output_mol2, VERBOSE: bool = False):
    
    # Set up OpenBabel conversion
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "mol2")
    # Read all molecules from the SDF file
    mol = ob.OBMol()
    if not obConversion.ReadFile(mol, sdf_file):
        print(f"Error reading SDF file: {sdf_file}")
        return
    molecule_count = 0
    with open(output_mol2, 'w') as out_file:
        while True:
            molecule_count += 1
            mol.SetAutomaticPartialCharge(False)

            # Convert molecule to MOL2 format and get as string
            mol2_str = obConversion.WriteString(mol)
            # Write to output file
            out_file.write(mol2_str)
            
            # Clear the molecule and read the next one
            mol.Clear()
            if not obConversion.Read(mol):
                break
    
    if VERBOSE: print(f"\tConverted and saved {molecule_count} conformations to {output_mol2}.")


def convert_sdf_mol2_str(sdf_io, VERBOSE: bool = False):
    # Set up OpenBabel conversion
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "mol2")
    # Read all molecules from the SDF file
    mol = ob.OBMol()
    sdf_data = sdf_io.getvalue()
    if not obConversion.ReadString(mol, sdf_data):
        print(f"Error reading SDF file: {sdf_data}")
        return
    molecule_count = 0
    output_mol2 = ""
    while True:
        molecule_count += 1
        mol.SetAutomaticPartialCharge(False)

        # Convert molecule to MOL2 format and get as string
        mol2_str = obConversion.WriteString(mol)
        
        # Write to output file
        output_mol2 += mol2_str
        
        # Clear the molecule and read the next one
        mol.Clear()
        if not obConversion.Read(mol):
            break
    
    if VERBOSE: print(f"\tConverted and saved {molecule_count} conformations.")

    return [line + '\n' for line in output_mol2.splitlines()]



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
        amide_linkages = smi2db2_utils.find_amide(mol)
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


def find_rigid_part(mol, rigid_rules):
    '''Find fused ring system of the molecule as the rigid part, 
    if none, use the hierarchical rules in rigid_part_rules.txt'''
    ssr = [set(ring) for ring in Chem.GetSymmSSSR(mol)]
    rigid_part = []
    rule_label = None
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


def embed_smiles_rdkit(smiles, name, randomSeed=42, VERBOSE=False):
    """
    Embeds a SMILES string into molecular conformers using RDKit,
    applying various quality checks and corrections inherited from the MMFF94s force field.
    Parameters:
        smiles (str): The SMILES representation of the molecule.
        name (str): The name to assign to the molecule.
        randomSeed (int, optional): Seed for the random number generator to ensure reproducibility. Defaults to 42.
        VERBOSE (bool, optional): If set to True, enables detailed logging of the embedding process. Defaults to False.
    Returns:
        tuple:
            - amsol_mol (rdkit.Chem.Mol): The molecule object containing the embedded conformers.
            - netcharge (int): The net formal charge of the molecule.
            - rigid_scaffolds (list of rdkit.Chem.Mol): A list of rigid scaffold molecules derived from the conformers.
            - sulfo_matches (list): A list of sulfonamide-like structures identified in the molecule.
    """

    mol_H = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mol_H.SetProp("_Name", name)
    amsol_mol = Chem.Mol(mol_H)
    empty_mol = Chem.Mol(mol_H)

    params = rdDistGeom.srETKDGv3()
    params.numThreads = 1  # Use all available threads
    params.pruneRmsThresh = 0.35  # Prune conformations that are too similar, not user-definable here
    params.randomSeed = randomSeed # For reproducibility
    params.useRandomCoords = True
    num_ring_confs = 1

    ssr = [set(ring) for ring in Chem.GetSymmSSSR(mol_H)]
    planar_rings, non_planar_rings = smi2db2_utils.get_flexible_ring(mol_H, ssr, planar_lib, non_planar_lib)
    sulfo_matches = smi2db2_utils.find_sulfonamide_like_scaffolds(mol_H)
    #In case only 1 ring is output, we don't need more embedding just for the flippable Ns
    if num_ring_confs > 1: flippable_Ns = smi2db2_utils.find_flipped_nitrogen(mol_H)
    else: flippable_Ns = []
    conjugated_substituted_Ns = smi2db2_utils.find_conjugated_substituted_nitrogen1(mol_H)
    additional_conjugated_substituted_Ns = smi2db2_utils.find_conjugated_substituted_nitrogen2(mol_H)
    amide_linkages = smi2db2_utils.find_amide(mol_H)

    if VERBOSE:
        if planar_rings:
            print('\t Found planar rings:')
            print(f'\t{planar_rings}')
        if non_planar_rings:
            print('\t Found non_planar rings:')
            for ring in non_planar_rings: print(f'\t {ring[0]} {ring[1]}')
        if flippable_Ns: 
            print('\tFound flippable N structures')
            for match in flippable_Ns: print(f'\t {match}')
        if sulfo_matches: 
            print('\tFound sulfonamide-like structures')
            for match in sulfo_matches: print(f'\t {match}')
        if conjugated_substituted_Ns:
            print('\tFound conjugated substituted N structures')
            for match in conjugated_substituted_Ns: print(f'\t {match}')
            for match in additional_conjugated_substituted_Ns: print(f'\t {match}')
    
    if sulfo_matches or non_planar_rings or flippable_Ns: numConfs = 100
    else: numConfs = 10

    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(1) #1 means vacumn, 80 means water, 20 is the compromised value (still arbitrary)
    conf_ring_descriptors_df = pd.DataFrame()
    try:
        for cid in rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId=cid)
            if conjugated_substituted_Ns:
                for a, b, c, d in conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5)
                for a, b, c, d in additional_conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 178, 182, 1)
            if amide_linkages:
                for a, b, c, d, e in amide_linkages: # O=C-N(-C)-H should be coplanar.
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 1)
                    ff.MMFFAddTorsionConstraint(a, b, c, e, False, 180, 180, 1)
            ff.Minimize()
            conformer = mol_H.GetConformer(cid)
            energy = ff.CalcEnergy()
            conf_ring_descriptors_df = smi2db2_utils.classify_confs(conformer, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df)
    except: pass
    if len(conf_ring_descriptors_df) == 0:
        # In case where srETKDGv3 failed in embedding the molecule, 
        # we have to use the macrocyclic version.
        logger.warning(f"srETKDGv3 failed for {name}, using macrocyclic version")
        params = rdDistGeom.ETKDGv3()
        params.numThreads = 1  # Use all available threads
        params.pruneRmsThresh = 0.35  # Prune conformations that are too similar, not user-definable here
        params.randomSeed = randomSeed # For reproducibility
        params.useRandomCoords = True
        for cid in rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId=cid)
            if conjugated_substituted_Ns:
                for a, b, c, d in conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5)
                for a, b, c, d in additional_conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 178, 182, 1)
            if amide_linkages:
                for a, b, c, d, e in amide_linkages: # O=C-N(-C)-H should be coplanar.
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 1)
                    ff.MMFFAddTorsionConstraint(a, b, c, e, False, 180, 180, 1)
            ff.Minimize()
            conformer = mol_H.GetConformer(cid)
            energy = ff.CalcEnergy()
            conf_ring_descriptors_df = smi2db2_utils.classify_confs(conformer, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df)

    conf_ring_descriptors_df.sort_values('Energy', inplace=True)
    # Keep a reservoir as the lowest energy possible conformer in case no good ring conformers are found.
    reservoir = conf_ring_descriptors_df.iloc[0, 0]
    # Keep the top 10 conformers for further processing (e.g., AMSOL)
    for idx in range(min(10, len(conf_ring_descriptors_df))):
        amsol_mol.AddConformer(conf_ring_descriptors_df.iloc[idx, 0], assignId=True)
    if VERBOSE:
        print(f"\tamsol_mol contains: {amsol_mol.GetNumConformers()}")
   
            
    #conf_ring_descriptors_df.to_csv(f"{name}_conf_ring_descriptors.csv", index=False)
    conf_ring_descriptors_df = smi2db2_utils.remove_unfavorable_confs(conf_ring_descriptors_df, name)
    

    rigid_scaffolds = []

    if len(conf_ring_descriptors_df) == 0:
        # No conformer is found to compromise all the non-planar rings. 
        # Use the lowest energy conformer``
        scaffold = Chem.Mol(empty_mol)
        scaffold.AddConformer(reservoir, assignId=True)
        rigid_scaffolds.append(scaffold)
        netcharge = sum(atom.GetFormalCharge() for atom in mol_H.GetAtoms())
        return amsol_mol, netcharge, rigid_scaffolds, False
    
    # Process rigid scaffolds based on sulfo descriptors
    
    if sulfo_matches:
        for sulfo_match in conf_ring_descriptors_df['sulfo_descriptors'].unique():
            temp_list = conf_ring_descriptors_df[conf_ring_descriptors_df['sulfo_descriptors'] == sulfo_match].values.tolist()
            initial_len = len(temp_list)
            num_confs_per_regioisomers = 0
            while num_confs_per_regioisomers < num_ring_confs and temp_list:
                lowest_energy_entry = temp_list.pop(0)
                conformer, current_descriptors = lowest_energy_entry[0], lowest_energy_entry[2:-1]
                scaffold = Chem.Mol(empty_mol)
                conf_id = scaffold.AddConformer(conformer, assignId=True)
                rigid_scaffolds.append(Chem.Mol(scaffold, conf_id))
                num_confs_per_regioisomers += 1
                temp_list = smi2db2_utils.ring_conf_clusters(current_descriptors, temp_list)
            if VERBOSE: print(f'\tBefore: {initial_len}, after: {conf_id+1}')


    else:
        temp_list = conf_ring_descriptors_df.values.tolist()
        while len(rigid_scaffolds) < num_ring_confs and temp_list:
            lowest_energy_entry = temp_list.pop(0)
            conformer, current_descriptors = lowest_energy_entry[0], lowest_energy_entry[2:-1]
            scaffold = Chem.Mol(empty_mol)
            scaffold.AddConformer(conformer, assignId=True)
            rigid_scaffolds.append(scaffold)
            temp_list = smi2db2_utils.ring_conf_clusters(current_descriptors, temp_list)
        if VERBOSE: print(f'\tBefore: {len(mol_H.GetConformers())}, after: {len(rigid_scaffolds)}')

    netcharge = sum(atom.GetFormalCharge() for atom in mol_H.GetAtoms())

    return amsol_mol, netcharge, rigid_scaffolds, sulfo_matches


def generate_conformation(queue, smiles, name, randomSeed, VERBOSE):
    """
    Generates initial 3D conformations for a given SMILES string and enqueues the results.
    This function leverages multiprocessing to generate 3D molecular conformations using RDKit.
    It is primarily used to handle operations that may require timeouts.
    Parameters:
        queue (multiprocessing.Queue): The queue to which the results will be put.
        smiles (str): The SMILES string representing the molecule.
        name (str): The name identifier for the molecule.
        randomSeed (int): Seed value for random number generation to ensure reproducibility.
        VERBOSE (bool): If True, prints progress messages.
    Enqueues:
        tuple: A tuple containing:
            - amsol_mol (RDKit mol): The RDKit molecule object with upto 10 conformers to try in AMSOL.
            - netcharge (int): The net charge of the molecule.
            - rigid_scaffolds (list): List of rigid scaffold structures.
            - sulfo_matches (list): List of sulfonate group matches.
            - error (str or None): Error message if an exception occurred, otherwise None.

    """
    if VERBOSE:
        print("Generating initial 3D conformations...")
    try:
        amsol_mol, netcharge, rigid_scaffolds, sulfo_matches = embed_smiles_rdkit(smiles, name, randomSeed=randomSeed, VERBOSE=VERBOSE)
        queue.put((amsol_mol, netcharge, rigid_scaffolds, sulfo_matches, None))
    except Exception as e:
        queue.put((None, None, None, None, str(e)))

def embed_smiles_babel(smiles, name, VERBOSE=False):
    # Step 1: Generate 3D conformer in Open Babel
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("smi", "sdf")
    
    mol = ob.OBMol()
    obConversion.ReadString(mol, smiles)
    mol.AddHydrogens()
    
    builder = ob.OBBuilder()
    builder.Build(mol)  # 3D coordinate generation

    # Step 2: Minimize energy using MMFF94 force field
    ff = ob.OBForceField.FindForceField("MMFF94s")
    ff.Setup(mol)
    
    # Step 3: Apply multi-step minimization
    # Steepest Descent
    ff.SteepestDescent(250, 1.0e-4)
    #ff.FastRotorSearch(True) # permute central bonds
    ff.WeightedRotorSearch(100, 25) # 100 cycles, each with 25 forcefield ops
    # Conjugate Gradients
    ff.ConjugateGradients(250, 1.0e-4)
    #feel free to tweak these to your balance of time / quality

    #update the coordinates
    ff.GetCoordinates(mol)

    # Step 4: Convert to SDF in-memory and read into RDKit
    sdf_data = obConversion.WriteString(mol)
    mol_rdkit = Chem.MolFromMolBlock(sdf_data, removeHs=False)
    
    # Set molecule properties
    mol_rdkit.SetProp("_Name", name)
    net_charge = Chem.GetFormalCharge(mol_rdkit)

    conjugated_substituted_Ns = smi2db2_utils.find_conjugated_substituted_nitrogen1(mol_rdkit)
    additional_conjugated_substituted_Ns = smi2db2_utils.find_conjugated_substituted_nitrogen2(mol_rdkit)
    amide_linkages = smi2db2_utils.find_amide(mol_rdkit)
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_rdkit, mmffVariant="MMFF94s")
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_rdkit, mp, confId=0)
    if conjugated_substituted_Ns:
        for a, b, c, d in conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5)
        for a, b, c, d in additional_conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 178, 182, 1)
    if amide_linkages:
        for a, b, c, d, e in amide_linkages: # O=C-N(-C)-H should be coplanar.
            ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 1)
            ff.MMFFAddTorsionConstraint(a, b, c, e, False, 180, 180, 1)
    ff.Minimize()
    rigid_scaffolds = [Chem.Mol(mol_rdkit)] # Replicate the output from embed_rdkit
    return mol_rdkit, net_charge, rigid_scaffolds, list()



def stochastic_sampling(mol, tolerance_level, match_torlib, numConfs,
                        total_possible_solutions, window = 25, max_attempts=15000,
                        product=list(), unvisited = None, visited=None):
    """
    Perform stochastic sampling of the conformational space of a molecule.

    Args:
        mol: The molecule for which to generate conformers.
        tolerance_level: Specifies the tolerance level for dihedral angles.
        match_torlib: List of torsion matches (rotatable bonds).
        numConfs: Maximum number of conformers to generate.
        total_possible_solutions: The total number of possible dihedral combinations.
        window: Energy window for filtering conformers.
        max_attempts: Maximum number of attempts to generate a valid conformer.
        product: List of conformers and their energies.
        unvisited: List of unvisited conformer combinations.
        visited: Set of visited conformers.

    Returns:
        product: List of generated conformers and their energies.
        visited: Updated set of visited conformers.
        unvisited: Updated list of unvisited conformer combinations.
    """
    # Initialize the force field using MMFF94s
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")

    if product: min_energy = min([conf[1] for conf in product])
    else: min_energy = 1e6

    attempts = 0

    bonded_pairs, same_parent_pairs = precompute_bonded_and_same_parent_pairs(mol)
    # Condition to switch between visited matrix and unvisited set approaches
    if total_possible_solutions > 2 * numConfs:
        # Use visited matrix approach for large spaces
        n_transform = len(match_torlib)  # Number of rotatable bonds
        visitting = [0 for _ in range(n_transform)]
        if visited is None: visited = set()

        while len(product) < numConfs:
            for idx in range(n_transform):
                bond_idx = random.randint(0, len(match_torlib) - 1)
                bond = match_torlib[bond_idx]
                peaks = bond[2]  # Extract peaks
                peak_idx = random.choices(range(len(peaks)), weights=[peak[3] for peak in peaks], k=1)[0]
                visitting[bond_idx] = peak_idx
                peak = peaks[peak_idx]
                rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *bond[1], value=get_random_angle(peak[0], peak[tolerance_level]))

            # Check if the conformation is valid and not already visited
            if check_too_close_nonbonded_atoms(mol.GetConformer(0), mol, bonded_pairs, same_parent_pairs)\
            or tuple(visitting) in visited:
                attempts += 1
                if attempts > max_attempts:
                    break
                continue
            
            visited.add(tuple(visitting.copy()))
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mp, confId=0)
            energy = ff.CalcEnergy()
            if energy < min_energy: min_energy = energy
            if energy <= min_energy + window: product.append((Chem.Conformer(mol.GetConformer(0)), energy))


    else:
        # Use unvisited set approach for smaller spaces
        if unvisited is None:
            combination_ranges = [range(len(peaks)) for _, _, peaks in match_torlib]
            unvisited = list(itertools.product(*combination_ranges))

        while len(product) < numConfs and len(unvisited) > 0:
            choice = random.choice(unvisited)

            # Set the dihedrals based on the chosen combination
            for idx, (_, bond, _) in enumerate(match_torlib):
                peak = match_torlib[idx][2][choice[idx]]
                value = get_random_angle(peak[0], peak[tolerance_level])
                rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *bond, value=value)

            # If atoms are too close or if we already visited this conformation
            if check_too_close_nonbonded_atoms(mol.GetConformer(0), mol, bonded_pairs, same_parent_pairs):
                attempts += 1
                if attempts > max_attempts:
                    break
                continue

            unvisited.remove(choice)  # Remove the chosen combination from the unvisited set
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mp, confId=0)
            energy = ff.CalcEnergy()
            if energy < min_energy: min_energy = energy
            if energy <= min_energy + window: product.append((Chem.Conformer(mol.GetConformer(0)), energy))
            #product.append((Chem.Conformer(mol.GetConformer(0)), energy))
    return product, visited, unvisited

def conf_sampling(rigid_scaffolds, name, smiles, numConfs, sulfo_matches,
                  energywindow, ignoreTorlib=False, cleanup=True, VERBOSE=False):
    """
    A function to prepare the input and process the output from stochastic sampling function.
    
    It counts the possible conformations based on rotatable bonds, adjusts the number of requested
    conformations accordingly, and generates conformers within a specified energy window.
    The resulting conformers are aligned, filtered, and stored as Mol2 objects. Optionally,
    SDF files can be written for each conformer.
    Args:
        rigid_scaffolds (List[rdkit.Chem.Mol]): A list of RDKit molecule objects representing rigid scaffolds.
        name (str): The base name used for output files and directories.
        smiles (str): The SMILES string of the molecule to generate conformers for.
        numConfs (int): The initial number of conformations to generate.
        sulfo_matches (Any): Matches related to sulfonamide groups, determining if special handling is needed.
        energywindow (float): The maximum allowed energy difference (in appropriate units) for conformers.
        ignoreTorlib (bool, optional): If True, torsion library is ignored during conformer sampling. Defaults to False.
        cleanup (bool, optional): If True, cleans up intermediate conformers. Defaults to True.
        VERBOSE (bool, optional): If True, prints detailed processing information. Defaults to False.
    Returns:
        List[List[mol2.Mol2]]: A nested list containing Mol2 objects for each rigid scaffold, each with their respective conformers.
    
    Notes:
        For assymetric sulfonamides, there would be two versions of rigid scaffolds handled by EirVS.
    """
    num_confs_by_rotbonds, match_torlib = count_confs_by_rotbonds(rigid_scaffolds[0], ignoreTorlib, VERBOSE)
    requested_num_confs = numConfs

    mol2_obj = mol2.Mol2(mol2fileName=f'../../3d/{name}/{name}.mol2')
    mol2_obj.smiles = smiles
    mol2_obj.cleanConfs()  # Clean previous conformations if any

    #if VERBOSE: print(f"\t{num_confs_by_rotbonds} {num_confs_H} {num_rotatable_H} {numConfs}")
    if VERBOSE: print(f"\tTheory: {num_confs_by_rotbonds} possible conformations")

    # Find the rigid part only once outside the loop to save processing time
    atom_maps, label_map = find_rigid_part(rigid_scaffolds[0], rigid_rules)

    # Molecules which don't have rings are not of interest --> only sample limitedly.
    if label_map is not None: numConfs = 30
    # For very flexible molecules, we need to sample more, then filter by energy later
    else:
        if numConfs*10 < num_confs_by_rotbonds: numConfs = min(int(numConfs * 1.5), num_confs_by_rotbonds)
        elif numConfs*5 < num_confs_by_rotbonds: numConfs = min(int(numConfs * 1.25), num_confs_by_rotbonds)
        else: numConfs = min(numConfs, num_confs_by_rotbonds)


    if VERBOSE:
        rigid_info = f'\tFound {atom_maps} ({label_map}) as a rigid part' if label_map else f'\tFound {atom_maps} (rings) as rigid parts'
        print(rigid_info)

    mol2_per_rigid_scaffold = []

    if (len(match_torlib) == 0):
        clean_mol2_obj = copy.deepcopy(mol2_obj)

        for idx, mol in enumerate(rigid_scaffolds):
            clean_mol2_obj.atomXyz.append([])  # Initialize a list for atom coordinates
            for atom_idx in range(mol.GetNumAtoms()):
                pos = mol.GetConformer(0).GetAtomPosition(atom_idx)
                clean_mol2_obj.atomXyz[-1].append((float(pos.x), float(pos.y), float(pos.z)))
            clean_mol2_obj.xyzCount = 1
            mol2_per_rigid_scaffold.append([clean_mol2_obj])
            with Chem.SDWriter(f"{name}_mol{idx}.sdf") as sdwriter:
                sdwriter.write(mol, 0)
        return mol2_per_rigid_scaffold
    

    for idx, mol in enumerate(rigid_scaffolds):
        if VERBOSE: print(f"\tHandling ring/sulfonamide conformation {idx+1}/{len(rigid_scaffolds)}")
        original_mol = Chem.Mol(mol)
        # Only remap the match_torlib when sulfo_matches is found
        if sulfo_matches: num_confs_by_rotbonds, match_torlib = count_confs_by_rotbonds(mol, ignoreTorlib, VERBOSE)
        if VERBOSE: print('\tRunning stochastic torsional sampling')
        
        product, visited, unvisited = stochastic_sampling(mol, 1, match_torlib, numConfs, num_confs_by_rotbonds, energywindow, 15000, list(), visited = None, unvisited=None)

        if len(product) <= min(numConfs, num_confs_by_rotbonds) // 3: 
            if VERBOSE: print(f'Failed for stochastic sampling (generated {len(product)} confs), use the 2nd tolerance level')
            product, visited, unvisited = stochastic_sampling(mol, 2, match_torlib, numConfs, num_confs_by_rotbonds, energywindow, 15000, product, visited = visited, unvisited = unvisited)

        if len(product) == 0:
            # A backup when no conformer is found for given tolerance levels and SMILES
            clean_mol2_obj = copy.deepcopy(mol2_obj)
            clean_mol2_obj.atomXyz.append([])  # Initialize a list for atom coordinates
            for atom_idx in range(mol.GetNumAtoms()):
                pos = mol.GetConformer(0).GetAtomPosition(atom_idx)
                clean_mol2_obj.atomXyz[-1].append((float(pos.x), float(pos.y), float(pos.z)))
            clean_mol2_obj.xyzCount = 1 
            mol2_per_rigid_scaffold.append([clean_mol2_obj])
            with Chem.SDWriter(f"{name}_mol{idx}.sdf") as sdwriter:
                sdwriter.write(mol, 0)
            continue

        product.sort(key=lambda x: x[1]) #Sort by energy
        product = product[:requested_num_confs]
        before_energy = len(product)
        min_energy = product[0][1]
        result_mol = Chem.Mol(mol)
        result_mol.RemoveAllConformers()
        product = [x[0] for x in product if x[1] - min_energy <= energywindow]


        for conf in product:
            result_mol.AddConformer(conf, assignId=True)
        if VERBOSE: print(f"\tEnergy filter: {before_energy} -> {len(product)}")

        mol2_objs = []

        for atom_map in atom_maps:
            # Load the Mol2 object for the specified molecule
            clean_mol2_obj = copy.deepcopy(mol2_obj)

            # Process each conformer in the result molecule
            for conf_id in range(result_mol.GetNumConformers()):
                clean_mol2_obj.atomXyz.append([])  # Initialize a list for atom coordinates

                # Align the conformer to the original molecule using the atom map
                rdMolAlign.AlignMol(result_mol, original_mol, conf_id, 0, atomMap=[(i, i) for i in atom_map])
                conf = result_mol.GetConformer(conf_id)

                # Extract the atomic positions for the aligned conformer
                for atom_idx in range(result_mol.GetNumAtoms()):
                    pos = conf.GetAtomPosition(atom_idx)
                    clean_mol2_obj.atomXyz[-1].append((float(pos.x), float(pos.y), float(pos.z)))

            # Set the number of conformations for this Mol2 object
            clean_mol2_obj.xyzCount = result_mol.GetNumConformers()
            mol2_objs.append(clean_mol2_obj)

        if not(cleanup): 
            with Chem.SDWriter(f"{name}_mol{idx}.sdf") as sdwriter:
                for confid in range(result_mol.GetNumConformers()):
                    sdwriter.write(result_mol, confId=confid)

        mol2_per_rigid_scaffold.append(mol2_objs)

    return mol2_per_rigid_scaffold

def gen_conf_chunk(df: pd.DataFrame, args, input_file='0'):
    """
    Generate conformations for a chunk of molecules and process them through solvation, 
    3D conformation generation, and DB2 format conversion.

    This function processes a subset of molecules from a DataFrame, handles restarting from
    previously interrupted jobs, and generates DB2 files for conformations. It uses RDKit
    for initial conformer generation, AMSOL for (de)solvation, and RDKit for torsional sampling.

    Args:
        df (pd.DataFrame): 
            A DataFrame containing molecule data with at least the following columns:
            - 'ids': Unique identifiers for each molecule.
            - 'smiles': SMILES strings representing the molecules.

        args (Namespace): 
            Parsed arguments containing various configuration options, including:
            - randomSeed (int): Seed for random number generation.
            - numconfs (int): Number of conformations to generate.
            - debug (bool): Verbose output for debugging.
            - cleanup (bool): Whether to remove intermediate files after processing.
            - energywindow (float): Energy window for conformer sampling.
            - timeout (int): Timeout (in minutes) for RDKit-based conformation generation.
            - ignoretorlib (bool): Whether to ignore torsion library constraints.
            - timing (bool): If enabled, logs timing information for each step.
            - smiles (bool): If True, skips restarting logic.
            - enrichment (bool): If True, writes individual DB2 files directly instead of to a tarball.

        input_file (str): 
            The base name of the input file being processed (default is '0' for the SMILES input).

    Workflow:
        1. **Restart Handling**:
           - Checks if an output tarball (`.db2.tgz`) already exists.
           - If so, transfer previously processed molecules to a new tarball and remove the old one.

        2. **Conformer Generation**:
           - Uses RDKit to generate initial conformations with a timeout.
           - Fallback to OpenBabel if RDKit fails or times out.

        3. **Solvation with AMSOL**:
           - Processes conformers through AMSOL for solvation in water and hexane.
           - Handles failures by iterating over available conformers.

        4. **3D Conformation Sampling**:
           - Performs torsional sampling using a torsion library.

        5. **DB2 Conversion**:
           - Converts the final processed molecules to DB2 format.
           - Handles both individual file writing (enrichment mode) and tarball aggregation.

        6. **Cleanup**:
           - Removes intermediate directories (`3d`, `solv`) and temporary files if cleanup is enabled.

        7. **Timing and Logging**:
           - Logs timing information for each molecule in a CSV file if `timing` is enabled.
           - Logs errors and skips problematic molecules to ensure continuity.

    Notes:
        - This function ensures resilience by handling errors at each step and skipping problematic molecules.
        - It manages restarting from incomplete jobs to avoid redundant computation.
    """
    randomSeed, numConfs, VERBOSE, cleanup, energywindow, timeout, ignoreTorlib = args.randomSeed, \
        args.numconfs, args.debug, args.cleanup, args.energywindow, args.timeout , args.ignoretorlib
    env = setup_env()
    if args.timing: 
        if not(os.path.exists('eirvs_timing.csv')): 
            with open('eirvs_timing.csv', 'w') as f: f.write('Name,Initial embedding,AMSOL,Torsional sampling,Mol2DB2,Total\n')
        logging_time = ""
    
    # Test mode in unittest, not to produce redundant files here
    if args.test: 
        os.chdir(args.prefix)

    processed_mols = set()
    os.makedirs(f"db2", exist_ok=True)
    output_tgz = f"db2/{input_file}.db2.tgz"

    # Check if the output file already exists. A sign of unfinished job
    restart_flag = False
    if not(args.smiles) and os.path.exists(output_tgz):
        
        logger.info(f"Output file {output_tgz} already exists, restarting from the last processed molecule")
        restart_tgz = f"db2/restart_{input_file}.db2.tgz"
        shutil.copy2(output_tgz, restart_tgz)
        restart_flag = True
        

    with tarfile.open(output_tgz, mode='w:gz') as output:
        # Write previously processed DB2 files to the tarball
        if restart_flag:
            try: 
                with tarfile.open(restart_tgz, mode='r:gz') as restart_file:
                    for member in restart_file.getmembers():
                        if member.isfile() and member.name.endswith(".db2"):
                            # Extract file content and keep track of processed molecules
                            processed_mols.add(member.name.split(".db2")[0])
                            output.addfile(member, restart_file.extractfile(member))
                os.remove(restart_tgz)
            except:
                logger.error(f"Error in reading the restart file {restart_tgz}. Start from the beginning")
                os.remove(restart_tgz)
        # Process the unprocessed molecules
        for idx, row in df.iterrows():
            if args.timing: start = time.time()
            random.seed(randomSeed)
            name = row['ids']
            smiles = row['smiles']
            if name in processed_mols:
                print(f"Skipping {name} as it already exists")
                logger.info(f"Skipping {name} as it already exists")
                continue
            logger.info(f"Handling {name}")
            if VERBOSE: print(f"Handling {name}")
            
            # Embed smiles into initial conformation using RDKit with a timeout
            queue = multiprocessing.Queue()
            process = multiprocessing.Process(target=generate_conformation, args=(queue, smiles, name, randomSeed, VERBOSE))
            process.start()
            process.join(timeout=timeout*60)  # default 2 minutes timeout

            # Check if process is still alive (meaning it exceeded timeout)
            if process.is_alive():
                logger.warning(f"Timeout occurred while generating conformation for {name}, using OpenBabel.")
                process.terminate()
                process.join()
                try:
                    amsol_mol, netcharge, rigid_scaffolds, sulfo_matches = embed_smiles_babel(smiles, name, VERBOSE)
                except Exception as e:
                    logger.error(f"Error in generating initial conformation using OpenBabel for {name}, skipping it {e}")
                    log_error(smiles, name)
                    continue
                if amsol_mol is None:
                    logger.error(f"Error in generating initial conformation using OpenBabel for {name}, skipping it")
                    log_error(smiles, name)
                    continue

            # Retrieve result from queue
            elif not queue.empty():
                amsol_mol, netcharge, rigid_scaffolds, sulfo_matches, error = queue.get()
                if error:
                    logger.error(f"Error in generating initial conformation using RDKIT for {name}, skipping it: {error}")
                    log_error(smiles, name)
                    continue
                amsol_mol.SetProp("_Name", name) #By somehow this implementation loses the _Name props
                for rigid_scaffold in rigid_scaffolds: rigid_scaffold.SetProp("_Name", name)
            else:
                logger.error(f"Unknown error in generating initial conformation for {name}, skipping it.")
                log_error(smiles, name)
                continue

            if args.timing: embed_time = time.time() # Time for embedding
            
            # Solvation using AMSOL
            if VERBOSE: print("Solvating...")
            os.makedirs(f"solv/{name}", exist_ok=True)
            os.chdir(f"solv/{name}")
            for conf_id in range(amsol_mol.GetNumConformers()):
                try:
                    # Idea: try from the energy minimum conformer if AMSOL fails -> next conformer until reach the last
                    if VERBOSE: print(f"\tTrying conformer: {conf_id}")
                    error_signal = 0

                    cp = Chem.Mol(amsol_mol, confId=conf_id) #Retrieve the conf_id-th conformer of mol object
                    #mol2_block = convert(Chem.MolToMolBlock(cp), "mol", "mol2")
                    #write_to_file(mol2_block, f"{name}.mol2")
                    
                    #print(Chem.MolToSmiles(cp))
                    mol2_obj = smi2db2_utils.Mol2Writer(cp)
                    mol2_obj.write_mol2(f"{name}.mol2")

                    run_amsol.prepare(f"{name}.mol2", name, netcharge)
                    error_signal = run_amsol.run('temp.in-hex', 'temp.o-hex', env)
                    if error_signal == -1: continue
                    error_signal = run_amsol.run('temp.in-wat', 'temp.o-wat', env)
                    if error_signal == -1: continue
                    error_signal = run_amsol.process_output('temp.o-wat', 'temp.o-hex', "temp.mol2", "output")#, VERBOSE=VERBOSE)
                    if error_signal == -1: continue
                    break
                except Exception as e:
                    error_signal = -1
                    continue
            os.chdir("../..")
            if error_signal == -1 and conf_id+1 == amsol_mol.GetNumConformers(): # AMSOL failed
                logger.error(f"AMSOL failed for {name}, skipping it")
                log_error(smiles, name)
                try: # Clean up the folders if error occurs. This help to not overfill the disk
                    shutil.rmtree(f"solv/{name}", ignore_errors=True)
                except: pass
                continue
            shutil.copy(f"solv/{name}/output.mol2", f"solv/{name}/{name}_solv.mol2")
            shutil.move(f"solv/{name}/output.solv", f"solv/{name}/{name}_solv.solv")
            if args.timing: amsol_time = time.time()

            # 3D generation
            if VERBOSE: print("3D generation...")
            try:
                os.makedirs(f"3d/{name}", exist_ok=True)
                shutil.copy2(os.path.join("solv", name, f"{name}_solv.mol2"), os.path.join("3d", name, f"{name}.mol2"))
                os.chdir(f"3d/{name}")
                mol2_per_rigid_scaffold = conf_sampling(rigid_scaffolds,
                                                                name, smiles, 
                                                                numConfs, sulfo_matches, 
                                                                energywindow, ignoreTorlib, 
                                                                cleanup, VERBOSE)
            except Exception as e:
                logger.error(f"Error in torsional sampling for {name}: {e}")
                os.chdir("../..")
                log_error(smiles, name)
                try: # Clean up the folders if error occurs. This help to not overfill the disk
                    shutil.rmtree(f"solv/{name}", ignore_errors=True)
                    shutil.rmtree(f"3d/{name}", ignore_errors=True)
                except: pass
                continue
            os.chdir("../..")
            if args.timing: sampling_time = time.time()

            # Mol2DB2
            if VERBOSE: print("Converting to DB2 format...")
            try:
                os.makedirs(f"db2/{name}", exist_ok=True)
                shutil.move(os.path.join("solv", name, f"{name}_solv.solv"), os.path.join("db2", name, f"{name}.solv"))
                os.chdir(f"db2/{name}")
                db2_data_all = ""
                solv_obj = solv.Solv(f"{name}.solv")
                for mol2objs in mol2_per_rigid_scaffold:
                    for mol2obj in mol2objs:
                        if args.synthon: mol2obj.longname = row['highlights']
                        db2_data = mol2db2.mol2db2_quick_ver2(mol2obj, solv_obj)
                        db2_data_all += db2_data
                if args.enrichment: write_to_file(db2_data_all, f"../{name}.db2") #Write directly to db2 files if in enrichment mode
                else: write_to_tarball(output, db2_data_all.encode('utf-8'), name=f"{name}.db2")
                os.chdir("../..")
                smi2db2_utils.remove_folders([f"solv/{name}"])
                if cleanup:
                    smi2db2_utils.remove_folders([f"3d/{name}", f"db2/{name}"])
            except Exception as e:
                logger.error(f"Error in converting {name} to DB2 format: {e}")
                os.chdir("../..")
                try: # Clean up the folders if error occurs. This help to not overfill the disk
                    shutil.rmtree(f"solv/{name}", ignore_errors=True)
                    shutil.rmtree(f"3d/{name}", ignore_errors=True)
                except: pass
                log_error(smiles, name)
                continue
            if args.timing: 
                mol2db2_time = time.time()
                logging_time += f'{name},{embed_time-start},{amsol_time-embed_time},{sampling_time-amsol_time},{mol2db2_time-sampling_time},{mol2db2_time-start}\n'
            processed_mols.add(name)

    # Use this method to remove the tarball if it is empty. 
    # The "with open" method is better to handle unexpected error that lead to corrupted files
    # (in the enrichment mode as the DB2 files are written directly)    
    if args.enrichment and os.path.exists(f"db2/{input_file}.db2.tgz"):
        os.remove(f"db2/{input_file}.db2.tgz")
    
    if not(smi2db2_utils.is_slurm_job()):
        folders_to_remove = ['3d', 'solv'] if cleanup else ['solv']
        for folder in folders_to_remove:
            if os.path.exists(folder) and os.path.isdir(folder):
                try:
                    os.rmdir(folder)
                except: pass
    if args.timing:
        with open('eirvs_timing.csv', 'a') as f:
            f.write(logging_time)

def gen_conf_chunk_corina(df: pd.DataFrame, args, input_file='0'):
    randomSeed, numConfs, VERBOSE, cleanup, energywindow, timeout, ignoreTorlib = args.randomSeed, \
        args.numconfs, args.debug, args.cleanup, args.energywindow, args.timeout , args.ignoretorlib
    env = setup_env()
    if args.timing: 
        if not(os.path.exists('eirvs_timing.csv')): 
            with open('eirvs_timing.csv', 'w') as f: f.write('Name,Initial embedding,AMSOL,Torsional sampling,Mol2DB2,Total\n')
        logging_time = ""
    processed_mols = set()
    os.makedirs(f"db2", exist_ok=True)

    processed_mols = set()
    os.makedirs(f"db2", exist_ok=True)
    output_tgz = f"db2/{input_file}.db2.tgz"

    # Check if the output file already exists. A sign of unfinished job
    restart_flag = False
    if not(args.smiles) and os.path.exists(output_tgz):      
        logger.info(f"Output file {output_tgz} already exists, restarting from the last processed molecule")
        restart_tgz = f"db2/restart_{input_file}.db2.tgz"
        shutil.copy2(output_tgz, restart_tgz)
        restart_flag = True
        

    with tarfile.open(output_tgz, mode='w:gz') as output:
        # Write previously processed DB2 files to the tarball
        if restart_flag:
            try: 
                with tarfile.open(restart_tgz, mode='r:gz') as restart_file:
                    for member in restart_file.getmembers():
                        if member.isfile() and member.name.endswith(".db2"):
                            # Extract file content and keep track of processed molecules
                            processed_mols.add(member.name.split(".db2")[0])
                            output.addfile(member, restart_file.extractfile(member))
                os.remove(restart_tgz)
            except:
                logger.error(f"Error in reading the restart file {restart_tgz}. Start from the beginning")
                os.remove(restart_tgz)
        # Process the unprocessed molecules
        for idx, row in df.iterrows():
            if args.timing: start = time.time()
            random.seed(randomSeed)
            name = row['ids']
            smiles = row['smiles']
            if name in processed_mols:
                print(f"Skipping {name} as it already exists")
                logger.info(f"Skipping {name} as it already exists")
                continue
            logger.info(f"Handling {name}")
            if VERBOSE: print(f"Handling {name}")
            os.makedirs(f"solv/{name}", exist_ok=True)
            os.chdir(f"solv/{name}")
            try:
                amsol_mol, netcharge, rigid_scaffolds, sulfo_matches = smi2db2_utils.embed_smiles_corina(smiles, name, VERBOSE)
                
            except Exception as e:
                logger.error(f"Error in generating initial conformation using Corina for {name}, skipping it {e}")
                os.chdir("../..")
                try: # Clean up the folders if error occurs. This help to not overfill the disk
                    shutil.rmtree(f"solv/{name}", ignore_errors=True)
                except: pass
                log_error(smiles, name)     
                continue
            

            if args.timing: embed_time = time.time() # Time for embedding
            
            # Solvation using AMSOL
            if VERBOSE: print("Solvating...")
            
            for conf_id in range(amsol_mol.GetNumConformers()):
                # Idea: try from the energy minimum conformer if AMSOL fails -> next conformer until reach the last
                try:
                    if VERBOSE: print(f"\tTrying conformer: {conf_id}")
                    error_signal = 0
                    run_amsol.prepare(f"{name}.mol2", name, netcharge)
                    error_signal = run_amsol.run('temp.in-hex', 'temp.o-hex', env)
                    if error_signal == -1: continue
                    error_signal = run_amsol.run('temp.in-wat', 'temp.o-wat', env)
                    if error_signal == -1: continue
                    error_signal = run_amsol.process_output('temp.o-wat', 'temp.o-hex', "temp.mol2", "output")#, VERBOSE=VERBOSE)
                    if error_signal == -1: continue
                    break
                except Exception as e:
                    error_signal = -1
                    continue
            os.chdir("../..")

            if error_signal == -1 and conf_id+1 == amsol_mol.GetNumConformers(): # AMSOL failed
                logger.error(f"AMSOL failed for {name}, skipping it")
                try: # Clean up the folders if error occurs. This help to not overfill the disk
                    shutil.rmtree(f"solv/{name}", ignore_errors=True)
                except: pass
                log_error(smiles, name)
                continue
            shutil.copy(f"solv/{name}/output.mol2", f"solv/{name}/{name}_solv.mol2")
            shutil.move(f"solv/{name}/output.solv", f"solv/{name}/{name}_solv.solv")
            if args.timing: amsol_time = time.time()

            # 3D generation
            if VERBOSE: print("3D generation...")
            try:
                os.makedirs(f"3d/{name}", exist_ok=True)
                shutil.copy2(os.path.join("solv", name, f"{name}_solv.mol2"), os.path.join("3d", name, f"{name}.mol2"))
                os.chdir(f"3d/{name}")
                
                mol2_per_rigid_scaffold = conf_sampling(rigid_scaffolds,
                                                                name, smiles, 
                                                                numConfs, sulfo_matches, 
                                                                energywindow, ignoreTorlib, 
                                                                cleanup, VERBOSE)
            except Exception as e:
                logger.error(f"Error in torsional sampling for {name}: {e}")
                os.chdir("../..")
                try: # Clean up the folders if error occurs. This help to not overfill the disk
                    shutil.rmtree(f"solv/{name}", ignore_errors=True)
                    shutil.rmtree(f"3d/{name}", ignore_errors=True)
                except: pass
                log_error(smiles, name)
                continue
            os.chdir("../..")
            if args.timing: sampling_time = time.time()

            # Mol2DB2
            if VERBOSE: print("Converting to DB2 format...")
            os.makedirs(f"db2/{name}", exist_ok=True)
            
            try:
                shutil.move(os.path.join("solv", name, f"{name}_solv.solv"), os.path.join("db2", name, f"{name}.solv"))
                os.chdir(f"db2/{name}")
                db2_data_all = ""
                solv_obj = solv.Solv(f"{name}.solv")
                for mol2objs in mol2_per_rigid_scaffold:
                    for mol2obj in mol2objs:
                        db2_data = mol2db2.mol2db2_quick_ver2(mol2obj, solv_obj)
                        db2_data_all += db2_data
                if args.enrichment: write_to_file(db2_data_all, f"../{name}.db2") #Write directly to db2 files if in enrichment mode
                else: write_to_tarball(output, db2_data_all.encode('utf-8'), name=f"{name}.db2")
                os.chdir("../..")
                smi2db2_utils.remove_folders([f"solv/{name}"])
                if cleanup:
                    smi2db2_utils.remove_folders([f"3d/{name}", f"db2/{name}"])
            except Exception as e:
                logger.error(f"Error in converting {name} to DB2 format: {e}")
                os.chdir("../..")
                try: # Clean up the folders if error occurs. This help to not overfill the disk
                    shutil.rmtree(f"solv/{name}", ignore_errors=True)
                    shutil.rmtree(f"3d/{name}", ignore_errors=True)
                except: pass
                log_error(smiles, name)
                continue
            if args.timing: 
                mol2db2_time = time.time()
                logging_time += f'{name},{embed_time-start},{amsol_time-embed_time},{sampling_time-amsol_time},{mol2db2_time-sampling_time},{mol2db2_time-start}\n'
            processed_mols.add(name)
        
    if args.enrichment and os.path.exists(f"db2/{input_file}.db2.tgz"):
        os.remove(f"db2/{input_file}.db2.tgz")
    
    if not(smi2db2_utils.is_slurm_job()):
        folders_to_remove = ['3d', 'solv'] if cleanup else ['solv']
        for folder in folders_to_remove:
            if os.path.exists(folder) and os.path.isdir(folder):
                try:
                    os.rmdir(folder)
                except: pass
    if args.timing:
        with open('eirvs_timing.csv', 'a') as f:
            f.write(logging_time)