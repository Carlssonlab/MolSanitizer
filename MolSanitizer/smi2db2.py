"""Based partly on: 
    https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py
    http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
    http://pubs.acs.org/doi/abs/10.1021/ci2004658
    https://greglandrum.github.io/rdkit-blog/posts/2024-02-11-more-multithreading.html

    Should try to sample all possible conformations based on dihedral angles sampling based on: https://github.com/dkoes/rdkit-scripts/blob/master/rdallconf.py
"""
# Author: Thua-Phong Lam, Jens Carlsson lab, Uppsala University
import pandas as pd
import itertools
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolTransforms, rdDistGeom, rdMolAlign, rdMolDescriptors
import os, glob, shutil
import subprocess
from pathlib import Path
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
from MolSanitizer.amsol import run_amsol
from MolSanitizer.db2 import mol2db2, hydrogens, mol2
from MolSanitizer import strain_filter, smi2db2_utils
import random
import time

import logging
logger = logging.getLogger('molsani')

#rotatable_pattern=r'[*]~[*;!$(*#*)!$([!#6&X2H])!$([!#6&X3H2])]-&!@[*;!$(*#*)!$([!#6&X2H])!$([!#6&X3H2])]~[*]'
rotatable_pattern=r'[*]~[*;!$(*#*)]-&!@[*;!$(*#*)]~[*]'

Torlib = strain_filter.parse_torlib()
rigid_rule_files = Path(__file__).parent / 'Data' / 'rigid_part_rules.txt'
rigid_rules = pd.read_csv(rigid_rule_files, header=None, sep ='\s+', names=['SMARTS','label'])
rigid_rules['mol'] = rigid_rules['SMARTS'].apply(lambda x: Chem.MolFromSmarts(x))


planar_lib, non_planar_lib = strain_filter.parse_sr_confs_library()

def setup_env():
    env = os.environ.copy()
    script_dir = Path(__file__).parent
    extra_libs_path = script_dir / "libs" / "extralibs-2"

    if 'LD_LIBRARY_PATH' in env:
        env['LD_LIBRARY_PATH'] += f":{script_dir}:{extra_libs_path}"
    else:
        env['LD_LIBRARY_PATH'] = f"{script_dir}:{extra_libs_path}"
    return env

def write_to_file(content, file):
    with open(file, "w") as f:
        f.write(content)

def convert(data, inf, otf):
    obConversion = ob.OBConversion()
    obMol = ob.OBMol()
    obConversion.SetInAndOutFormats(inf, otf)
    obConversion.ReadString(obMol, data)

    return obConversion.WriteString(obMol)

def rmsd_filter(mol, ref_conf, conf_energies, threshold):
    """
    Ref:    https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.GetConformerRMS
            https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py"""
    # we use heavy atoms RMSD; not all atoms (Peter Gedeck's suggestion)
    mol_noH = Chem.Mol(mol)
    mol_noH = Chem.RemoveHs(mol_noH)
    ref_conf_id = ref_conf.GetId()
    res = []
    for e, curr_conf in conf_energies:
        curr_conf_id = curr_conf.GetId()
        rms = rdMolAlign.GetBestRMS(mol_noH, mol_noH, ref_conf_id, curr_conf_id, numThreads=1, maxMatches=10000) #Use this to avoid the symmetrical problem
        if rms > threshold:
            res.append((e, curr_conf))
    return res

def get_num_confs_for_mol(mol):
    """
    Ref:    https://pubs.acs.org/doi/full/10.1021/ci2004658"""
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True)
    if rb <= 7: return 50
    elif 8 <= rb <= 12: return 200
    else: return 300
    
def embed_smiles(smiles, name, rmsd=0.25, randomSeed=42, VERBOSE=False):
    """
    Embed SMILES into multiple conformations, minimize using MMFF94s, and return the MOL2 format.
    
    Args:
    smiles (str): The SMILES string of the molecule.
    name (str): The name of the molecule.
    rmsd (float): The RMSD threshold for pruning conformations.
    randomSeed (int): The random seed for reproducibility.
    numConfs (int): The number of conformations to generate.
    
    Returns:
    str: The lowest energy conformation in MOL2 format.
    """
    
    params = rdDistGeom.srETKDGv3()
    params.numThreads = 0  # Use all available threads
    params.pruneRmsThresh = 0.5  # Prune conformations that are too similar, not user-definable here
    params.randomSeed = randomSeed # For reproducibility
    #params.useMacrocycleTorsions = True
    #params.useSmallRingTorsions = True

    mol = Chem.MolFromSmiles(smiles)
    numConfs = 300
    #numConfs = get_num_confs_for_mol(mol)
    #if VERBOSE: print(f"\tTry with {numConfs} conformations")
    mol_H = Chem.AddHs(mol)
    res = Chem.Mol(mol_H) # res = result molecule with conformations
    res.RemoveAllConformers() # An empty conformer list

    conf_energies = []
    
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(1) #1 means vacumn, 80 means water, 20 is the compromised value (still arbitrary)

    for cid in rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId=cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        conformer = mol_H.GetConformer(cid)
        conf_energies.append((energy, conformer))
    
    conf_energies = sorted(conf_energies, key=lambda x: x[0]) # sort by increasing E

    while (res.GetNumConformers() < 10) and (len(conf_energies) > 0): # Attempts to reduce computational power!
        energy, conformer = conf_energies.pop(0) # get the lowest energy conformer
        res.AddConformer(conformer, assignId = True) # add it to the conformations list
        conf_energies = rmsd_filter(mol_H, conformer, conf_energies, rmsd) # remove all conformers that are too similar to it

    
    if VERBOSE: print(f'\tBefore: {len(mol_H.GetConformers())}, after: {res.GetNumConformers()}')
    res.SetProp("_Name", name)
    netcharge = sum(atom.GetFormalCharge() for atom in res.GetAtoms())
    return res, netcharge

def load_molecule(file_path_mol: str) -> Chem.rdchem.Mol:
    """
    Read a MOL2 molecule file and return an RDKit Mol object.
    Parameters:
    - file_path_mol: str or Path, path to the molecule file.
    Returns:
    - mol: RDKit Mol object, or None if the format is unsupported or loading fails.
    """

    # Extract the file extension
    # Convert file_path to a Path object if it's not already one
    return Chem.MolFromMol2File(file_path_mol, removeHs=False, sanitize=False, cleanupSubstructures=True)

def write_to_sdf(mol, filename):
    with Chem.SDWriter(filename+'.sdf') as writer:
        for i, conf in enumerate(mol.GetConformers()):
            rdMolAlign.AlignMol(mol, mol, prbCid = conf.GetId(), refCid = 0)
            mol.SetProp("_Name", f"Conformer_{i+1}")
            writer.write(mol, confId=conf.GetId())


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


def count_confs_by_rotbonds(mol, VERBOSE=False):
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

    # We need to put the terminal rot bonds at the beginning as they dont 
    # contribute much to the diversity of the conformations
    for rotatable_bond in rotatable_bonds:
        if is_terminal(mol, rotatable_bond[1:3]):
            if is_symmetric(mol, rotatable_bond[1:3]):
                #Exclude meaningless angles for symmetric terminal rotatable bonds (eg. CH3, CF3 only rotate onces)
                for rule in match_torlib:
                    if set(rotatable_bond[1:3]) == set(rule[1][1:3]):
                        while len(rule[2]) > 2: # rotate these angles by 120 degrees is the same
                            while (len(rule[2]) >= 2) and ((rule[2][0][0] - rule[2][1][0]) % 120 == 0): rule[2].pop(1)
                            if len(rule[2]) > 2 : rule[2].pop(-1)
                            


    num_confs = 1
    match_torlib_clean = []
    for bond in rotatable_bonds:
        for rule in match_torlib:
            if set(bond[1:3]) == set(rule[1][1:3]):
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
    while True:
        random_angle = random.gauss(mean, tolerance)
        if mean-tolerance < random_angle < mean+tolerance: 
            break

    # Normalize to the [-180, 180] range
    normalized_angle = (random_angle + 180) % 360 - 180
    
    return normalized_angle

# Deprecated Aug 22, 2024
def already_sampled(sampled_mol, current_mol, map_alignment, threshold=0.5):    
    """
    https://greglandrum.github.io/rdkit-blog/posts/2023-03-02-clustering-conformers.html
    Check if the current conformer is already sampled in the product_list.

    Parameters:
    sampled_mol (rdkit.Chem.Mol): The molecule with the sampled conformers.
    current_mol (rdkit.Chem.Mol): The molecule with the current conformer.
    map_alignment (list): The atom mapping between the original and current conformers.

    Returns:
    bool: True if the current conformer is already sampled, False otherwise.
    """
    # Remove hydrogens from the molecule
    #print(map_alignment)
    sampled_no_H = Chem.RemoveHs(sampled_mol)
    current_no_H = Chem.RemoveHs(current_mol)
    # Iterate through each conformer in the product list
    for conf_id in range(sampled_no_H.GetNumConformers()):
        # Compute RMSD between current conformer and each conformer in the product_list
        rmsd = rdMolAlign.GetBestRMS(current_no_H, sampled_no_H, 0, conf_id, maxMatches=1000)
        if rmsd <= threshold:
            return True
    return False

def torsional_scan_rand(mol, conf, i, matches, match_torlib, sdwriter, product, original,
                        numConfs, atom_maps, visited, visitting):
    '''
    Recursively enumerates all angles for matching dihedrals.
    
    Parameters:
    - mol: The molecule object.
    - conf: The conformer object.
    - i: The index of the dihedral being enumerated.
    - matches: The list of matching dihedrals.
    - match_torlib: The torsion library for matching dihedrals.
    - sdwriter: The SD writer object.
    - product: The list of conformers.
    - original: The original molecule object (for alignment).
    - numConfs: The maximum number of conformers to generate.
    - atom_maps: The list of atom maps.
    - visited: The set of visited dihedrals.
    
    Returns:
    - product: The list of conformers.
    '''

    if len(product) >= numConfs: return product, visited
    if i >= len(matches): #base case, torsions should be set in conf
        #print(check_too_close_nonbonded_atoms(mol.GetConformer(conf), mol))
        if check_too_close_nonbonded_atoms(mol.GetConformer(conf), mol): return product, visited
        product.append(Chem.Conformer(mol.GetConformer(conf))) 
        visited.add(tuple(visitting.copy()))
        rdMolAlign.AlignMol(mol, original, conf, 0, atomMap=[(i, i) for i in atom_maps])
        sdwriter.write(mol, conf)
        return product, visited
    else:
        peaks = strain_filter.extract_peaks(match_torlib, matches[i][1:3])
        dihedral_4_atoms = peaks[0]
        for peakidx, (prefered, tolerance, _ , _) in enumerate(peaks[1]):
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(conf),*dihedral_4_atoms,value = get_random_angle(prefered, tolerance))
            visitting[i] = peakidx
            product, visited = torsional_scan_rand(mol, conf, i+1, matches, match_torlib, sdwriter, 
                                                   product, original, numConfs, atom_maps, visited, visitting)        
        return product, visited
    

def stochastic_sampling(mol, tolerance_level, reordered_rot_bonds, match_torlib, sdwriter, original, numConfs, rmsd, atom_maps, max_attempts = 100, product = [], visited = set()):
    """
    Perform stochastic sampling of the conformational space of a molecule using a Monte Carlo method.

    This function generates a specified number of conformers for a given molecule by iteratively modifying
    its torsional angles according to predefined torsion library rules. The method utilizes a Monte Carlo 
    approach to explore the conformational space, checking for already sampled or invalid conformers 
    (e.g., those with non-bonded atoms too close to each other) and aligning the resulting conformers to 
    a reference structure.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object to sample.
        tolerance_level (int): The tolerance level for sampling torsional angles:
                               1 for relaxed sampling, 2 for more tolerable sampling.
        reordered_rot_bonds (list): A list of rotatable bonds in the molecule, reordered for processing.
        match_torlib (list): The torsion library rules that dictate the allowed torsional angles and their probabilities.
        sdwriter (rdkit.Chem.SDWriter): An SDWriter object to output the generated conformers to an SD file.
        original (rdkit.Chem.Mol): The original conformation of the molecule used for alignment reference.
        numConfs (int): The number of unique conformers to generate.
        rmsd (float): The RMSD threshold for pruning conformers.
        atom_maps (list of tuples): A list of tuples representing the atom mapping between the original
                                    conformation and the generated conformers, used for rigid part alignment.
        max_attempts (int, optional): The maximum number of attempts to find a valid, unique conformer before stopping. 
                                      Defaults to 100.
        product (list, optional): A list to store the successfully generated conformers. Defaults to an empty list.
        visited (set, optional): A set to store the visited dihedral angles (especially from the torsional scan) to avoid redundant conformers. 
                                 Defaults to an empty

    Returns:
        list: A list of RDKit Conformer objects representing the successfully sampled conformers.
    """
    sampled_mol = Chem.Mol(mol)
    sampled_mol.RemoveAllConformers()
    for conf in product:
        sampled_mol.AddConformer(conf, assignId = True)

    max_angles = 0 
    for rule in match_torlib:
        max_angles = max(len(rule[2]), max_angles)
    #visit_matrix = np.zeros((len(reordered_rot_bonds), max_angles))
    n_transform = len(reordered_rot_bonds)
    visitting = [0 for _ in range(n_transform)]
    attempts = 0
    while (len(product) < numConfs):
        for idx in range(n_transform):
            # Each rotatable bond has equally likely chance to be selected
            bond_idx = random.randint(0, len(reordered_rot_bonds)-1)
            bond = reordered_rot_bonds[bond_idx]
            
            # Which peak to be selected is based on the weights of the peaks (defined by the "score" in TorLib)
            peaks = strain_filter.extract_peaks(match_torlib, bond[1:3])
            #peak_idx = random.randint(0, len(peaks[1])-1)
            peak_idx = random.choices(range(len(peaks[1])), weights = [peak[3] for peak in peaks[1]], k=1)[0]
            #print(peak_idx)
            visitting[bond_idx] = peak_idx
            peak = peaks[1][peak_idx]
            #visit_matrix[bond_idx][peak_idx] += 1
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(0),*peaks[0],value = get_random_angle(peak[0], peak[tolerance_level]))

        if check_too_close_nonbonded_atoms(mol.GetConformer(0), mol) or \
        tuple(visitting) in visited:
        #already_sampled(sampled_mol, mol, rmsd): 
            attempts += 1
            if attempts > max_attempts: break
            continue
        attempts = 0
        rdMolAlign.AlignMol(mol, original, 0, 0, atomMap=[(i, i) for i in atom_maps])
        #print(visitting)
        visited.add(tuple(visitting.copy()))
        product.append(Chem.Conformer(mol.GetConformer(0)))
        sampled_mol.AddConformer(mol.GetConformer(0), assignId = True)
        sdwriter.write(mol,confId=0)
    #print(visit_matrix)
    return product, visited


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
    

def choose_sampling_method(mol, name, numConfs, rmsd, VERBOSE=False):
   
    # Count number of rotatable hydrogens and number of conformations contributed by them
    mol2_countH = mol2.Mol2(mol2fileName=f"{name}.mol2", nameFileName=None, mol2text=None)
    num_confs_H = hydrogens.count_confs_by_H(mol2_countH)
    num_rotatable_H = mol2_countH.hydrogensToRotate

    # Divide the number of conformations by that contributed by rotatable hydrogens
    # This adopts the same strategy from previous DB2 pipeline from UCSF
    if num_rotatable_H >= 6:  
        logger.warning(f"{name} has too many rotatable hydrogens, will reduce by 60")
        numConfs = numConfs // 60
    elif num_rotatable_H >= 4: numConfs = numConfs // 30
    elif num_rotatable_H >= 2: numConfs = numConfs // 3  
    num_confs_by_rotbonds, reordered_rot_bonds, match_torlib = count_confs_by_rotbonds(mol, VERBOSE)
    
    if VERBOSE: print(f"\t{num_confs_by_rotbonds} {num_confs_H} {num_rotatable_H} {numConfs}")

    # Create an SD writer object to output the conformers
    rotated_file = f"{name}_rotated.sdf"
    sdwriter = Chem.SDWriter(rotated_file)
    original_mol = Chem.Mol(mol)

    # If no rotatable bonds are found, use the conformations generated by RDKit
    # "." in name is a flag to indicate that the SMILES has been stereoisomerically expanded --> no need to use RDKit conformations
    if (len(reordered_rot_bonds) == 0 or num_confs_by_rotbonds == 1) and ('.' not in name):
        #atom_maps = find_rigid_part(mol, one_conf_rules)
        if VERBOSE: print(f'\tFound {atom_maps} as rigid parts')
        if VERBOSE: logger.warning(f"No rotatable bonds found for {name}, use the conformations from Rdkit")
        for conf_id in range(mol.GetNumConformers()):
            rdMolAlign.AlignMol(mol, mol, conf_id, 0, atomMap=[(i, i) for i in atom_maps])
            sdwriter.write(mol, confId=conf_id)
        sdwriter.close()
        convert_sdf_mol2(rotated_file, f"{name}_rotated.mol2", VERBOSE)
        return 2
    
    # Find rigid parts as anchor points for the molecules
    atom_maps = find_rigid_part(mol, rigid_rules)
    if VERBOSE: print(f'\tFound {atom_maps} as rigid parts')
    
    # Initialize the visited and visitting list for the systematic scan
    visitting = [-1 for _ in range(len(reordered_rot_bonds))]

    if num_confs_by_rotbonds <= numConfs:
        if VERBOSE: print('Running systematic torsional scan')
        product, visited = torsional_scan_rand(mol, conf=0, i=0, matches = reordered_rot_bonds, match_torlib = match_torlib,
                                      sdwriter=sdwriter, product=list(), original=original_mol, numConfs = numConfs, 
                                      atom_maps = atom_maps, visited = set(), visitting = visitting)
        #   If the systematic scan is not enough, do stochastic sampling
        #   This part is to prevent the case when only small torsional rotation could prevent the clashes
        #   Only produce 1/3 of the desired conformations is an indicator of clashes
        if len(product) <= num_confs_by_rotbonds // 3: 
            if VERBOSE: print(f'Failed for systematic scan (generated {len(product)} confs), use stochastic method instead')
            #second arg = 1 is using the 1st tolerance level (relaxed)
            product, visited = stochastic_sampling(mol, 1, reordered_rot_bonds, match_torlib, sdwriter, original_mol,
                                           num_confs_by_rotbonds, rmsd, atom_maps, 250, product, visited) 
            if len(product) <= num_confs_by_rotbonds // 3:
                if VERBOSE: print(f'Failed even for stochastic scan (generated {len(product)} confs), use the 2nd tolerance level')
                product, visited = stochastic_sampling(mol, 2, reordered_rot_bonds, match_torlib, sdwriter, original_mol, 
                                              num_confs_by_rotbonds, rmsd, atom_maps, 500, product, visited)
    else:
        if VERBOSE: print('Running stochastic torsional sampling')
        product, visited = stochastic_sampling(mol, 1, reordered_rot_bonds, match_torlib, sdwriter, original_mol, numConfs,
                                       rmsd, atom_maps, 250, list(), set())
    sdwriter.close()
    convert_sdf_mol2(rotated_file, f"{name}_rotated.mol2", VERBOSE)
    return 0

def log_error(smiles, name):
    with open('msani_error.log', 'a') as f:
        f.write(f"{smiles} \t {name}\n")

def gen_conf_chunk(df: pd.DataFrame, randomSeed = 42, numConfs = 10000, rmsd = 0.25, VERBOSE = False, cleanup=False):
        
        env = setup_env()
        #if VERBOSE: print(df)
        #if VERBOSE: print(env['LD_LIBRARY_PATH'])
        for idx, row in df.iterrows():
            random.seed(randomSeed)
            name = row['ids']
            # Embed smiles into initial conformation 
            # (The number of initial confs will be estimated from https://pubs.acs.org/doi/abs/10.1021/ci2004658
            # then we only use the minimal energy ones)
            if VERBOSE: print(f"\nHandling {name} \nGenerating initial 3D conformations...")
            mol, netcharge = embed_smiles(row['smiles'], name, rmsd = rmsd, 
                                          randomSeed = randomSeed, VERBOSE=VERBOSE)

            # Solvation using AMSOL
            if VERBOSE: print("Solvating...")
            subprocess.run(f"mkdir -p solv/{name}", shell=True)
            os.chdir(f"solv/{name}")

            for conf_id in range(mol.GetNumConformers()):
                # Idea: try from the energy minimum conformer if AMSOL fails -> next conformer until reach the last
                if VERBOSE: print(f"\tTrying conformer: {conf_id}")
                error_signal = 0

                cp = Chem.Mol(mol, confId=conf_id) #Retrieve the conf_id-th conformer of mol object
                mol2_block = convert(Chem.MolToMolBlock(cp), "mol", "mol2")
                write_to_file(mol2_block, f"{name}.mol2")

                run_amsol.prepare(f"{name}.mol2", name, netcharge)
                error_signal = run_amsol.run('temp.in-hex', 'temp.o-hex', env)
                if error_signal == -1: continue
                error_signal = run_amsol.run('temp.in-wat', 'temp.o-wat', env)
                if error_signal == -1: continue
                error_signal = run_amsol.process_output('temp.o-wat', 'temp.o-hex', "temp.mol2", "output")#, VERBOSE=VERBOSE)
                if error_signal == -1: continue
                break
            os.chdir("../..")
            if error_signal == -1 and conf_id+1 == mol.GetNumConformers(): # AMSOL failed
                logger.error(f"AMSOL failed for {name}, skipping it")
                log_error(row['smiles'], name)
                continue
            subprocess.run(f"cp solv/{name}/output.mol2 solv/{name}/{name}_solv.mol2", shell=True)
            subprocess.run(f"mv solv/{name}/output.solv solv/{name}/{name}_solv.solv", shell=True)
            

            # 3D generation
            if VERBOSE: print("3D generation...")
            subprocess.run(f"mkdir -p 3d/{name}", shell=True)
            subprocess.run(f"cp solv/{name}/{name}_solv.mol2 3d/{name}/{name}.mol2", shell=True)
            os.chdir(f"3d/{name}")

            sampling_signal = choose_sampling_method(mol, name, numConfs, rmsd, VERBOSE)
            os.chdir("../..")

            # Mol2DB2
            if VERBOSE: print("Converting to DB2 format...")
            subprocess.run(f"mkdir -p db2/{name}", shell=True)
            subprocess.run(f"mv 3d/{name}/{name}_rotated.mol2 db2/{name}/{name}.mol2", shell=True)
            subprocess.run(f"mv solv/{name}/{name}_solv.solv db2/{name}/{name}.solv", shell=True)    
            os.chdir(f"db2/{name}")
            if not os.path.isfile(f"{name}.mol2") or not os.path.isfile(f"{name}.solv"):
                logger.error(f"Not found mol2 files and solv for db2 generation of {name}")
                os.chdir("../..")
                log_error(row['smiles'], name)
                continue
            else:
                try:
                    if sampling_signal == 0:
                        db2_data = mol2db2.mol2db2_quick(f"{name}.mol2", f"{name}.solv")
                    elif sampling_signal == 2: 
                        # As we use different ring conformer from rdkit, we need a wider tolerance for small variations
                        # in ring-atoms (in previous case, we use exactly the same ring conformer, but not in this case)
                        db2_data = mol2db2.mol2db2_quick(f"{name}.mol2", f"{name}.solv", disttol=0.05)

                    write_to_file(db2_data, f"../{name}.db2")
                    os.chdir("../..")
                    if cleanup:
                        subprocess.run(f"rm -rf 3d/{name} solv/{name} db2/{name}", shell=True)
                except Exception as e:
                    logger.error(f"Error in converting {name} to DB2 format: {e}")
                    os.chdir("../..")
                    log_error(row['smiles'], name)
                    continue
        if cleanup:
            subprocess.run("find 3d -type d -empty -delete", shell=True)
            subprocess.run("find solv -type d -empty -delete", shell=True)

def embed_smiles_ver2(smiles, name, randomSeed=42, VERBOSE=False):
    mol_H = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mol_H.SetProp("_Name", name)
    amsol_mol = Chem.Mol(mol_H)
    empty_mol = Chem.Mol(mol_H)

    params = rdDistGeom.srETKDGv3()
    params.numThreads = 0  # Use all available threads
    params.pruneRmsThresh = 0.35  # Prune conformations that are too similar, not user-definable here
    params.randomSeed = randomSeed # For reproducibility
    params.useRandomCoords = True
    num_ring_confs = 1

    ssr = [set(ring) for ring in Chem.GetSymmSSSR(mol_H)]
    planar_rings, non_planar_rings = smi2db2_utils.get_flexible_ring(mol_H, ssr, planar_lib, non_planar_lib)
    sulfo_matches = smi2db2_utils.find_sulfonamide_like_scaffolds(mol_H)
    flippable_Ns = smi2db2_utils.find_flipped_nitrogen(mol_H)

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

    if sulfo_matches or non_planar_rings or flippable_Ns: numConfs = 100
    else: numConfs = 10


    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(1) #1 means vacumn, 80 means water, 20 is the compromised value (still arbitrary)
    conf_ring_descriptors_df = pd.DataFrame()
    try:
        for cid in rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId=cid)
            ff.Minimize()
            conformer = mol_H.GetConformer(cid)
            energy = ff.CalcEnergy()
            conf_ring_descriptors_df = smi2db2_utils.classify_confs(conformer, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df)
    except: pass

    if len(conf_ring_descriptors_df) == 0:
        # In case where srETKDGv3 failed in embedding the molecule, 
        # we have to use the macrocyclic version.
        params = rdDistGeom.ETKDGv3()
        params.numThreads = 0  # Use all available threads
        params.pruneRmsThresh = 0.35  # Prune conformations that are too similar, not user-definable here
        params.randomSeed = randomSeed # For reproducibility
        params.useRandomCoords = True
        for cid in rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId=cid)
            ff.Minimize()
            conformer = mol_H.GetConformer(cid)
            energy = ff.CalcEnergy()
            conf_ring_descriptors_df = smi2db2_utils.classify_confs(conformer, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df)


    conf_ring_descriptors_df.sort_values('Energy', inplace=True)
    # Keep a reservoir as the lowest energy possible conformer in case no confor
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
        # Use the lowest energy conformer
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

def stochastic_sampling_v3(mol, tolerance_level, match_torlib, sdwriters, original, numConfs, total_possible_solutions, atom_maps, max_attempts=100, product=0, unvisited = None, visited=None):
    """
    Perform stochastic sampling of the conformational space of a molecule using a hybrid approach.
    It switches between a visited matrix approach for large solution spaces and an unvisited set approach for smaller spaces.

    This function generates conformers of a molecule and prunes conformers that are either too close in terms of atomic distances or have already been visited (i.e., previously generated conformers). The method of tracking unvisited or visited conformations depends on the size of the search space relative to the number of allowed conformers.

    Args:
        mol (rdkit.Chem.Mol): The molecule for which to generate conformers.
        tolerance_level (int): The index in the torsion library specifying the tolerance for angle deviations.
        match_torlib (list): List of tuples containing torsion matching information for dihedral angles.
        sdwriter (rdkit.Chem.SDWriter): Writer object to output generated conformers to an SDF file.
        original (rdkit.Chem.Mol): The reference molecule to align the generated conformers to.
        numConfs (int): The maximum number of conformers to generate.
        total_possible_solutions (int): The total number of possible dihedral combinations (solution space).
        atom_maps (list): List of atom indices to use for alignment when generating conformers.
        max_attempts (int, optional): Maximum number of attempts to generate a valid conformer. Defaults to 100.
        product (int, optional): Number of conformers generated so far. Defaults to 0.
        unvisited (list or None, optional): List of unvisited conformer combinations. Defaults to None. Used when the solution space is smaller than twice the maximum allowed conformers.
        visited (set or None, optional): Set of visited conformers. Defaults to None. Used when the solution space is larger than twice the maximum allowed conformers.

    Returns:
        tuple: (product, visited, unvisited) - 
            - `product`: Number of generated conformers.
            - `visited`: Set of visited conformations if the large-space approach is used.
            - `unvisited`: List of unvisited conformer combinations if the small-space approach is used.
    
    Logic:
    - If the total number of possible conformers exceeds twice the maximum allowed conformers (`total_possible_solutions > 2 * numConfs`), a **visited matrix approach** is used. This method randomly samples from the search space and tracks visited conformers to avoid generating duplicates.
    - If the total number of possible conformers is smaller or equal to twice the maximum allowed conformers, an **unvisited set approach** is used. This method randomly selects from a set of all unvisited combinations and generates conformers until the desired number is reached or all possibilities are exhausted.

    """
    
    bonded_pairs, same_parent_pairs = precompute_bonded_and_same_parent_pairs(mol)
    # Condition to switch between visited matrix and unvisited set approaches
    if total_possible_solutions > 2 * numConfs:
        # Use visited matrix approach for large spaces
        n_transform = len(match_torlib)  # Number of rotatable bonds
        visitting = [0 for _ in range(n_transform)]
        attempts = 0
        if visited is None: visited = set()

        while product < numConfs:
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
            
            attempts = 0
            visited.add(tuple(visitting.copy()))
            product+=1
            for sdwriter, sdf_file, mol2_file, atom_map in sdwriters:
                rdMolAlign.AlignMol(mol, original, 0, 0, atomMap=atom_map)
                sdwriter.write(mol, confId=0)

    else:
        # Use unvisited set approach for smaller spaces
        if unvisited is None:
            combination_ranges = [range(len(peaks)) for _, _, peaks in match_torlib]
            unvisited = list(itertools.product(*combination_ranges))

        attempts = 0
        while product < numConfs and len(unvisited) > 0:
            # Select a random combination from the unvisited set
            # As the number of possible solutions <= allowance, we don't need to prioritize any combination
            #print(unvisited)
            #print(unvisited.dtype)
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
            attempts = 0
            product+=1                        
            for sdwriter, sdf_file, mol2_file, atom_map in sdwriters:
                rdMolAlign.AlignMol(mol, original, 0, 0, atomMap=atom_map)
                sdwriter.write(mol, confId=0)
    return product, visited, unvisited

def choose_sampling_method_ver2(rigid_scaffolds, name, numConfs, sulfo_matches, VERBOSE=False):
    """
    
    """
    num_confs_by_rotbonds, match_torlib = count_confs_by_rotbonds(rigid_scaffolds[0], VERBOSE)
    
    #if VERBOSE: print(f"\t{num_confs_by_rotbonds} {num_confs_H} {num_rotatable_H} {numConfs}")
    if VERBOSE: print(f"\tTheory: {num_confs_by_rotbonds} possible conformations")

    # Find the rigid part only once outside the loop to save processing time
    atom_maps, label_map = find_rigid_part(rigid_scaffolds[0], rigid_rules)
    
    if VERBOSE:
        rigid_info = f'\tFound {atom_maps} ({label_map}) as a rigid part' if label_map else f'\tFound {atom_maps} (rings) as rigid parts'
        print(rigid_info)

    if (len(match_torlib) == 0 or num_confs_by_rotbonds == 1):
        for idx, mol in enumerate(rigid_scaffolds):
            sdf_file, mol2_file = smi2db2_utils.get_sdf_mol2_filename(name, idx, 0)
            with Chem.SDWriter(sdf_file) as sdwriter: sdwriter.write(mol)
            convert_sdf_mol2(sdf_file, mol2_file, VERBOSE)
        return 1
    

    for idx, mol in enumerate(rigid_scaffolds):
        if VERBOSE: print(f"\tHandling ring/sulfonamide conformation {idx+1}/{len(rigid_scaffolds)}")
        original_mol = Chem.Mol(mol)
        # Only remap the match_torlib when sulfo_matches is found
        if sulfo_matches: num_confs_by_rotbonds, match_torlib = count_confs_by_rotbonds(mol, VERBOSE)
        sdwriters = []
        for align_copy, atom_map in enumerate(atom_maps):
            sdf_file, mol2_file = smi2db2_utils.get_sdf_mol2_filename(name, idx, align_copy)
            sdwriter = Chem.SDWriter(sdf_file)
            sdwriters.append((sdwriter, sdf_file, mol2_file, [(i,i) for i in atom_map]))

        if VERBOSE: print('\tRunning stochastic torsional sampling')
        
        product, visited, unvisited = stochastic_sampling_v3(mol, 1, match_torlib, sdwriters, original_mol, 
                                      numConfs, num_confs_by_rotbonds, atom_maps, 250, 0, visited = None, unvisited=None)

        if product <= min(numConfs, num_confs_by_rotbonds) // 3: 
            if VERBOSE: print(f'Failed for stochastic scan (generated {len(product)} confs), use the 2nd tolerance level')
            product, visited, unvisited = stochastic_sampling_v3(mol, 2, match_torlib, sdwriters, original_mol, 
                                          numConfs, num_confs_by_rotbonds, atom_maps, 500, product, visited = visited, unvisited = unvisited)
    
        for sdwriter, sdf_file, mol2_file, atom_map in sdwriters:
            if product == 0: sdwriter.write(mol)
            sdwriter.close()
            convert_sdf_mol2(sdf_file, mol2_file, VERBOSE)
            
    return len(sdwriters)

def stochastic_sampling_v4(mol, tolerance_level, match_torlib, numConfs, total_possible_solutions, max_attempts=100, product=list(), unvisited = None, visited=None):
    """
    Perform stochastic sampling of the conformational space of a molecule.

    Args:
        mol: The molecule for which to generate conformers.
        tolerance_level: Specifies the tolerance level for dihedral angles.
        match_torlib: List of torsion matches (rotatable bonds).
        numConfs: Maximum number of conformers to generate.
        total_possible_solutions: The total number of possible dihedral combinations.
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

    bonded_pairs, same_parent_pairs = precompute_bonded_and_same_parent_pairs(mol)
    # Condition to switch between visited matrix and unvisited set approaches
    if total_possible_solutions > 2 * numConfs:
        # Use visited matrix approach for large spaces
        n_transform = len(match_torlib)  # Number of rotatable bonds
        visitting = [0 for _ in range(n_transform)]
        attempts = 0
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
            
            attempts = 0
            visited.add(tuple(visitting.copy()))
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mp, confId=0)
            product.append((Chem.Conformer(mol.GetConformer(0)), ff.CalcEnergy()))

    else:
        # Use unvisited set approach for smaller spaces
        if unvisited is None:
            combination_ranges = [range(len(peaks)) for _, _, peaks in match_torlib]
            unvisited = list(itertools.product(*combination_ranges))

        attempts = 0
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
            attempts = 0
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mp, confId=0)
            product.append((Chem.Conformer(mol.GetConformer(0)), ff.CalcEnergy()))
    return product, visited, unvisited

def choose_sampling_method_ver3(rigid_scaffolds, name, numConfs, sulfo_matches, energywindow, VERBOSE=False):
    """
    
    """
    num_confs_by_rotbonds, match_torlib = count_confs_by_rotbonds(rigid_scaffolds[0], VERBOSE)
    
    #if VERBOSE: print(f"\t{num_confs_by_rotbonds} {num_confs_H} {num_rotatable_H} {numConfs}")
    if VERBOSE: print(f"\tTheory: {num_confs_by_rotbonds} possible conformations")

    # Find the rigid part only once outside the loop to save processing time
    atom_maps, label_map = find_rigid_part(rigid_scaffolds[0], rigid_rules)

    # Molecules which don't have rings are not of interest --> only sample limitedly.
    if label_map is not None: numConfs = 30

    if VERBOSE:
        rigid_info = f'\tFound {atom_maps} ({label_map}) as a rigid part' if label_map else f'\tFound {atom_maps} (rings) as rigid parts'
        print(rigid_info)

    if (len(match_torlib) == 0 or num_confs_by_rotbonds == 1):
        for idx, mol in enumerate(rigid_scaffolds):
            sdf_file, mol2_file = smi2db2_utils.get_sdf_mol2_filename(name, idx, 0)
            with Chem.SDWriter(sdf_file) as sdwriter: sdwriter.write(mol)
            convert_sdf_mol2(sdf_file, mol2_file, VERBOSE)
        return 1
    

    for idx, mol in enumerate(rigid_scaffolds):
        if VERBOSE: print(f"\tHandling ring/sulfonamide conformation {idx+1}/{len(rigid_scaffolds)}")
        original_mol = Chem.Mol(mol)
        # Only remap the match_torlib when sulfo_matches is found
        if sulfo_matches: num_confs_by_rotbonds, match_torlib = count_confs_by_rotbonds(mol, VERBOSE)
 
        if VERBOSE: print('\tRunning stochastic torsional sampling')
        
        product, visited, unvisited = stochastic_sampling_v4(mol, 1, match_torlib, numConfs, num_confs_by_rotbonds, 250, list(), visited = None, unvisited=None)

        if len(product) <= min(numConfs, num_confs_by_rotbonds) // 3: 
            if VERBOSE: print(f'Failed for stochastic scan (generated {len(product)} confs), use the 2nd tolerance level')
            product, visited, unvisited = stochastic_sampling_v4(mol, 2, match_torlib, numConfs, num_confs_by_rotbonds, 500, product, visited = visited, unvisited = unvisited)

        if len(product) == 0:
            # A backup when no conformer is found for given tolerance levels and SMILES
            for align_copy, atom_map in enumerate(atom_maps):
                sdf_file, mol2_file = smi2db2_utils.get_sdf_mol2_filename(name, idx, align_copy)
                with Chem.SDWriter(sdf_file) as sdwriter: sdwriter.write(mol)
                convert_sdf_mol2(sdf_file, mol2_file, VERBOSE)
            continue

        product.sort(key=lambda x: x[1]) #Sort by energy
        before_energy = len(product)
        min_energy = product[0][1]
        result_mol = Chem.Mol(mol)
        result_mol.RemoveAllConformers()
        product = [x[0] for x in product if x[1] - min_energy <= energywindow]

        #TODO: RMSD clustering if bad

        for conf in product:
            result_mol.AddConformer(conf, assignId=True)
        if VERBOSE: print(f"\tEnergy filter: {before_energy} -> {len(product)}")

        for align_copy, atom_map in enumerate(atom_maps):
            sdf_file, mol2_file = smi2db2_utils.get_sdf_mol2_filename(name, idx, align_copy)
            with Chem.SDWriter(sdf_file) as sdwriter:
                for conf_id in range(result_mol.GetNumConformers()):
                    rdMolAlign.AlignMol(result_mol, original_mol, conf_id, 0, atomMap=[(i,i) for i in atom_map])
                    sdwriter.write(result_mol, confId=conf_id)
            convert_sdf_mol2(sdf_file, mol2_file, VERBOSE)
            
    return len(atom_maps)


def gen_conf_chunk_ver2(df: pd.DataFrame, args):
    randomSeed, numConfs, VERBOSE, cleanup, energywindow = args.randomSeed, args.numconfs, args.debug, args.cleanup, args.energywindow 
    env = setup_env()
    if args.timing: 
        if not(os.path.exists('msani_timing.csv')): 
            with open('msani_timing.csv', 'w') as f: f.write('Name, Initial embedding, AMSOL, Torsional sampling, Mol2DB2, Total\n')
        logging_time = ""
    for idx, row in df.iterrows():
        if args.timing: start = time.time()
        random.seed(randomSeed)
        name = row['ids']
        if os.path.exists(f"db2/{name}/{name}.db2"):
            print(f"Skipping {name} as it already exists")
            continue
        print(f"Handling {name}")

        # Embed smiles into initial conformation using RDKit        
        if VERBOSE: print("Generating initial 3D conformations...")
        try:
            amsol_mol, netcharge, rigid_scaffolds, sulfo_matches = embed_smiles_ver2(row['smiles'], name, 
                                        randomSeed = randomSeed, VERBOSE=VERBOSE)
        except Exception as e:
                logger.error(f"Error in generating initial conformation for {name}, skipping it {e}")
                log_error(row['smiles'], name)
                continue
        if args.timing: embed_time = time.time()

        # Solvation using AMSOL
        if VERBOSE: print("Solvating...")
        os.makedirs(f"solv/{name}", exist_ok=True)
        os.chdir(f"solv/{name}")
        
        for conf_id in range(amsol_mol.GetNumConformers()):
            # Idea: try from the energy minimum conformer if AMSOL fails -> next conformer until reach the last
            if VERBOSE: print(f"\tTrying conformer: {conf_id}")
            error_signal = 0

            cp = Chem.Mol(amsol_mol, confId=conf_id) #Retrieve the conf_id-th conformer of mol object
            mol2_block = convert(Chem.MolToMolBlock(cp), "mol", "mol2")
            write_to_file(mol2_block, f"{name}.mol2")

            run_amsol.prepare(f"{name}.mol2", name, netcharge)
            error_signal = run_amsol.run('temp.in-hex', 'temp.o-hex', env)
            if error_signal == -1: continue
            error_signal = run_amsol.run('temp.in-wat', 'temp.o-wat', env)
            if error_signal == -1: continue
            error_signal = run_amsol.process_output('temp.o-wat', 'temp.o-hex', "temp.mol2", "output")#, VERBOSE=VERBOSE)
            if error_signal == -1: continue
            break
        os.chdir("../..")
        if error_signal == -1 and conf_id+1 == amsol_mol.GetNumConformers(): # AMSOL failed
            logger.error(f"AMSOL failed for {name}, skipping it")
            log_error(row['smiles'], name)
            continue
        subprocess.run(f"cp solv/{name}/output.mol2 solv/{name}/{name}_solv.mol2", shell=True)
        subprocess.run(f"mv solv/{name}/output.solv solv/{name}/{name}_solv.solv", shell=True)
        if args.timing: amsol_time = time.time()

        # 3D generation
        if VERBOSE: print("3D generation...")
        os.makedirs(f"3d/{name}", exist_ok=True)
        shutil.copy2(os.path.join("solv", name, f"{name}_solv.mol2"), os.path.join("3d", name, f"{name}.mol2"))
        os.chdir(f"3d/{name}")
        numPossibleRigidAlignments = choose_sampling_method_ver3(rigid_scaffolds, name, numConfs, sulfo_matches, energywindow, VERBOSE)
        os.chdir("../..")
        if args.timing: sampling_time = time.time()

        # Mol2DB2
        if VERBOSE: print("Converting to DB2 format...")
        os.makedirs(f"db2/{name}", exist_ok=True)
        shutil.move(os.path.join("solv", name, f"{name}_solv.solv"), os.path.join("db2", name, f"{name}.solv"))
        smi2db2_utils.move_and_rename_mol2_files(name, len(rigid_scaffolds), numPossibleRigidAlignments, VERBOSE)
        os.chdir(f"db2/{name}")
        if not any(glob.glob(f"{name}_mol*.mol2")) or not os.path.isfile(f"{name}.solv"):
            logger.error(f"Not found any mol2 files or solv for db2 generation of {name}")
            os.chdir("../..")
            log_error(row['smiles'], name)
            continue
        else:
            try:
                db2_data_all = ""
                for idx in range(len(rigid_scaffolds)):
                    for rigid_alignment in range(numPossibleRigidAlignments):
                        db2_data = mol2db2.mol2db2_quick(f"{name}_mol{idx}_align{rigid_alignment}.mol2", f"{name}.solv")
                        db2_data_all += db2_data
                write_to_file(db2_data_all, f"../{name}.db2")
                os.chdir("../..")
                if cleanup:
                    subprocess.run(f"rm -rf 3d/{name} solv/{name} db2/{name}", shell=True)
            except Exception as e:
                logger.error(f"Error in converting {name} to DB2 format: {e}")
                os.chdir("../..")
                log_error(row['smiles'], name)
                continue
        if args.timing: 
            mol2db2_time = time.time()
            logging_time += f'{name}, {embed_time-start}, {amsol_time-embed_time}, {sampling_time-amsol_time}, {mol2db2_time-sampling_time}, {mol2db2_time-start} \n'

    if cleanup:
        if not(smi2db2_utils.is_slurm_job()):
            # Remove empty directories if not running on SLURM
            smi2db2_utils.remove_empty_directories("3d")
            smi2db2_utils.remove_empty_directories("solv")
    if args.timing:
        with open('msani_timing.csv', 'a') as f:
            f.write(logging_time)