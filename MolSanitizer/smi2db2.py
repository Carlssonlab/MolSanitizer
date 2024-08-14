"""Based partly on: 
    https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py
    http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
    http://pubs.acs.org/doi/abs/10.1021/ci2004658
    https://greglandrum.github.io/rdkit-blog/posts/2024-02-11-more-multithreading.html

    Should try to sample all possible conformations based on dihedral angles sampling based on: https://github.com/dkoes/rdkit-scripts/blob/master/rdallconf.py
"""
# Author: Thua-Phong Lam, Jens Carlsson lab, Uppsala University
import pandas as pd
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolTransforms, rdDistGeom, rdMolAlign
import os
import subprocess
from pathlib import Path
import numpy as np
from MolSanitizer.amsol import run_amsol
from MolSanitizer.db2 import mol2db2, hydrogens, mol2
from MolSanitizer import strain_filter
import random

import logging
logger = logging.getLogger('molsani')

rotatable_pattern=r'[*]~[*;!$(*#*)!$([!#6&X2H])!$([!#6&X3H2])]-&!@[*;!$(*#*)!$([!#6&X2H])!$([!#6&X3H2])]~[*]'
Torlib = strain_filter.parse_torlib()
rigid_rule_files = Path(__file__).parent / 'Data' / 'rigid_part_rules.txt'
rigid_rules = pd.read_csv(rigid_rule_files, header=None, sep ='\s+', names=['SMARTS','label'])
rigid_rules['mol'] = rigid_rules['SMARTS'].apply(lambda x: Chem.MolFromSmarts(x))

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

def embed_smiles(smiles, name, rmsd=0.5, randomSeed=42, numConfs=300):
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
    params = rdDistGeom.ETKDGv3()
    params.numThreads = 0  # Use all available threads
    params.pruneRmsThresh = rmsd  # Prune conformations that are too similar
    params.randomSeed = randomSeed # For reproducibility

    mol = Chem.MolFromSmiles(smiles)
    mol_H = Chem.AddHs(mol)
    res = Chem.Mol(mol_H) # res = result molecule with conformations
    res.RemoveAllConformers() # An empty conformer list

    conf_energies = []
    
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(20) #1 means vacumn, 80 means water, 20 is the compromised value (still arbitrary)

    for cid in rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId=cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        conformer = mol_H.GetConformer(cid)
        conf_energies.append((energy, conformer))
    
    conf_energies = sorted(conf_energies, key=lambda x: x[0]) # sort by increasing E
    energy, conformer = conf_energies.pop(0) # get the lowest energy conformer
    res.AddConformer(conformer, assignId=True) # add it to the conformations list
    res.SetProp("_Name", name)
    netcharge = sum(atom.GetFormalCharge() for atom in res.GetAtoms())
    return res, Chem.MolToMolBlock(res), netcharge

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


def check_too_close_nonbonded_atoms(conformer, mol, threshold=1.6):
    """
    Check if any non-bonded atoms in a given conformer are too close to each other.

    Parameters:
    conformer (rdkit.Chem.rdchem.Conformer): The RDKit conformer object.
    mol (rdkit.Chem.Mol): The RDKit molecule object.
    threshold (float): The distance threshold below which atoms are considered too close.

    Returns:
    bool: True if any non-bonded atoms are too close, False otherwise.
    """
    positions = conformer.GetPositions()
    
    # Get the number of atoms in the molecule
    num_atoms = mol.GetNumAtoms()
    
    # Iterate over all pairs of atoms
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            # Check if the atoms are bonded
            if not mol.GetBondBetweenAtoms(i, j):
                # Calculate the distance between the non-bonded atoms
                distance = np.linalg.norm(positions[i] - positions[j])
                if distance < threshold:
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
    return(len(set(second_atom_symbols))==1) or (len(set(third_atom_symbols))==1)
   

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
    
    if VERBOSE: print(f"Converted and saved {molecule_count} conformations to {output_mol2} with assigned charges")


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
    reordered_rot_bonds = [] 

    match_torlib = strain_filter.get_match_dihedral(mol, Torlib)

    # We need to put the terminal rot bonds at the beginning as they dont 
    # contribute much to the diversity of the conformations
    for rotatable_bond in rotatable_bonds:
        if is_terminal(mol, rotatable_bond[1:3]):
            reordered_rot_bonds.append(rotatable_bond)
            if is_symmetric(mol, rotatable_bond[1:3]):
                #Exclude meaningless angles for symmetric terminal rotatable bonds (eg. CH3, CF3 only rotate onces)
                for rule in match_torlib:
                    if set(rotatable_bond[1:3]) == set(rule[1][1:3]):
                        while len(rule[2]) > 2: # rotate these angles by 120 degrees is the same
                            while (len(rule[2]) >= 2) and ((rule[2][0][0] - rule[2][1][0]) % 120 == 0): rule[2].pop(1)
                            if len(rule[2]) > 2 : rule[2].pop(-1)
                            

    for bond in rotatable_bonds:
        if bond not in reordered_rot_bonds:
            reordered_rot_bonds.append(bond)
    if VERBOSE: 
        print(reordered_rot_bonds)
        for i in match_torlib: print(i)
    num_confs = 1

    for bond in reordered_rot_bonds:
        for rule in match_torlib:
            if set(bond[1:3]) == set(rule[1][1:3]):
                num_confs *= (len(rule[2]))
                break
    
    return num_confs, reordered_rot_bonds, match_torlib

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

def torsional_scan_rand(mol, conf, i, matches, match_torlib, sdwriter, product, original, numConfs, atom_maps):
    '''recursively enumerate all angles for matches dihedrals.  i is where is
    which dihedral we are enumerating by degree to output conformers to out'''
    if len(product) >= numConfs: return product
    if i >= len(matches): #base case, torsions should be set in conf
        if check_too_close_nonbonded_atoms(mol.GetConformer(conf), mol): return product
        product.append(Chem.Conformer(mol.GetConformer(conf))) 
        rdMolAlign.AlignMol(mol, original, conf, 0, atomMap=[(i, i) for i in atom_maps])
        sdwriter.write(mol, conf)
        return product
    else:
        peaks = strain_filter.extract_peaks(match_torlib, matches[i][1:3])
        dihedral_4_atoms = peaks[0]
        for (prefered, tolerance) in peaks[1]:
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(conf),*dihedral_4_atoms,value = get_random_angle(prefered, tolerance))
            product = torsional_scan_rand(mol, conf, i+1, matches, match_torlib, sdwriter, product, original, numConfs, atom_maps)        
        return product
    

def stochastic_sampling(mol, name, reordered_rot_bonds, match_torlib, sdwriter, original, numConfs, atom_maps):
    product = []
    max_angles = 0 
    for rule in match_torlib:
        max_angles = max(len(rule[2]), max_angles)
    #visit_matrix = np.zeros((len(reordered_rot_bonds), max_angles))
    n_transform = len(reordered_rot_bonds)
    max_attempts = 1000
    attempts = 0
    while (len(product)< numConfs):
        for idx in range(n_transform):
            bond_idx = random.randint(0, len(reordered_rot_bonds)-1)
            bond = reordered_rot_bonds[bond_idx]

            peaks = strain_filter.extract_peaks(match_torlib, bond[1:3])
            peak_idx = random.randint(0, len(peaks[1])-1)
            peak = peaks[1][peak_idx]
            #visit_matrix[bond_idx][peak_idx] += 1
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(0),*peaks[0],value = get_random_angle(peak[0], peak[1]))
        #TODO: prune by RMSD?
        #TODO: if more than eg. 1000 failed attempts, break 
        if check_too_close_nonbonded_atoms(mol.GetConformer(0), mol): 
            attempts += 1
            if attempts > max_attempts: break
            continue
        attempts = 0
        rdMolAlign.AlignMol(mol, original, 0, 0, atomMap=[(i, i) for i in atom_maps])
        product.append(Chem.Conformer(mol.GetConformer(0)))
        sdwriter.write(mol,confId=0)
    #print(visit_matrix)
    return product

def find_rigid_part(mol):
    '''Find rigid parts of the molecule'''
    for rule in rigid_rules.itertuples():
        matches = mol.GetSubstructMatches(rule.mol)
        if len(matches) > 0:
            return matches[0]
    

def choose_sampling_method(mol, name, numConfs, VERBOSE=False):
   
    # Find rigid parts as anchor points for the molecules
    atom_maps = find_rigid_part(mol)
    if VERBOSE: print(f'Uses {atom_maps} as rigid part')
    # Count number of rotatable hydrogens and number of conformations contributed by them
    mol2_countH = mol2.Mol2(mol2fileName=f"{name}.mol2", nameFileName=None, mol2text=None)
    num_confs_H = hydrogens.count_confs_by_H(mol2_countH)

    # Divide the number of conformations by that contributed by rotatable hydrogens
    # This is a way to mimick the method of the previous approach of UCSF
    if num_confs_H > 30: num_confs_H = 30
    elif num_confs_H > 3: num_confs_H = 3
    numConfs = numConfs // num_confs_H 
    num_confs_by_rotbonds, reordered_rot_bonds, match_torlib = count_confs_by_rotbonds(mol, VERBOSE)
    
    if VERBOSE: print(num_confs_by_rotbonds, num_confs_H)
    rotated_file = f"{name}_rotated.sdf"
    sdwriter = Chem.SDWriter(rotated_file)
    original_mol = Chem.Mol(mol)
    if num_confs_by_rotbonds <= numConfs:
        if VERBOSE: print('Running systematic torsional scan')
        torsional_scan_rand(mol, conf=0, i=0, matches = reordered_rot_bonds, match_torlib = match_torlib,
                            sdwriter=sdwriter, product=list(), original=original_mol, numConfs = numConfs, atom_maps = atom_maps)
    else:
        if VERBOSE: print('Running stochastic torsional sampling')
        stochastic_sampling(mol, name, reordered_rot_bonds, match_torlib, sdwriter, original_mol, numConfs, atom_maps)
    sdwriter.close()
    convert_sdf_mol2(rotated_file, f"{name}_rotated.mol2", VERBOSE)
    #subprocess.run(f"rm {rotated_file}", shell=True)

def gen_conf_chunk(df: pd.DataFrame, randomSeed = 42, numConfs = 10000, VERBOSE = False, cleanup=False):
    
        env = setup_env()
        random.seed(randomSeed)
        if VERBOSE: print(env['LD_LIBRARY_PATH'])
        for idx, row in df.iterrows():
            # Embed smiles into initial conformation 
            # (300  conformers is inspired from https://pubs.acs.org/doi/abs/10.1021/ci2004658, then we only use the minimal energy one)
            if VERBOSE: print("Generating initial 3D conformations...")
            mol, molblock, netcharge = embed_smiles(row['smiles'], row['ids'], randomSeed = randomSeed, numConfs = 300)
            name = row['ids']
            mol2_block = convert(molblock, "mol", "mol2")

            # Solvation using AMSOL
            if VERBOSE: print("Solvating...")
            subprocess.run(f"mkdir -p solv/{name}", shell=True)
            write_to_file(mol2_block, f"solv/{name}/{name}.mol2")
            os.chdir(f"solv/{name}")
            run_amsol.prepare(f"{name}.mol2", name, netcharge)
            run_amsol.run('temp.in-hex', 'temp.o-hex', env)
            run_amsol.run('temp.in-wat', 'temp.o-wat', env)
            run_amsol.process_output('temp.o-wat', 'temp.o-hex', "temp.mol2", "output")
            subprocess.run(f"cp output.mol2 {name}_solv.mol2", shell=True)
            subprocess.run(f"mv output.solv {name}_solv.solv", shell=True)
            os.chdir("../..")

            # 3D generation
            if VERBOSE: print("3D generation...")
            subprocess.run(f"mkdir -p 3d/{name}", shell=True)
            subprocess.run(f"cp solv/{name}/{name}_solv.mol2 3d/{name}/{name}.mol2", shell=True)
            os.chdir(f"3d/{name}")

            choose_sampling_method(mol ,name, numConfs, VERBOSE)
            os.chdir("../..")

            # Mol2DB2
            if VERBOSE: print("Converting to DB2 format...")
            subprocess.run(f"mkdir -p db2/{name}", shell=True)
            subprocess.run(f"mv 3d/{name}/{name}_rotated.mol2 db2/{name}/{name}.mol2", shell=True)
            subprocess.run(f"mv solv/{name}/{name}_solv.solv db2/{name}/{name}.solv", shell=True)
            os.chdir(f"db2/{name}")

            if not os.path.isfile(f"{name}.mol2") or not os.path.isfile(f"{name}.solv"):
                logger.error("Not found mol2 files and solv for db2 generation")
                os.chdir("../..")
                continue
            else:
                try:
                    db2_data = mol2db2.mol2db2_quick(f"{name}.mol2", f"{name}.solv")
                    write_to_file(db2_data, f"../{name}.db2")
                    os.chdir("../..")
                    if cleanup:
                        subprocess.run(f"rm -rf 3d/{name} solv/{name}", shell=True)
                except Exception as e:
                    logger.error(f"Error in converting {name} to DB2 format {e}")
                    os.chdir("../..")
                    continue
                