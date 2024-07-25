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
from MolSanitizer.db2 import mol2db2

import logging
logger = logging.getLogger('molsani')

internal_pattern=r'[*;!$([#1])]~[*;!$(*#*)]-&!@[*;!$(*#*)]~[*;!$([#1])]'
terminal_pattern=r'*~[*;!$(*#*)]-[$(A([#1])([#1])[#1]),$(A(F)(F)F),$(A(Cl)(Cl)Cl),$(A(Br)(Br)Br),$(A(I)(I)I)]~[#1,F,Cl,Br,I,CH3]'

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

def embed_smiles(smiles, name, rmsd=0.25, randomSeed=42, numConfs=50):
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
    mp.SetMMFFDielectricConstant(80)

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
    matches = mol.GetSubstructMatches(qmol);
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

def genConformer_terminal_r(mol, conf, i, matches,  sdwriter, product, original, degree = 60, maxconf = 10000):
    '''recursively enumerate all angles for matches dihedrals.  i is where is
    which dihedral we are enumerating by degree to output conformers to out'''
    if len(product) > maxconf: return product
    if i >= len(matches): #base case, torsions should be set in conf
        if check_too_close_nonbonded_atoms(mol.GetConformer(conf), mol): return product
        product.append(Chem.Conformer(mol.GetConformer(conf)))
        sdwriter.write(mol,conf)
        return product
    else:
        deg = rdMolTransforms.GetDihedralDeg(original,*matches[i])
        for idx in range(int(120/degree)): # 120 is the maximum angle for symmetrical terminal dihedrals
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(conf),*matches[i],value = deg + degree*idx)
            product = genConformer_terminal_r(mol, conf, i+1, matches, sdwriter, product, original, degree, maxconf)
        return product
    
def genConformer_internal_r(mol, conf, i, internal_matches, terminal_matches, sdwriter, product, original, degree = 60, maxconf = 10000):
    '''recursively enumerate all angles for matches dihedrals.  i is where is
    which dihedral we are enumerating by degree to output conformers to out'''
    if i >= len(internal_matches): #base case, torsions should be set in conf
        product = genConformer_terminal_r(mol, conf, 0, terminal_matches, sdwriter, product, original, degree, maxconf)
        return product
    else:
        deg = rdMolTransforms.GetDihedralDeg(original,*internal_matches[i])
        for idx in range(int(360/degree)):
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(conf),*internal_matches[i],value = deg + degree*idx)
            product = genConformer_internal_r(mol, conf, i+1, internal_matches, terminal_matches, sdwriter, product, original, degree, maxconf)
        return product
    
def convert_file(input_file, output_file, input_format, output_format):
    with open(input_file, "r") as f:
        data = f.read()
    converted_data = convert(data, input_format, output_format)
    with open(output_file, "w") as f:
        f.write(converted_data)

def extract_partial_charges(mol2_file, output_file='partialcharges.txt', VERBOSE=False):
    charges = []
    with open(mol2_file, 'r') as f:
        atom_section = False
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                atom_section = True
                continue
            if ('@<TRIPOS>' in line) and ('ATOM' not in line) and (atom_section == True):
                atom_section = False
                break
            if atom_section:
                parts = line.split()
                if len(parts) >= 9:
                    charges.append(float(parts[8]))
    with open(output_file, 'w') as f:
        for charge in charges:
            f.write(f"{charge}\n")
    
    if VERBOSE: print(f"Partial charges extracted and saved to {output_file}")

def assign_charges_and_convert(sdf_file, charges_file, output_mol2):
    # Read charges
    with open(charges_file, 'r') as f:
        charges = [float(line.strip()) for line in f]
    
    # Set up OpenBabel conversion
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "mol2")
        
        # Read all molecules from the SDF file
    mol = ob.OBMol()
    if not obConversion.ReadFile(mol, sdf_file):
        print(f"Error reading SDF file: {sdf_file}")
        return
    unprocessed_mol2_str = ""
    molecule_count = 0
    with open(output_mol2, 'w') as out_file:
        while True:
            molecule_count += 1
            # Assign charges
            for i, atom in enumerate(ob.OBMolAtomIter(mol)):
                if i < len(charges):
                    atom.SetPartialCharge(charges[i])
            mol.SetAutomaticPartialCharge(False)
            mol.SetPartialChargesPerceived()
            
            # Convert molecule to MOL2 format and get as string
            mol2_str = obConversion.WriteString(mol)
            unprocessed_mol2_str+=(mol2_str)
            # Write to output file
            out_file.write(mol2_str)
            
            # Clear the molecule and read the next one
            mol.Clear()
            if not obConversion.Read(mol):
                break
    
    print(f"Converted and saved {molecule_count} conformations to {output_mol2} with assigned charges")
    return unprocessed_mol2_str

def torsional_scan(mol, name, maxconf = 10000, partial_charges='partialcharges.txt'):
    extract_partial_charges(f"{name}.mol2", partial_charges)
    res = Chem.Mol(mol); 
    res.RemoveAllConformers()
    internal_matches = getDihedralMatches(mol, internal_pattern)
    terminal_matches = getDihedralMatches(mol, terminal_pattern)
    #print(rotmatches)
    rotated_file = f"{name}_rotated.sdf"
    sdwriter = Chem.SDWriter(rotated_file)
    product = genConformer_internal_r(mol, conf=0, i=0, internal_matches = internal_matches, terminal_matches = terminal_matches, sdwriter=sdwriter, product=list(), original=mol.GetConformer(0), degree=60, maxconf=maxconf)
    for conf in product:
        res.AddConformer(conf, assignId = True)
    unprocessed_mol2_str = assign_charges_and_convert(rotated_file, partial_charges, f"{name}_rotated.mol2")
    subprocess.run(f"rm {rotated_file}", shell=True)

    return res, unprocessed_mol2_str


def gen_conf_chunk(df: pd.DataFrame, rmsd = 0.25, randomSeed = 42, numConfs = 10000, VERBOSE = False):
        mol2_data = []

        env = setup_env()
        db2_all_data = ""
        if VERBOSE: print(env['LD_LIBRARY_PATH'])
        for idx, row in df.iterrows():
            # Embed smiles into initial conformation
            if VERBOSE: print("Embedding...")
            mol, molblock, netcharge = embed_smiles(row['smiles'], row['ids'], rmsd = rmsd, randomSeed = randomSeed, numConfs = numConfs)
            name = row['ids']
            mol2 = convert(molblock, "mol", "mol2")

            # Solvation using AMSOL
            if VERBOSE: print("Solvating...")
            subprocess.run(f"mkdir -p solv/{name}", shell=True)
            write_to_file(mol2, f"solv/{name}/{name}.mol2")
            os.chdir(f"solv/{name}")
            run_amsol.prepare(f"{name}.mol2", name, netcharge)
            run_amsol.run('temp.in-hex', 'temp.o-hex', env)
            run_amsol.run('temp.in-wat', 'temp.o-wat', env)
            run_amsol.process_output('temp.o-wat', 'temp.o-hex', "temp.mol2", "output")
            subprocess.run(f"cp output.mol2 {name}_solv.mol2", shell=True)
            subprocess.run(f"mv output.solv {name}_solv.solv", shell=True)
            os.chdir("../..")

            # Torsion scan
            if VERBOSE: print("Torsional scan...")
            subprocess.run(f"mkdir -p 3d/{name}", shell=True)
            subprocess.run(f"cp solv/{name}/{name}_solv.mol2 3d/{name}/{name}.mol2", shell=True)
            os.chdir(f"3d/{name}")
            torsional_scan(mol ,name, numConfs)
            
            # Torsional filtering
            # We need an rdkit Mol with generated conformers to filter the torsion outs.
            os.chdir("../..")
            # Mol2DB2
            if VERBOSE: print("Converting to DB2 format...")
            subprocess.run(f"mkdir -p db2/{name}", shell=True)
            subprocess.run(f"mv 3d/{name}/{name}_rotated.mol2 db2/{name}/{name}.mol2", shell=True)
            subprocess.run(f"mv solv/{name}/{name}_solv.solv db2/{name}/{name}.solv", shell=True)
            os.chdir(f"db2/{name}")
            db2_data = mol2db2.mol2db2_quick(f"{name}.mol2", f"{name}.solv")
            write_to_file(db2_data, f"{name}.db2")
            os.chdir("../..")
