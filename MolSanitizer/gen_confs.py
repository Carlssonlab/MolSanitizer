"""Based partly on: 
    https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py
    http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
    http://pubs.acs.org/doi/abs/10.1021/ci2004658
    https://greglandrum.github.io/rdkit-blog/posts/2024-02-11-more-multithreading.html

    Should try to sample all possible conformations based on dihedral angles sampling based on: https://github.com/dkoes/rdkit-scripts/blob/master/rdallconf.py
"""
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

import logging
logger = logging.getLogger('molsani')


def remove_nonpolar_hydrogens(mol):
    """
    Remove non-polar hydrogen atoms from the molecule.
    Returns a new molecule with only polar hydrogens retained.
    """
    # Create an editable mol object
    edit_mol = Chem.RWMol(mol)
    
    # Identify polar hydrogens
    polar_h_indices = set()
    for atom in edit_mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # It's a hydrogen
            neighbor = atom.GetNeighbors()[0]  # Get the atom it's bonded to
            if neighbor.GetAtomicNum() in [7, 8, 16]:  # N, O, or S
                polar_h_indices.add(atom.GetIdx())
    
    # Remove non-polar hydrogens
    atoms_to_remove = []
    for atom in edit_mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetIdx() not in polar_h_indices:
            atoms_to_remove.append(atom.GetIdx())
    
    # Remove atoms in reverse order to avoid index issues
    for idx in sorted(atoms_to_remove, reverse=True):
        edit_mol.RemoveAtom(idx)
    
    # Convert back to a regular mol object and return
    return edit_mol.GetMol()

def rmsd_filter(mol, ref_conf, conf_energies, threshold):
    """https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.GetConformerRMS"""
    # we use heavy atoms RMSD; not all atoms (Peter Gedeck's suggestion)
    mol_noH = remove_nonpolar_hydrogens(mol)
    ref_conf_id = ref_conf.GetId()
    res = []
    for e, curr_conf in conf_energies:
        curr_conf_id = curr_conf.GetId()
        rms = AllChem.GetBestRMS(mol_noH, mol_noH, ref_conf_id, curr_conf_id, numThreads=0) #Use this to avoid the symmetrical problem
        if rms > threshold:
            res.append((e, curr_conf))
    return res

def write_to_sdf(mol, filename):
    with Chem.SDWriter(filename+'.sdf') as writer:
        for i, conf in enumerate(mol.GetConformers()):
            rdMolAlign.AlignMol(mol, mol, prbCid = conf.GetId(), refCid = 0)
            mol.SetProp("_Name", f"Conformer_{i+1}")
            writer.write(mol, confId=conf.GetId())

def generate_conformers(smiles, name, rmsd = 0.25, randomSeed = 42, numConfs = 2000):

    params = AllChem.ETKDGv3()
    params.numThreads = 0  # Use all available threads
    params.pruneRmsThresh = rmsd  # Prune conformations that are too similar
    params.randomSeed = randomSeed # For reproducible
    params.maxIterations = 20000

    mol = Chem.MolFromSmiles(smiles)
    mol_H = Chem.AddHs(mol)
    res = Chem.Mol(mol_H) # res = result molecule with conformations
    res.RemoveAllConformers() # An empty conformer list

    conf_energies = []
    
    """ This part is for writing out sdf before optimization
    print(len(mol_H.GetConformers())
    write_to_sdf(mol_H, name+"_noopt.sdf")
    """

    mp = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(80)

    for cid in Chem.rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
        ff = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId = cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        conformer = mol_H.GetConformer(cid)
        conf_energies.append((energy, conformer))
    
    
    conf_energies = sorted(conf_energies, key = lambda x: x[0]) # sort by increasing E
    while conf_energies:
        energy, conformer = conf_energies.pop(0) # get the lowest energy conformer
        res.AddConformer(conformer, assignId = True) # add it to the conformations list
        conf_energies = rmsd_filter(mol_H, conformer, conf_energies, rmsd) # remove all conformers that are too similar to it

    #print(len(res.GetConformers()))
    write_to_sdf(res, name)
    

def gen_conf_chunk(df: pd.DataFrame, rmsd: float, randomSeed: int, numConfs: int):
        for idx, row in df.iterrows():
            generate_conformers(row['smiles'], row['ids'], rmsd = rmsd, randomSeed = randomSeed, numConfs = numConfs)