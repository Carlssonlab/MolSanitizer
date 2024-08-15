# %%
"""Based partly on: 
    https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py
    http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
    http://pubs.acs.org/doi/abs/10.1021/ci2004658
    https://greglandrum.github.io/rdkit-blog/posts/2024-02-11-more-multithreading.html

    Should try to sample all possible conformations based on dihedral angles sampling based on: https://github.com/dkoes/rdkit-scripts/blob/master/rdallconf.py
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms, rdDistGeom
import math
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

def getDihedralMatches(mol):
    '''return list of atom indices of dihedrals'''
    #this is rdkit's "strict" pattern
    pattern = r"*~[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]~*"
    qmol = Chem.MolFromSmarts(pattern)
    matches = mol.GetSubstructMatches(qmol);
    #these are all sets of 4 atoms, uniquify by middle two
    uniqmatches = []
    seen = set()
    for (a,b,c,d) in matches:
        if (b,c) not in seen:
            seen.add((b,c))
            uniqmatches.append((a,b,c,d))
    return uniqmatches

def genConformer_r(mol, conf, i, matches,  sdwriter, product, degree = 30):
    '''recursively enumerate all angles for matches dihedrals.  i is where is
    which dihedral we are enumerating by degree to output conformers to out'''
    if i >= len(matches): #base case, torsions should be set in conf
        product.append(conf)
        #sdwriter.write(mol,conf)
        return product
    else:
        total = 0
        deg = 0
        while deg < 360.0:
            rad = math.pi*deg / 180.0
            rdMolTransforms.SetDihedralRad(mol.GetConformer(conf),*matches[i],value=rad)
            total += genConformer_r(mol, conf, i+1, matches, sdwriter, product, degree)
            deg += degree
        return product
    
def generate_conformers(smiles, name, rmsd = 0.25, randomSeed = 42, numConfs = 2000):

    params = rdDistGeom.ETKDGv3()
    params.numThreads = 0  # Use all available threads
    params.pruneRmsThresh = rmsd  # Prune conformations that are too similar
    params.randomSeed = randomSeed # For reproducible
    params.maxIterations = 20000

    mol = Chem.MolFromSmiles(smiles)
    mol_H = Chem.AddHs(mol)
    res = Chem.Mol(mol_H) # res = result molecule with conformations
    res.RemoveAllConformers() # An empty conformer list

    conf_energies = []
    
    #This part is for writing out sdf before optimization
    #print(len(mol_H.GetConformers())
    #write_to_sdf(mol_H, name+"_noopt")
    

    mp = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(80)

    for cid in Chem.rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
        ff = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId = cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        conformer = mol_H.GetConformer(cid)
        conf_energies.append((energy, conformer))
    
    #write_to_sdf(mol_H, name+"_noopt")
    
    conf_energies = sorted(conf_energies, key = lambda x: x[0]) # sort by increasing E
    
    # This part is for outputing the minimized conformations
    while (len(conf_energies) > 0):
        energy, conformer = conf_energies.pop(0) # get the lowest energy conformer
        res.AddConformer(conformer, assignId = True) # add it to the conformations list
        conf_energies = rmsd_filter(mol_H, conformer, conf_energies, rmsd) # remove all conformers that are too similar to it

    print(f'{name}: before: {len(mol_H.GetConformers())}, after: {len(res.GetConformers())}')
    write_to_sdf(res, f'tricky_ligand_cleaned/{name}')

    """rotmatches = getDihedralMatches(res)
    total = 0
    while (len(conf_energies) > 0):
        energy, conformer = conf_energies.pop(0) # get the lowest energy conformer
        for m in rotmatches:
            rdMolTransforms.SetDihedralRad(conformer,*m,value=0) # set all dihedrals to 0
        res.AddConformer(conformer, assignId = True) # add it to the conformations list
        conf_energies = rmsd_filter(mol_H, conformer, conf_energies, rmsd) # remove all conformers that are too similar to it

    sdwriter = Chem.SDWriter(f"{name}_rotated.sdf")    
    for conf in range(len(res.GetConformers())):
        product = list(genConformer_r(res, conf, 0, rotmatches, sdwriter, degree=30))
        """    

def gen_conf_chunk(df: pd.DataFrame, rmsd: float, randomSeed: int, numConfs: int):
        for idx, row in df.iterrows():
            generate_conformers(row['smiles'], row['ids'], rmsd = rmsd, randomSeed = randomSeed, numConfs = numConfs)

# %%
df = pd.read_csv('tricky_ligands_new_clean.smi', sep=r'\s+', names=['smiles', 'ids'], header=None)
gen_conf_chunk(df, rmsd = 0.25, randomSeed=200, numConfs=2000)


# %% [markdown]
# 

# %%
def genConformer_r(mol, conf, i, matches,  sdwriter, product, degree = 30):
    '''recursively enumerate all angles for matches dihedrals.  i is where is
    which dihedral we are enumerating by degree to output conformers to out'''
    if i >= len(matches): #base case, torsions should be set in conf
        product.append(Chem.Conformer(mol.GetConformer(conf)))
        sdwriter.write(mol,conf)
        return product
    else:
        deg = 0
        while deg < 360.0:
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(conf),*matches[i],value=deg)
            product = genConformer_r(mol, conf, i+1, matches, sdwriter, product, degree)
            deg += degree
        return product

def rmsd_filter2(mol, ref_conf, conf_energies, threshold):
    """https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.GetConformerRMS"""
    # we use heavy atoms RMSD; not all atoms (Peter Gedeck's suggestion)
    #mol_noH = remove_nonpolar_hydrogens(mol)
    mol_noH = mol
    ref_conf_id = ref_conf.GetId()
    res = []
    for curr_conf in conf_energies:
        curr_conf_id = curr_conf.GetId()
        rms = AllChem.GetBestRMS(mol_noH, mol_noH, ref_conf_id, curr_conf_id, numThreads=0) #Use this to avoid the symmetrical problem
        print(rms)
        if rms > threshold:
            res.append((curr_conf))
    return res

def generate_conformers_2(smiles, name, rmsd = 0.25, randomSeed = 42, numConfs = 2000):

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
    
    #This part is for writing out sdf before optimization
    #print(len(mol_H.GetConformers())
    #write_to_sdf(mol_H, name+"_noopt")
    

    mp = Chem.rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(80)

    for cid in Chem.rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
        ff = Chem.rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId = cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        conformer = mol_H.GetConformer(cid)
        conf_energies.append((energy, conformer))
    
    #write_to_sdf(mol_H, name+"_noopt")
    
    conf_energies = sorted(conf_energies, key = lambda x: x[0]) # sort by increasing E
    
    # This part is for outputing the minimized conformations
    while (len(conf_energies) > 0):
        energy, conformer = conf_energies.pop(0) # get the lowest energy conformer
        res.AddConformer(conformer, assignId = True) # add it to the conformations list
        conf_energies = rmsd_filter(mol_H, conformer, conf_energies, rmsd) # remove all conformers that are too similar to it

    print(f'{name}: before: {len(mol_H.GetConformers())}, after: {len(res.GetConformers())}')
    #write_to_sdf(res, f'tricky_ligand/{name}')

    rotmatches = getDihedralMatches(res)
    total = 0
    while (len(conf_energies) > 0):
        energy, conformer = conf_energies.pop(0) # get the lowest energy conformer
        for m in rotmatches:
            rdMolTransforms.SetDihedralRad(conformer,*m,value=0) # set all dihedrals to 0
        res.AddConformer(conformer, assignId = True) # add it to the conformations list
        conf_energies = rmsd_filter(mol_H, conformer, conf_energies, rmsd) # remove all conformers that are too similar to it

    sdwriter = Chem.SDWriter(f"{name}_rotated.sdf")    
    for conf in range(len(res.GetConformers())):
        product = genConformer_r(res, conf, 0, rotmatches, sdwriter, list(), degree=30)
    print(product)
    print(len(product))
    

    res_rotated = Chem.Mol(mol_H)
    res_rotated.RemoveAllConformers()
    res_rotated_dedup=Chem.Mol(res_rotated)
    for conf in product:
        res_rotated.AddConformer(conf, assignId = True)

    while (product):
        conf = product.pop(0)
        
        res_rotated_dedup.AddConformer(conf, assignId = True)
        product = rmsd_filter2(mol_H, conf, product, rmsd)
    print(len(res_rotated_dedup.GetConformers()))
    #print(AllChem.GetBestRMS(res_rotated, res_rotated, 0, 32, numThreads=0))
    #write_to_sdf(res_rotated_dedup, 'methyl_aniline_rotated_deduplicated')

generate_conformers_2('c1ccccc1NC', 'methyl_aniline')

# %%



