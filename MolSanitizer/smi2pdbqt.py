try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
except ImportError:
    print("Please install the meeko package to use this script.")
from MolSanitizer import smi2db2_utils, strain_filter
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
import pandas as pd
import logging
import os
logger = logging.getLogger('molsani')

planar_lib, non_planar_lib = strain_filter.parse_sr_confs_library()

def embed_smiles_rdkit(smiles, name, randomSeed = 42, num_ring_confs = 1, VERBOSE = False) -> list:
    '''
    Adopt a similar method to smi2db2, but don't need to store redundant conformations just for AMSOL.
    If no aliphatic ring, embed 1 conformation. If multiple conformations, embed multiple ring conformations
    and use the ring conformation library to filter out "good" ones.
    '''
    mol_H = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mol_H.SetProp("_Name", name)
    empty_mol = Chem.Mol(mol_H)

    params = rdDistGeom.srETKDGv3()
    params.numThreads = 1  # Use all available threads
    params.pruneRmsThresh = 0.35  # Prune conformations that are too similar, not user-definable here
    params.randomSeed = randomSeed # For reproducibility
    params.useRandomCoords = True

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
    else: numConfs = 1
    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_H, mmffVariant="MMFF94s")
    mp.SetMMFFDielectricConstant(1) #1 means vacumn, 80 means water, 20 is the compromised value (still arbitrary)
    mp.SetMMFFEleTerm(False)

    conf_ring_descriptors_df = pd.DataFrame()
    try:
        for cid in rdDistGeom.EmbedMultipleConfs(mol_H, numConfs=numConfs, params=params):
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_H, mp, confId=cid)
            if conjugated_substituted_Ns:
                for a, b, c, d in conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5)
                for a, b, c, d in additional_conjugated_substituted_Ns: ff.MMFFAddTorsionConstraint(a, b, c, d, False, 178, 182, 1)
                for a, b, c, d, e in amide_linkages: # O=C-N(-C)-H should be coplanar.
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 1)
                    ff.MMFFAddTorsionConstraint(a, b, c, e, False, 180, 180, 1)
            ff.Minimize()
            conformer = mol_H.GetConformer(cid)
            energy = ff.CalcEnergy()
            conf_ring_descriptors_df = smi2db2_utils.classify_confs(conformer, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df)
    except: pass
    #conf_ring_descriptors_df.to_csv(f'{name}_confs.csv')
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
                for a, b, c, d, e in amide_linkages: # O=C-N(-C)-H should be coplanar (puckered because of MMFF94s).
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 1)
                    ff.MMFFAddTorsionConstraint(a, b, c, e, False, 180, 180, 1)
            ff.Minimize()
            conformer = mol_H.GetConformer(cid)
            energy = ff.CalcEnergy()
            conf_ring_descriptors_df = smi2db2_utils.classify_confs(conformer, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df)
    if num_ring_confs == 1:
        empty_mol.AddConformer(conformer, assignId=True)
        return [empty_mol]
    else:
        conf_ring_descriptors_df.sort_values('Energy', inplace=True)
        # Keep a reservoir as the lowest energy possible conformer in case no good ring conformers are found.
        rigid_scaffolds = []
        conf_ring_descriptors_df = smi2db2_utils.remove_unfavorable_confs(conf_ring_descriptors_df, name)
        temp_list = conf_ring_descriptors_df.values.tolist()
        while len(rigid_scaffolds) < num_ring_confs and temp_list:
            lowest_energy_entry = temp_list.pop(0)
            conformer, current_descriptors = lowest_energy_entry[0], lowest_energy_entry[2:-1]
            scaffold = Chem.Mol(empty_mol)
            scaffold.AddConformer(conformer, assignId=True)
            rigid_scaffolds.append(scaffold)
            temp_list = smi2db2_utils.ring_conf_clusters(current_descriptors, temp_list)
        if VERBOSE: print(f'\tBefore: {len(mol_H.GetConformers())}, after: {len(rigid_scaffolds)}')
        return rigid_scaffolds
def smi2pdbqt(smiles, name, randomSeed = 42, num_ring_confs = 1, VERBOSE = False):
    '''
    Convert a SMILES string to a PDBQT file.
    '''
    if VERBOSE: 
        print(f"Converting {name}...")
        print(f"\tSMILES: {smiles}")
        print(f"\tNumber of ring conformers: {num_ring_confs}")

    try:
        mols = embed_smiles_rdkit(smiles, name, randomSeed, num_ring_confs, VERBOSE)
        mkprep = MoleculePreparation()

        if not mols:
            logger.warning(f"Failed to embed {name}.")
            return
        is_multi = len(mols) > 1
        for i, mol in enumerate(mols):
            prepared_mol = mkprep(mol)
            pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(prepared_mol[0])
            if success:
                #print(pdbqt_string)
                if is_multi:
                    with open(f"{name}.nr{i}.pdbqt", 'w') as f:
                        for line in pdbqt_string:
                            f.write(line)
                else:
                    with open(f"{name}.pdbqt", 'w') as f:
                        for line in pdbqt_string:
                            f.write(line)
            else:
                print(error_msg)
            
    except Exception as e:
        logger.error(f"Failed to convert {name}. Error: {e}")
        return

def gen_conf_chunk(df: pd.DataFrame, args):
    '''
    Generate conformers for a chunk of the input dataframe.
    '''
    # Test mode in unittest, not to produce redundant files here
    if args.test: 
        os.chdir(args.prefix)
    os.makedirs('pdbqt', exist_ok=True)
    os.chdir('pdbqt')
    for i, row in df.iterrows():
        name = row['ids']
        smiles = row['smiles']
        smi2pdbqt(smiles, name, args.randomSeed, args.nringconfs, args.debug)
    os.chdir('..')