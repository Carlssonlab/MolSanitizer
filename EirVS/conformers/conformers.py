"""Based partly on: 
    https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py
    http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
    http://pubs.acs.org/doi/abs/10.1021/ci2004658
    https://greglandrum.github.io/rdkit-blog/posts/2024-02-11-more-multithreading.html

    Should try to sample all possible conformations based on dihedral angles sampling based on: https://github.com/dkoes/rdkit-scripts/blob/master/rdallconf.py
"""
# Author: Thua-Phong Lam, Jens Carlsson lab, Uppsala University
# Date: 2025-05-05

from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolAlign, rdMolTransforms, PropertyPickleOptions
from pathlib import Path

from EirVS.conformers import utils, mol2writer
from EirVS.filtering import strain_filter, filters
from EirVS.amsol import run_amsol
from EirVS.db2 import solv

import pandas as pd
import numpy as np
import logging
import os
import multiprocessing
import shutil
import random
import itertools
import tarfile, io
import time
import argparse

logger = logging.getLogger('eirvs')

planar_lib, non_planar_lib = strain_filter.parse_sr_confs_library()

class ConformerGenerator:
    '''
    Class to generate conformers from SMILES strings.
    
    Examples:
    >>> from EirVS.conformers.conformers import ConformerGenerator
    >>> smiles = 'CC(=O)C1=CC=CC=C1C(=O)O'
    >>> name = 'test'
    >>> confgen = ConformerGenerator(smiles, name, randomSeed=42, method='rdkit', num_ring_confs=1, numcores=1)
    >>> confgen.conf_sampling(numConfs=2000, energywindow=25, ignoreTorlib=False, AllowNonRing=False)
    >>> confgen.to_pdbqt()
    >>> confgen.to_sdf()
    >>> confgen.to_mol2()
    >>> confgen.to_db2()
    '''
    def __init__(self, 
                 smiles, 
                 name='test', 
                 forcefield='MMFF94s',
                 method='rdkit',
                 randomSeed=42, 
                 num_ring_confs=1, 
                 numcores=1, 
                 request_alignment=None,
                 ignoreTorlib=False,
                 threshold=1.6,
                 mode='vs',
                 tolerance=30,
                 VERBOSE=False):
        """
        Initialize the ConformerGenerator object.
        
        Args:
            smiles (str): SMILES string of the molecule
            name (str): Name of the molecule
            forcefield (str): Force field to use ('MMFF94', 'MMFF94s', or 'UFF')
            method (str): Embedding method to use ('rdkit', 'obabel', or 'corina')
            randomSeed (int): Random seed for reproducibility
            num_ring_confs (int): Number of ring conformers to generate
            numcores (int): Number of CPU cores to use
            request_alignment (str): Ring alignment in SMILES or SMARTS to support constrained docking
            ignoreTorlib (bool): Whether to ignore the torsion library
            threshold (float): Threshold in Angstrom for non-bonded atom distance.
            mode (str): Mode for conformer sampling ('vs' - virtual screening, 'extensive', 'ignoretorlib')
            tolerance (float): Tolerance for dihedral angle sampling.
            VERBOSE (bool): Whether to print verbose output
        """
        # Validate inputs
        self.smiles = smiles
        self.name = name
        self.method = method
        self.randomSeed = randomSeed
        random.seed(self.randomSeed)
        self.num_ring_confs = num_ring_confs
        self.numcores = numcores
        self.request_alignment = request_alignment
        self.ignoreTorlib = ignoreTorlib
        self.conf_sampled = False
        self.threshold = threshold
        self.mode = mode
        self.tolerance = tolerance
        self.VERBOSE = VERBOSE

        # Initialize molecule
        self._initialize_molecule()
        self.substructure_perception()
        
        # Calculate net charge
        self.netcharge = sum(atom.GetFormalCharge() for atom in self.mol.GetAtoms())
        
        # Initialize output containers
        self.ring_confs = [] # List of RDKit Mol with different ring conformations
        
        # Setup force field
        self._setup_forcefield(forcefield)

        if self.method == 'corina':
            self.embed_smiles_corina()
        elif self.method == 'obabel':
            self.embed_smiles_babel()
        elif self.method == 'rdkit':
            self.embed_smiles_rdkit()
        else:
            raise ValueError(f"Invalid embedding method: {self.method}. Supported methods are: rdkit, obabel, corina.")
        
    @classmethod
    def from_existing_data(cls, smiles, name, amsol_mol, ring_confs = None, mol2_str = None, request_alignment = None, mode:str = 'vs', tolerance = 30, VERBOSE=False):
        """Alternative constructor that initializes from existing data"""
        # Create a minimal instance
        instance = cls(smiles, name=name)
        
        # Override the instance attributes
        mol = Chem.Mol(amsol_mol)
        instance.amsol_mol = mol
        instance.mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(instance.amsol_mol, mmffVariant="MMFF94s")
        instance.method = 'rdkit'
        instance.ring_confs = [Chem.Mol(ring_conf) for ring_conf in ring_confs] if ring_confs else []
        instance.mol2_str = mol2_str
        instance.sulfo_matches = utils.find_sulfonamide_like_scaffolds(instance.amsol_mol)
        instance.request_alignment = request_alignment
        instance.atom_maps, instance.label_map = utils.find_rigid_part(mol, request_alignment)
        instance.VERBOSE = VERBOSE
        instance.mode = mode
        instance.tolerance = tolerance
        return instance
    
    def _initialize_molecule(self):
        """Initialize the molecule from SMILES"""
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError(f"Invalid SMILES string: {self.smiles}")
        self.mol_H = Chem.AddHs(self.mol)
        self.mol_H.SetProp("_Name", self.name)
        self.empty_mol = Chem.Mol(self.mol_H) # An RDKit Mol Object with no confs
        self.amsol_mol = Chem.Mol(self.mol_H) # An RDKit Mol Object with upto 10 confs for AMSOL
        
    def _setup_forcefield(self, forcefield):
        """Setup the force field for the molecule"""
        if forcefield not in ['MMFF94', 'MMFF94s', 'UFF']:
            raise ValueError(f"Invalid forcefield: {forcefield}. Supported forcefields are: MMFF94, MMFF94s, UFF.")
        
        if forcefield.startswith('MMFF'):
            self.mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(self.mol_H, mmffVariant=forcefield)
        # elif forcefield == 'UFF': Not supported now...
        #    self.mp = rdForceFieldHelpers.UFFGetMoleculeProperties(self.mol_H)
        
        # Store forcefield type for later use
        self.forcefield = forcefield

    def substructure_perception(self):
        """Identify important substructures in the molecule"""
        # Get ring information
        self.planar_rings, self.non_planar_rings = utils.get_flexible_ring(
            self.mol_H, planar_lib, non_planar_lib
        )
        
        # Get other important substructures.
        # Some of them are for correctures of MMFF94s
        self.sulfo_matches = utils.find_sulfonamide_like_scaffolds(self.mol_H)
        self.conjugated_substituted_nitrogen_5aro = utils.find_conjugated_substituted_nitrogen_5aro(self.mol_H)
        self.conjugated_substituted_nitrogen_6aro = utils.find_conjugated_substituted_nitrogen_6aro(self.mol_H)
        self.barbiturate_matches = utils.find_barbiturates(self.mol_H)
        self.hydantoin_matches = utils.find_hydantoins(self.mol_H)
        self.substituted_N_barbi_hydan_like = utils.find_substituted_N_barbi_hydan_like(self.mol_H, self.barbiturate_matches, self.hydantoin_matches)
        self.amide_linkages = utils.find_amide(self.mol_H)
        
        # Only find flippable Ns if we need multiple conformations
        self.flippable_Ns = utils.find_flipped_nitrogen(self.mol_H)
        self.flippable_Cs = utils.find_flipped_carbon(self.mol_H)
        
            
        if self.VERBOSE:
            if self.planar_rings:
                print('\t Found planar rings:')
                print(f'\t{self.planar_rings}')
            if self.non_planar_rings:
                print('\t Found non_planar rings:')
                for ring in self.non_planar_rings: 
                    print(f'\t {ring[0]} {ring[1]}')
            if self.flippable_Ns:
                print('\tFound flippable N structures')
                for match in self.flippable_Ns: 
                    print(f'\t {match}')
            if self.flippable_Cs:
                print('\tFound flippable C structures')
                for match in self.flippable_Cs: 
                    print(f'\t {match}')
            if self.sulfo_matches:
                print('\tFound sulfonamide-like structures')
                for match in self.sulfo_matches: 
                    print(f'\t {match}')
            if self.conjugated_substituted_nitrogen_5aro:
                print('\tFound conjugated substituted nitrogen in 5-membered aromatic rings')
                for match in self.conjugated_substituted_nitrogen_5aro: 
                    print(f'\t {match}')
            if self.conjugated_substituted_nitrogen_6aro:
                print('\tFound conjugated substituted nitrogen in 6-membered aromatic rings')
                for match in self.conjugated_substituted_nitrogen_6aro: 
                    print(f'\t {match}')
            if self.barbiturate_matches:
                print('\tFound barbiturate-like structures')
                for match in self.barbiturate_matches: 
                    print(f'\t {match}')
            if self.hydantoin_matches:
                print('\tFound hydatoin-like structures')
                for match in self.hydantoin_matches: 
                    print(f'\t {match}')
            if self.substituted_N_barbi_hydan_like:
                print('\tFound substituted N barbiturate/hydantoin-like structures')
                for match in self.substituted_N_barbi_hydan_like:
                    print(f'\t {match}')
            if self.amide_linkages:
                print('\tFound amide linkages')
                for match in self.amide_linkages: 
                    print(f'\t {match}')
            
        # Determine number of initial conformations needed
        self.num_initialConfs = 50 if (self.sulfo_matches or self.non_planar_rings or self.flippable_Ns) else 10

    #================= Generation of initial conformers =========================
    # There are three methods for the generation of the initial conformers:
    # 1. RDKit
    # 2. Open Babel
    # 3. Corina    
    def embed_smiles_rdkit(self):
        """
        Embeds a SMILES string into molecular conformers using RDKit,
        applying various quality checks and corrections inherited from the MMFF94s force field.
        """
        def embed_fix_ring_confs(method: str ='srETKDGv3'):
            """
            Generate initial conformers for the molecule using RDKit.
            Parameters:
                method (str): The method to use for embedding ('srETKDGv3' or 'ETKDGv3').
                """
            if method == 'srETKDGv3': params = rdDistGeom.srETKDGv3()
            else: params = rdDistGeom.ETKDGv3()
            params.numThreads = self.numcores
            if self.mol.GetNumHeavyAtoms() > 15: 
                params.pruneRmsThresh = 0.35 # An arbitrary threshold for small molecules and fragments, 
                                             # the RMSD pruning maynot be suitable anymore
            params.randomSeed = self.randomSeed # For reproducibility
            params.useRandomCoords = True
            conf_ring_descriptors_df = pd.DataFrame()
            for cid in rdDistGeom.EmbedMultipleConfs(self.mol_H, numConfs=self.num_initialConfs, params=params):
                ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(self.mol_H, self.mp, confId=cid)
                if self.conjugated_substituted_nitrogen_5aro:
                    # *-[nX3&+0:1]1[a:2][a:3][a:4][a:5]1 
                    # a-b-c-d -> 180; a-b-f-e -> 180
                    # b-c-d-e -> 0; d-e-f-b -> 0
                    for a, b, c, d, e, f in self.conjugated_substituted_nitrogen_5aro:
                        ff.MMFFAddTorsionConstraint(a, b, c, d, False, 180, 180, 1)
                        ff.MMFFAddTorsionConstraint(a, b, f, e, False, 180, 180, 1)
                        ff.MMFFAddTorsionConstraint(b, c, d, e, False, 0, 0, 5)
                        ff.MMFFAddTorsionConstraint(d, e, f, b, False, 0, 0, 5)
                if self.conjugated_substituted_nitrogen_6aro:
                    # *-[nX3&+0:1]1[a:2][a:3][a:4][a:5][a:6]1
                    # a-b-c-d -> 180; a-b-g-f -> 180
                    # b-c-d-e -> 0; e-f-g-b -> 0
                    for a, b, c, d, e, f, g in self.conjugated_substituted_nitrogen_6aro:
                        ff.MMFFAddTorsionConstraint(a, b, c, d, False, 180, 180, 1)
                        ff.MMFFAddTorsionConstraint(a, b, g, f, False, 180, 180, 1)
                        ff.MMFFAddTorsionConstraint(b, c, d, e, False, 0, 0, 5)
                        ff.MMFFAddTorsionConstraint(e, f, g, b, False, 0, 0, 5)
                if self.barbiturate_matches:
                    for match in self.barbiturate_matches:
                        n = len(match)
                        for i in range(n):
                            a, b, c, d = [match[(i + j) % n] for j in range(4)]
                            ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5)
                if self.hydantoin_matches:
                    for match in self.hydantoin_matches:
                        n = len(match)
                        for i in range(n):
                            a, b, c, d  = [match[(i + j) % n] for j in range(4)]
                            ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5)                
                if self.substituted_N_barbi_hydan_like:
                    for a, b, c, d in self.substituted_N_barbi_hydan_like:
                        ff.MMFFAddTorsionConstraint(a, b, c, d, False, 180, 180, 1)
                if self.amide_linkages:
                    for a, b, c, d, e in self.amide_linkages: # O=C-N(-C)-H should be coplanar.
                        ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 1)
                        ff.MMFFAddTorsionConstraint(a, b, c, e, False, 180, 180, 1)
                ff.Minimize()
                conformer = self.mol_H.GetConformer(cid)
                energy = ff.CalcEnergy()
                conf_ring_descriptors_df = utils.classify_confs(conformer, 
                                                                energy, 
                                                                self.non_planar_rings, 
                                                                self.flippable_Ns, 
                                                                self.flippable_Cs,
                                                                self.sulfo_matches, 
                                                                conf_ring_descriptors_df)
            return conf_ring_descriptors_df
        
        
        # Step 1: Generate initial conformers
        try:
            conf_ring_descriptors_df = embed_fix_ring_confs(method = 'srETKDGv3')
        except Exception as e: 
            print(e)

        if len(conf_ring_descriptors_df) == 0:
            # In case where srETKDGv3 failed in embedding the molecule, 
            # we have to use the macrocyclic version.
            print(f"srETKDGv3 failed for {self.name}, using macrocyclic version")
            logger.warning(f"srETKDGv3 failed for {self.name}, using macrocyclic version")
            conf_ring_descriptors_df = embed_fix_ring_confs(method = 'ETKDGv3')

        conf_ring_descriptors_df.sort_values(['equatorial_subs_Ns', 'equatorial_subs_Cs', 'Energy'],
                                            ascending=[False, False, True], inplace=True) 
        
        # Unlikely to have duplicate energy, but may happen for very small symmetric molecucles
        conf_ring_descriptors_df['Round_energy'] = conf_ring_descriptors_df['Energy'].round(4)
        conf_ring_descriptors_df.drop_duplicates(subset=['Round_energy'], keep='first', inplace=True) 

        # Keep a reservoir as the lowest energy possible conformer in case no good ring conformers are found.
        reservoir = conf_ring_descriptors_df.iloc[0, 0]
        # Keep the top 10 conformers for further processing (e.g., AMSOL)
        for idx in range(min(10, len(conf_ring_descriptors_df))):
            self.amsol_mol.AddConformer(conf_ring_descriptors_df.iloc[idx, 0], assignId=True)
        if self.VERBOSE:
            print(f"\tamsol_mol contains: {self.amsol_mol.GetNumConformers()}")
            conf_ring_descriptors_df.to_csv(f"{self.name}_confs.csv", index=False)
            temp_mol = Chem.Mol(self.empty_mol)
            with Chem.SDWriter(f"{self.name}_confs.sdf") as w:
                for _, row in conf_ring_descriptors_df.iterrows():
                    conf_idx = temp_mol.AddConformer(row[0], assignId=True)
                    w.write(temp_mol, confId=conf_idx)
        # Remove the conformers that do not compromise all the non-planar rings       
        #conf_ring_descriptors_df.to_csv(f"{self.name}_conf_ring_descriptors.csv", index=False)
        conf_ring_descriptors_df = utils.remove_unfavorable_confs(conf_ring_descriptors_df, self.name)

        if len(conf_ring_descriptors_df) == 0:
            # No conformer is found to compromise all the non-planar rings. 
            # Use the lowest energy conformer``
            logger.warning(f"No conformer is found to compromise all the non-planar rings for {self.name}. Using the lowest energy conformer.")
            scaffold = Chem.Mol(self.empty_mol)
            scaffold.AddConformer(reservoir, assignId=True)
            self.ring_confs.append(scaffold)
            mol2_obj = mol2writer.Mol2Writer(Chem.Mol(scaffold, confId = 0))
            self.mol2_str = mol2_obj.write_mol2()
            return

        # Process rigid scaffolds based on sulfo descriptors
        #conf_ring_descriptors_df.to_csv(f'{self.name}_confs.csv')
        if self.sulfo_matches:
            align_on = self.sulfo_matches[0]
            for sulfo_match in conf_ring_descriptors_df['sulfo_descriptors'].unique():
                temp_list = conf_ring_descriptors_df[conf_ring_descriptors_df['sulfo_descriptors'] == sulfo_match].values.tolist()
                initial_len = len(temp_list)
                num_confs_per_regioisomers = 0
                while num_confs_per_regioisomers < self.num_ring_confs and temp_list:
                    lowest_energy_entry = temp_list.pop(0)
                    conformer, current_descriptors = lowest_energy_entry[0], lowest_energy_entry[2:-1]
                    scaffold = Chem.Mol(self.empty_mol)
                    conf_id = scaffold.AddConformer(conformer, assignId=True)
                    # Need to align briefly so that coordinates are not too far apart and disrupt Mol2DB2
                    if self.ring_confs: rdMolAlign.AlignMol(scaffold, self.ring_confs[0], 0, 0, atomMap=[(i, i) for i in align_on])
                    self.ring_confs.append(Chem.Mol(scaffold, conf_id))
                    num_confs_per_regioisomers += 1
                    temp_list = utils.ring_conf_clusters(current_descriptors, temp_list)
                if self.VERBOSE: print(f'\tBefore: {initial_len}, after: {num_confs_per_regioisomers}')


        else:
            align_on = list(self.planar_rings)[0] if self.planar_rings else (1, 2, 3)
            temp_list = conf_ring_descriptors_df.values.tolist()
            while len(self.ring_confs) < self.num_ring_confs and temp_list:
                lowest_energy_entry = temp_list.pop(0)
                conformer, current_descriptors = lowest_energy_entry[0], lowest_energy_entry[2:-1]
                scaffold = Chem.Mol(self.empty_mol)
                scaffold.AddConformer(conformer, assignId=True)
                # Need to align briefly so that coordinates are not too far apart and disrupt Mol2DB2
                if self.ring_confs: rdMolAlign.AlignMol(scaffold, self.ring_confs[0], 0, 0, atomMap=[(i, i) for i in align_on])
                self.ring_confs.append(scaffold)
                temp_list = utils.ring_conf_clusters(current_descriptors, temp_list)
            if self.VERBOSE: print(f'\tBefore: {len(self.mol_H.GetConformers())}, after: {len(self.ring_confs)}')
        # Align the AMSOL coordinates to the ring conformations so that the coordinates are not too far apart
        for conf_id in range(self.amsol_mol.GetNumConformers()):
            rdMolAlign.AlignMol(self.amsol_mol, self.ring_confs[0], conf_id, 0, atomMap=[(i, i) for i in align_on])

        mol2_obj = mol2writer.Mol2Writer(Chem.Mol(self.amsol_mol, confId = 0))
        self.mol2_str = mol2_obj.write_mol2()

    def embed_smiles_babel(self):
        from openbabel import openbabel as ob
        # Step 1: Generate 3D conformer in Open Babel
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats("smi", "sdf")
        
        mol = ob.OBMol()
        obConversion.ReadString(mol, self.smiles)
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
        mol_rdkit.SetProp("_Name", self.name)

        conjugated_substituted_nitrogen_5aro = utils.find_conjugated_substituted_nitrogen_5aro(mol_rdkit)
        conjugated_substituted_nitrogen_6aro = utils.find_conjugated_substituted_nitrogen_6aro(mol_rdkit)

        amide_linkages = utils.find_amide(mol_rdkit)
        barbiturate_matches = utils.find_barbiturates(mol_rdkit)
        hydantoin_matches = utils.find_hydantoins(mol_rdkit)
        substituted_N_barbi_hydan_like = utils.find_substituted_N_barbi_hydan_like(mol_rdkit, barbiturate_matches, hydantoin_matches)
        self.mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_rdkit, mmffVariant="MMFF94s")
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_rdkit, self.mp, confId=0)
        if conjugated_substituted_nitrogen_5aro:
            for a, b, c, d, e, f in conjugated_substituted_nitrogen_5aro:
                ff.MMFFAddTorsionConstraint(a, b, c, d, False, 180, 180, 1)
                ff.MMFFAddTorsionConstraint(a, b, f, e, False, 180, 180, 1)
                ff.MMFFAddTorsionConstraint(b, c, d, e, False, 0, 0, 5)
                ff.MMFFAddTorsionConstraint(d, e, f, b, False, 0, 0, 5)
        if conjugated_substituted_nitrogen_6aro:
            for a, b, c, d, e, f, g in conjugated_substituted_nitrogen_6aro:
                ff.MMFFAddTorsionConstraint(a, b, c, d, False, 180, 180, 1)
                ff.MMFFAddTorsionConstraint(a, b, g, f, False, 180, 180, 1)
                ff.MMFFAddTorsionConstraint(b, c, d, e, False, 0, 0, 5)
                ff.MMFFAddTorsionConstraint(e, f, g, b, False, 0, 0, 5)
        if barbiturate_matches:
            for match in barbiturate_matches:
                n = len(match)
                for i in range(n):
                    a, b, c, d = [match[(i + j) % n] for j in range(4)]
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5) 
        if hydantoin_matches:
            for match in hydantoin_matches:
                n = len(match)
                for i in range(n):
                    a, b, c, d  = [match[(i + j) % n] for j in range(4)]
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 5)  
        if substituted_N_barbi_hydan_like:
            for a, b, c, d in substituted_N_barbi_hydan_like:
                ff.MMFFAddTorsionConstraint(a, b, c, d, False, 180, 180, 1)
        if amide_linkages:
            for a, b, c, d, e in amide_linkages: # O=C-N(-C)-H should be coplanar.
                ff.MMFFAddTorsionConstraint(a, b, c, d, False, 0, 0, 1)
                ff.MMFFAddTorsionConstraint(a, b, c, e, False, 180, 180, 1)
        ff.Minimize()
        self.ring_confs = [Chem.Mol(mol_rdkit)] # Replicate the output from embed_rdkit
        self.amsol_mol = Chem.Mol(mol_rdkit) # An RDKit Mol Object with upto 10 confs for AMSOL
        mol2_obj = mol2writer.Mol2Writer(mol_rdkit)
        self.mol2_str = mol2_obj.write_mol2()

    def embed_smiles_corina(self):
        '''
        Embed the SMILES string using CORINA and return the mol, net_charge,
        rigid_scaffolds, and flexible_scaffolds
        '''
        self.mol2_str, self.ring_confs = utils.embed_smiles_corina(self.smiles, self.name, self.num_ring_confs, self.VERBOSE)
        self.mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(self.ring_confs[0], mmffVariant="MMFF94s")
        self.amsol_mol = Chem.Mol(self.ring_confs[0]) # An RDKit Mol Object with upto 10 confs for AMSOL
        self.sulfo_matches = [] # No sulfonamide flipping in CORINA

    # ========== Torsional sampling =========================
    def stochastic_sampling(self, mol, tolerance_level, match_torlib, numConfs,
                        total_possible_solutions, window = 25, max_attempts=15000,
                        product=list(), unvisited = None, visited=None):
        """
        
        """
        if product: min_energy = min([conf[1] for conf in product])
        else: min_energy = 1e6

        attempts = 0

        bonded_pairs, same_parent_pairs = utils.precompute_bonded_and_same_parent_pairs(mol)
        # Condition to switch between visited matrix and unvisited set approaches
        if self.mode == 'extensive':
                n_transform = len(match_torlib)  # Number of rotatable bonds
                visitting = [0 for _ in range(n_transform)]
                if visited is None: visited = np.empty((0, n_transform))

                # New approach: use angles
                while len(product) < numConfs:
                    for idx in range(n_transform):
                        bond_idx = random.randint(0, len(match_torlib) - 1)
                        bond = match_torlib[bond_idx]
                        peaks = bond[2]  # Extract peaks
                        peak_idx = random.choices(range(len(peaks)), weights=[peak[3] for peak in peaks], k=1)[0]
                        peak = peaks[peak_idx]
                        random_angle = utils.get_random_angle(peak[0], peak[tolerance_level], method = 'uniform')
                        visitting[bond_idx] = random_angle
                        rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *bond[1], value=random_angle)

                    if utils.is_similar_conformer(np.array(visitting), visited, tol = self.tolerance) or \
                        utils.check_too_close_nonbonded_atoms(mol.GetConformer(0), mol, bonded_pairs, same_parent_pairs, threshold = self.threshold):
                        attempts += 1
                        if attempts > max_attempts:
                            break
                        continue
                    visited = np.vstack((visited, visitting))
                    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, self.mp, confId=0)
                    energy = ff.CalcEnergy()
                    if energy < min_energy: min_energy = min(energy, min_energy)
                    if energy <= min_energy + window: product.append((Chem.Conformer(mol.GetConformer(0)), energy))
        else:
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
                        rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *bond[1], value=utils.get_random_angle(peak[0], peak[tolerance_level]))

                    # Check if the conformation is valid and not already visited
                    if tuple(visitting) in visited or utils.check_too_close_nonbonded_atoms(mol.GetConformer(0), mol, bonded_pairs, same_parent_pairs, threshold = self.threshold):
                        attempts += 1
                        if attempts > max_attempts:
                            break
                        continue
                    
                    visited.add(tuple(visitting.copy()))
                    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, self.mp, confId=0)
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
                        value = utils.get_random_angle(peak[0], peak[tolerance_level])
                        rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *bond, value=value)

                    # If atoms are too close or if we already visited this conformation
                    if utils.check_too_close_nonbonded_atoms(mol.GetConformer(0), mol, bonded_pairs, same_parent_pairs, threshold = self.threshold):
                        attempts += 1
                        if attempts > max_attempts:
                            break
                        continue

                    unvisited.remove(choice)  # Remove the chosen combination from the unvisited set
                    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, self.mp, confId=0)
                    energy = ff.CalcEnergy()
                    if energy < min_energy: min_energy = energy
                    if energy <= min_energy + window: product.append((Chem.Conformer(mol.GetConformer(0)), energy))
                    #product.append((Chem.Conformer(mol.GetConformer(0)), energy))
        return product, visited, unvisited

    def conf_sampling(self, numConfs=2000, energywindow = 25, ignoreTorlib=False, AllowNonRing=False, request_alignment=None):
        """
        Perceive the allowed dihedral angles and call stochastic sampling to generate conformers.
        Args:
            numConfs (int): Number of conformers to generate.
            energywindow (float): Energy window for conformer generation.
            ignoreTorlib (bool): Whether to ignore the torsion library.
            AllowNonRing (bool): Whether to allow the full sampling of non-ring compounds.
            request_alignment (list): List of atom indices for alignment.
        
        Notes:
            For assymetric sulfonamides, there would be two versions of rigid scaffolds handled by EirVS.
        """
        self.conf_sampled = True
        if self.request_alignment is None and request_alignment is not None:
            self.request_alignment = request_alignment
        
        if self.mode == 'extensive2':
            # This is still experimental, call conf_samplingv2, where Torlib is read differently
            # and the angle chosen would be deterministic. By default, peak +- 30 degrees, if tol2 >= 30.
            self.conf_samplingv2(numConfs = numConfs, energywindow = energywindow, AllowNonRing=False, request_alignment=request_alignment)
            return
        
        num_confs_by_rotbonds, match_torlib = utils.count_confs_by_rotbonds(self.ring_confs[0], ignoreTorlib, self.VERBOSE)
        requested_num_confs = numConfs

        #if VERBOSE: print(f"\t{num_confs_by_rotbonds} {num_confs_H} {num_rotatable_H} {numConfs}")
        if self.VERBOSE: print(f"\tTheory: {num_confs_by_rotbonds} possible conformations")

        # Find the rigid part only once outside the loop to save processing time
        self.atom_maps, self.label_map = utils.find_rigid_part(self.ring_confs[0], request_alignment)
        # Molecules which don't have rings are not of interest --> only sample limitedly.
        if (self.label_map) and not (AllowNonRing): numConfs = 30

        if request_alignment and not self.atom_maps:
            utils.log_error(self.smiles, self.name)
            return
        # For very flexible molecules, we need to sample more, then filter by energy later
        else:
            if numConfs*10 < num_confs_by_rotbonds: numConfs = min(int(numConfs * 1.5), num_confs_by_rotbonds)
            elif numConfs*5 < num_confs_by_rotbonds: numConfs = min(int(numConfs * 1.25), num_confs_by_rotbonds)
            else: numConfs = numConfs#min(numConfs, num_confs_by_rotbonds)


        if self.VERBOSE:
            rigid_info = f'\tFound {self.atom_maps} ({self.label_map}) as a rigid part' if self.label_map else f'\tFound {self.atom_maps} (rings) as rigid parts'
            print(rigid_info)


        if (len(match_torlib) == 0):
            print('No rotatable bonds found, returning the original conformation')
            return

        for idx, mol in enumerate(self.ring_confs):
            if self.VERBOSE: print(f"\tHandling ring/sulfonamide conformation {idx+1}/{len(self.ring_confs)}")
            original_mol = Chem.Mol(mol)
            processing_mol = Chem.Mol(mol)
              
            # Only remap the match_torlib when sulfo_matches is found
            if self.sulfo_matches: num_confs_by_rotbonds, match_torlib = utils.count_confs_by_rotbonds(mol, ignoreTorlib, self.VERBOSE)
            if self.VERBOSE: print('\tRunning stochastic torsional sampling')
            
            product, visited, unvisited = self.stochastic_sampling(processing_mol, 1, match_torlib, numConfs, num_confs_by_rotbonds, energywindow, 15000, list(), visited = None, unvisited=None)
            if len(product) <= min(numConfs, num_confs_by_rotbonds) // 3: 
                if self.VERBOSE: print(f'Failed for stochastic sampling (generated {len(product)} confs), use the 2nd tolerance level')
                product, visited, unvisited = self.stochastic_sampling(processing_mol, 2, match_torlib, numConfs, num_confs_by_rotbonds, energywindow, 15000, product, visited = visited, unvisited = unvisited)

            if len(product) == 0: # No conformers are generated, use initial conformation instead
                print(f'Failed for stochastic sampling (generated {len(product)} confs), use the original conformation')
                continue
            
            product.sort(key=lambda x: x[1]) #Sort by energy
            before_energy = len(product)
            min_energy = product[0][1]
            filtered_product = [x[0] for x in product if x[1] - min_energy <= energywindow]
            filtered_product = filtered_product[:requested_num_confs]

            mol.RemoveAllConformers()
            largest_ring = max(self.atom_maps, key=len)
            for conf in filtered_product:
                confId = mol.AddConformer(conf, assignId=True)
                rdMolAlign.AlignMol(mol, original_mol, confId, 0, atomMap=[(i, i) for i in largest_ring])

            if self.VERBOSE: print(f"\tEnergy filter: {before_energy} -> {len(filtered_product)}")

    # Deterministic sampling
    # This is still experimental, call conf_samplingv2, where Torlib is read differently
    # and the angle chosen would be deterministic. By default, peak +- 30 degrees, if tol2 >= 30.
    def stochastic_sampling_v2(self, mol, angle_map, score_map, numConfs, possible_numConfs, importance_order,
                        window = 25, max_attempts=50_000, product=list()):
        """"""
        bonded_pairs, same_parent_pairs = utils.precompute_bonded_and_same_parent_pairs(mol)
        attempts = 0
        min_energy = 1e6

        
        if possible_numConfs <= max_attempts: # Try 
            # Generate all combinations, then randomly taken from them, only valid for small combinatorial space
            # If the number of conformations is manageable, we can enumerate all combinations
            
            # Generate all possible combinations of dihedral angles
            unvisited = list(itertools.product(*(angle_map[bond_idx][2] for bond_idx in range(len(angle_map)))))
            if self.VERBOSE: 
                print(f"\tEnumerated all possible combinations of dihedral angles")
                print(f"\tTotal number of unique conformations: {len(unvisited)}")
            random.shuffle(unvisited)

            while len(product) < numConfs and len(unvisited) > 0:
                choice = unvisited.pop()  # Randomly select a combination of dihedral angles

                # Set the dihedrals based on the chosen combination
                for idx, val in enumerate(angle_map.values()): #Iterate through the rotatable bonds
                    dihedral_atoms = val[1]  # Get the atom indices for the dihedral
                    #print(dihedral_atoms)
                    angle = choice[idx]  # Get the angle for this dihedral from the chosen combination
                    # bond_idx corresponds to the index of the bond in the list of rotatable bonds
                    # Set the dihedral angle for the corresponding bond
                    rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *dihedral_atoms, angle)
                # If atoms are too close or if we already visited this conformation
                if utils.check_too_close_nonbonded_atoms(mol.GetConformer(0), mol, bonded_pairs, same_parent_pairs, threshold = self.threshold):
                    continue

                ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, self.mp, confId=0)
                energy = ff.CalcEnergy()
                if energy < min_energy: min_energy = energy
                if energy <= min_energy + window: product.append((Chem.Conformer(mol.GetConformer(0)), energy))
        else:
            # Reweight the importance of the bonds
            # If the number of conformations is too large, we can randomly sample
            if self.VERBOSE:
                print("\tStochastic sampling with importance-based weights")
                print("\tImportance order of bonds:", importance_order)
                print(f"\tNumber of conformations {possible_numConfs}. Random sampling will be performed.")
            
            # Initialize tracking variables for angles
            k = len(angle_map)
            visited = set()
            visitting = [0] * len(angle_map)

            # Adaptive sampling parameters
            max_stagnation = min(max_attempts // 10, 5000)  # Stop if no progress
            stagnation_counter = 0
            last_product_size = 0

            # First, reset all the dihedral to a default angle_peak 0.
            for bond_idx, rule in enumerate(angle_map.values()):
                visitting[bond_idx] = rule[2][0]
                rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *rule[1], rule[2][0])

            while len(product) < numConfs and attempts < max_attempts:
                # Use importance-based weights for rotation selection
                to_rotate = set(random.choices(range(len(angle_map)), weights=importance_order, k=k))
                # For each selected bond, choose a random angle
                for bond_idx in to_rotate:
                    if bond_idx >= len(angle_map):
                        continue
                        
                    # Get the atom indices and possible angles
                    dihedral_atoms = angle_map[bond_idx][1]
                    possible_angles = angle_map[bond_idx][2]
                    angle_scores = score_map[bond_idx]
                    angle_idx = random.choices(range(len(possible_angles)), 
                                                    weights=angle_scores, k=1)[0]
                            
                    angle = possible_angles[angle_idx]
                    visitting[bond_idx] = angle
                    # Set the dihedral angle
                    rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), *dihedral_atoms, angle)
                
                state_tuple = tuple(visitting)    
                
                if state_tuple in visited:
                    attempts += 1
                    continue
                
                if utils.check_too_close_nonbonded_atoms(mol.GetConformer(0), mol, bonded_pairs, same_parent_pairs, threshold = self.threshold):
                    attempts += 1
                    visited.add(state_tuple)
                    continue

                visited.add(state_tuple)
                ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, self.mp, confId=0)
                energy = ff.CalcEnergy()
                if energy < min_energy: min_energy = energy
                if energy <= min_energy + window: product.append((Chem.Conformer(mol.GetConformer(0)), energy))

                # Check for early stopping conditions
                if len(product) == last_product_size:
                    stagnation_counter += 1
                    if stagnation_counter >= max_stagnation:
                        if self.VERBOSE:
                            print(f"Stopping due to stagnation after {attempts} attempts. Generated {len(product)} conformers.")
                        break
                else:
                    stagnation_counter = 0
                    last_product_size = len(product)
                
        return product
    
    def conf_samplingv2(self, numConfs=2000, energywindow = 25, AllowNonRing=False, request_alignment=None):
        rot_bonds = utils.getDihedralMatches_v2(self.ring_confs[0])
        possible_numConfs, angle_map, score_map = utils.count_confs_by_rotbonds_v2(self.ring_confs[0], rot_bonds)#, VERBOSE=True)
        importance_order = utils.get_importance_order(self.ring_confs[0], rot_bonds)
        requested_num_confs = numConfs

        if self.VERBOSE:
            print(f"\tTheory: {possible_numConfs} possible conformations")
            print(f"\tRotatable bonds: {rot_bonds}")
            print(f"\tPossible angles:")
            for idx in range(len(angle_map)):
                print(f"\t{angle_map[idx]}")
                print(f"\t{score_map[idx]}")

        # Find the rigid part only once outside the loop to save processing time
        self.atom_maps, self.label_map = utils.find_rigid_part(self.ring_confs[0], request_alignment)



        if request_alignment and not self.atom_maps:
            utils.log_error(self.smiles, self.name)
            return
        # For very flexible molecules, we need to sample more, then filter by energy later
        else:
            if numConfs*10 < possible_numConfs: numConfs = min(int(numConfs * 1.5), possible_numConfs, requested_num_confs + 1000)
            elif numConfs*5 < possible_numConfs: numConfs = min(int(numConfs * 1.25), possible_numConfs, requested_num_confs + 1000)
            else: numConfs = numConfs#min(numConfs, num_confs_by_rotbonds)

        # Molecules which don't have rings are not of interest --> only sample limitedly.
        if (self.label_map) and not (AllowNonRing): numConfs = 30

        if (len(rot_bonds) == 0):
            print('No rotatable bonds found, returning the original conformation')
            return
        
        for idx, mol in enumerate(self.ring_confs):
            if self.VERBOSE: print(f"\tHandling ring/sulfonamide conformation {idx+1}/{len(self.ring_confs)}")
            original_mol = Chem.Mol(mol)
            processing_mol = Chem.Mol(mol)
              
            # Only remap the match_torlib when sulfo_matches is found
            if self.sulfo_matches: possible_numConfs, angle_map, score_map = utils.count_confs_by_rotbonds_v2(mol, rot_bonds)
            if self.VERBOSE: print('\tRunning stochastic torsional sampling')
            
            product = self.stochastic_sampling_v2(processing_mol, angle_map, score_map, numConfs, possible_numConfs, importance_order, energywindow, 50_000, list())
            
            if len(product) == 0:
                print(f'Failed to find any confs (generated {len(product)} confs), using the random dihedral angles approach as a fallback')
                num_confs_by_rotbonds, match_torlib = utils.count_confs_by_rotbonds(mol = self.ring_confs[0], VERBOSE = self.VERBOSE)
                product, _, _ = self.stochastic_sampling(processing_mol, 2, match_torlib, numConfs, num_confs_by_rotbonds, energywindow, 15000, list(), visited = None, unvisited=None)
                
                if len(product) == 0: 
                    print(f'Failed for stochastic sampling (generated {len(product)} confs), use the original conformation')
                    continue

            product.sort(key=lambda x: x[1]) #Sort by energy
            before_energy = len(product)
            min_energy = product[0][1]
            filtered_product = [x[0] for x in product if x[1] - min_energy <= energywindow]
            filtered_product = filtered_product[:requested_num_confs]
            mol.RemoveAllConformers()
            largest_ring = max(self.atom_maps, key=len)
            for conf in filtered_product:
                confId = mol.AddConformer(conf, assignId=True)
                rdMolAlign.AlignMol(mol, original_mol, confId, 0, atomMap=[(i, i) for i in largest_ring])
            if self.VERBOSE: print(f"\tEnergy filter: {before_energy} -> {len(filtered_product)}")



    # ========== Output to different file formats ===========

    def to_sdf(self, filename = None):
        """
        Write the conformers to an SDF file.
        Args:
            filename (str): The name of the output SDF file. If None, defaults to self.name.sdf.
        """
        if self.conf_sampled == False:
            raise ValueError("Conformers have not been generated yet. Please call conf_sampling() first.")
        
        if filename is None:
            filename = self.name

        if len(self.ring_confs) == 1:
            with Chem.SDWriter(f"sdf/{filename}.sdf") as writer:
                for confid in range(self.ring_confs[0].GetNumConformers()):
                    writer.write(self.ring_confs[0], confId=confid)
        else:
            for idx, ring_conf in enumerate((self.ring_confs)):
                with Chem.SDWriter(f"sdf/{filename}.nr{idx}.sdf") as writer:
                    for confid in range(ring_conf.GetNumConformers()):
                        writer.write(ring_conf, confId=confid)

    def to_mol2(self, filename = None):
        """
        Write the conformers to an Mol2 file.
        Args:
            filename (str): The name of the output Mol2 file. If None, defaults to self.name.sdf.
        """
        if self.conf_sampled == False:
            raise ValueError("Conformers have not been generated yet. Please call conf_sampling() first.")
        
        if filename is None:
            filename = self.name

        if len(self.ring_confs) == 1:
            mol2_obj = mol2writer.Mol2Writer(self.ring_confs[0], mol2_template=self.mol2_str)
            mol2_obj.write_mol2(filename=f"mol2/{filename}.mol2")
        else:
            for idx, ring_conf in enumerate(self.ring_confs, mol2_template=self.mol2_str):
                mol2_obj = mol2writer.Mol2Writer(ring_conf)
                mol2_obj.write_mol2(filename=f"mol2/{filename}.nr{idx}.mol2")

    def to_pdbqt(self, filename = None):
        """
        Write the conformers to a PDBQT file using Meeko.
            
        Args:
            filename (str): The name of the output PDBQT file. If None, defaults to self.name.pdbqt.
        """
        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy
        except ImportError:
            print("""Please install the meeko package using "pip install meeko" to use this script.""")
            exit(1)
        pass

        if filename is None:
            filename = self.name
            
        mkprep = MoleculePreparation()
        if not self.ring_confs:
            logger.warning(f"Failed to embed {self.name}.")
            return
        is_multi = len(self.ring_confs) > 1
        for i, mol in enumerate(self.ring_confs):
            prepared_mol = mkprep(Chem.Mol(mol, confId=0))
            pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(prepared_mol[0])
            if success:
                #print(pdbqt_string)
                if is_multi:
                    with open(f"pdbqt/{filename}.nr{i}.pdbqt", 'w') as f:
                        for line in pdbqt_string:
                            f.write(line)
                else:
                    with open(f"pdbqt/{filename}.pdbqt", 'w') as f:
                        for line in pdbqt_string:
                            f.write(line)
            else:
                print(error_msg)

    def to_db2(self, 
               numConfs = 2000,
               energywindow = 25,
               ignoreTorlib = False,
               AllowNonRing = False,
               request_alignment = None,
               longname = "fake",
               env = None,
               cleanup = True,
               tarfile = None
               ):
        """
        Convert the conformers to DB2 format and save them to a file.
        If the ConformerGenerator object has undergone conformational sampling (confgen.conf_sampled == True),
        the arguments numConfs, energywindow, ignoreTorlib, AllowNonRing, and request_alignment are ignored.

        Args:
            numConfs (int): Number of conformers to generate.
            energywindow (float): Energy window for conformer generation.
            ignoreTorlib (bool): Whether to ignore the torsion library.
            AllowNonRing (bool): Whether to allow the full sampling of non-ring compounds.
            request_alignment (list): List of atom indices for alignment.
            longname (str): Long name for the molecule.
            env (str): Environment for AMSOL.
            cleanup (bool): Whether to clean up the temporary files.
            tarfile (tarball object): The tarball object to write the DB2 data to.
        """
        
        if self.VERBOSE: print("Solvating...")
        if self.request_alignment is None and request_alignment is not None:
            self.request_alignment = request_alignment
        else: request_alignment = self.request_alignment # Need this so that if not determined, use the class attributes.
        ### Solvation ###
        if not(env): env = setup_env() 
        os.makedirs(f"solv/{self.name}", exist_ok=True)
        os.chdir(f"solv/{self.name}")
        for conf_id in range(self.amsol_mol.GetNumConformers()):
            try:
                # Idea: try from the energy minimum conformer if AMSOL fails -> next conformer until reach the last
                if self.VERBOSE: print(f"\tTrying conformer: {conf_id}")
                error_signal = 0
                if self.method == 'rdkit':
                    cp = Chem.Mol(self.amsol_mol, confId=conf_id) #Retrieve the conf_id-th conformer of mol object
                    mol2_obj = mol2writer.Mol2Writer(cp, mol2_template = self.mol2_str, atom_attributes = True)
                    mol2_obj.write_mol2(f"{self.name}.mol2")
                else:
                    # Babel and CORINA, use the mol2_str as only 1 conformer is needed
                    write_to_file(self.mol2_str, f"{self.name}.mol2")

                run_amsol.prepare(f"{self.name}.mol2", self.name, self.netcharge)
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
        if error_signal == -1 and conf_id + 1 == self.amsol_mol.GetNumConformers(): # AMSOL failed
            logger.error(f"AMSOL failed for {self.name}, skipping it")
            utils.log_error(self.smiles, self.name)
            try: # Clean up the folders if error occurs. This help to not overfill the disk
                shutil.rmtree(f"solv/{self.name}", ignore_errors=True)
            except: pass
            return
        self.amsol_time = time.time()
        shutil.copy(f"solv/{self.name}/output.mol2", f"solv/{self.name}/{self.name}_solv.mol2")
        shutil.move(f"solv/{self.name}/output.solv", f"solv/{self.name}/{self.name}_solv.solv")

        ### Torsional sampling ###
        if self.VERBOSE: print("Torsional sampling...")

        if not(self.conf_sampled):
            self.conf_sampling(numConfs=numConfs,
                               energywindow = energywindow,
                               ignoreTorlib = ignoreTorlib,
                               AllowNonRing = AllowNonRing,
                               request_alignment = request_alignment)

        ### Output to DB2 ###
        if self.VERBOSE: print("Output to DB2...")
        os.makedirs(f"db2/{self.name}", exist_ok=True)
        try:
            shutil.move(os.path.join("solv", self.name, f"{self.name}_solv.solv"), os.path.join("db2", self.name, f"{self.name}.solv"))
            shutil.move(os.path.join("solv", self.name, f"{self.name}_solv.mol2"), os.path.join("db2", self.name, f"{self.name}.mol2"))

            os.chdir(f"db2/{self.name}")
            db2_data_all = ""
            solv_obj = solv.Solv(f"{self.name}.solv")
            for ring_conf in self.ring_confs:
                for rigid_scaffold in self.atom_maps:
                    db2_data = utils.Align_ConvertToDb2(ring_conf, rigid_scaffold, solv_obj, self.name, self.smiles, longname) 
                    db2_data_all += db2_data
            if tarfile: write_to_tarball(tarfile, db2_data_all.encode('utf-8'), name=f"{self.name}.db2")
            else: write_to_file(db2_data_all, f"../{self.name}.db2")
            os.chdir("../..")
            if not self.VERBOSE: utils.remove_folders([f"solv/{self.name}"])
            if cleanup:
                utils.remove_folders([f"db2/{self.name}"])
                
        except Exception as e:
            logger.error(f"Error in converting {self.name} to DB2 format: {e}")
            os.chdir("../..")
            try: # Clean up the folders if error occurs. This help to not overfill the disk
                shutil.rmtree(f"solv/{self.name}", ignore_errors=True)
                shutil.rmtree(f"db2/{self.name}", ignore_errors=True)
            except: pass
            utils.log_error(self.smiles, self.name)
            return


def initial_embedding(queue, smiles, name, randomSeed, nr = 1, numcores = 1, VERBOSE = False):
    """
    A wrapper for multiprocessing to call so that timeout works.
    Embed the SMILES string using RDKit, then return the 3D coordinates in the binary format.
    The coordinates will be used for recovery of the confgen object as many of the other class variables are not picklable.
    Args:
        queue (multiprocessing.Queue): The queue to put the results into.
        smiles (str): The SMILES string to embed.
        name (str): The name of the molecule.
        randomSeed (int): The random seed for embedding.
        nr (int): The number of ring conformers to generate.
        numcores (int): The number of CPU cores to use for parallel processing.
        VERBOSE (bool): If True, print verbose output.
    returns:
        bin_amsol_mol (bytes): The binary representation of the AMSOL molecule.
        bin_conf_rings (list): A list of binary representations of different ring conformers.
        mol2_str (str): The mol2 block of the molecule.
    """
    if VERBOSE:
        print("Generating initial 3D conformations...")
    try:
        confgen = ConformerGenerator(smiles, 
                                     name, 
                                     numcores=numcores, 
                                     num_ring_confs=nr, 
                                     randomSeed=randomSeed,
                                     method='rdkit', 
                                     VERBOSE=VERBOSE)
        property_flags = (
                PropertyPickleOptions.MolProps |
                PropertyPickleOptions.PrivateProps
                )
        bin_amsol_mol = confgen.amsol_mol.ToBinary(propertyFlags=property_flags)
        bin_conf_rings = [ringconf.ToBinary(propertyFlags=property_flags) for ringconf in confgen.ring_confs]
        mol2_str = confgen.mol2_str
        queue.put((bin_amsol_mol, bin_conf_rings, mol2_str,  None))
    except Exception as e:
        print(e)
        queue.put((None, None, None,  str(e)))

def setup_env():
    '''
    Set up the environment variables for AMSOL libraries.
    '''
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

def gen_conf_chunk(df: pd.DataFrame, args, input_file='0'):
    """
    Generate conformers for a given DataFrame of SMILES strings and save them in different formats.
    Args:
        df (pd.DataFrame): DataFrame containing SMILES strings and other relevant information.
        args (Namespace): 
            Parsed arguments containing various configuration options, including:
            - randomSeed (int): Seed for random number generation.
            - numconfs (int): Number of conformations to generate.
            - debug (bool): Verbose output for debugging.
            - cleanup (bool): Whether to remove intermediate files after processing.
            - energywindow (float): Energy window for conformer sampling.
            - timeout (int): Timeout (in minutes) for RDKit-based conformation generation.
            - ignoretorlib (bool): Whether to ignore torsion library constraints.
            - rigid (str): SMILES/SMARTS for conformers to be aligned to.
            - nringconfs (int): Number of ring conformers to generate.
            - numcores (int): Number of CPU cores to use for parallel processing.
            - method (str): Method for initial conformation generation ('rdkit', 'obabel', 'corina').
            - timing (bool): If enabled, logs timing information for each step.
            - smiles (bool): If True, skips restarting logic.
            - out (list): List of output formats to generate (e.g., 'pdbqt', 'sdf', 'mol2', 'db2').

        input_file (str): Name of the input file (default is '0').
    """
    if 'mol' not in df.columns:
        df['mol'] = df['smiles'].apply(Chem.MolFromSmiles)
    df = filters.Filters.remove_exotic_chem_to_db2(df)
    randomSeed, numConfs, VERBOSE, cleanup, energywindow, timeout, request_alignment, nr, numcores, mode, tolerance = \
        args.randomSeed, args.numconfs, args.debug, args.cleanup, args.energywindow, args.timeout, args.rigid, args.nringconfs, args.numcores, args.mode, args.tolerance
    
    ignoreTorlib = (args.mode == 'ignoretorlib')
    request_alignment = Chem.MolFromSmarts(utils.canonicalize_if_smiles(request_alignment)) if request_alignment else None
    
    # Test mode in unittest, not to produce redundant files here
    if args.test: 
        os.chdir(args.prefix)
    processed_mols = set()
    os.makedirs(f"db2", exist_ok=True)
    if 'pdbqt' in args.format: os.makedirs(f"pdbqt", exist_ok=True)
    if 'sdf' in args.format: os.makedirs(f"sdf", exist_ok=True)
    if 'mol2' in args.format: os.makedirs(f"mol2", exist_ok=True)
    output_tgz = f"db2/{input_file}.db2.tgz"
    env = setup_env()

    if args.timing: 
        if not(os.path.exists('eirvs_timing.csv')): 
            with open('eirvs_timing.csv', 'w') as f: f.write('Name,Initial embedding,AMSOL,Torsional sampling,Mol2DB2,Total\n')
        logging_time = ""
    
    # Check if the output file already exists. A sign of unfinished job
    restart_flag = False
    if not(args.smiles) and ('db2.tgz' in args.format) and os.path.exists(output_tgz):
        
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
            random.seed(randomSeed)
            smiles = row['smiles']
            name = row['ids']
            longname = row['highlights'] if args.synthon else None
            if name in processed_mols:
                print(f"Skipping {name} as it already exists")
                logger.info(f"Skipping {name} as it already exists")
                continue
            logger.info(f"Handling {name}")
            if VERBOSE: print(f"Handling {name}")
            if args.timing: start = time.time() 
            try:
                if args.method == 'corina':
                    confgen = ConformerGenerator(smiles, name, method='corina', mode = mode, tolerance=tolerance, VERBOSE=VERBOSE)
                elif args.method == 'obabel':
                    confgen = ConformerGenerator(smiles, name, method='obabel', mode = mode, tolerance=tolerance, VERBOSE=VERBOSE)
                else:
                    queue = multiprocessing.Queue()
                    process = multiprocessing.Process(target=initial_embedding, args=(queue, smiles, name, randomSeed, nr, numcores, VERBOSE))
                    process.start()
                    process.join(timeout=timeout*60)  # default 2 minutes timeout
                    # Check if process is still alive (meaning it exceeded timeout)
                    if process.is_alive():
                        logger.warning(f"Timeout occurred while generating conformation for {name}, using OpenBabel.")
                        process.terminate()
                        process.join()
                        try:
                            confgen = ConformerGenerator(smiles, name, num_ring_confs=nr, method='babel', tolerance=tolerance, VERBOSE=VERBOSE)
                        except Exception as e:
                            logger.error(f"Error in generating initial conformation using OpenBabel for {name}, skipping it {e}")
                            utils.log_error(smiles, name)
                            continue
                        if confgen.amsol_mol is None:
                            logger.error(f"Error in generating initial conformation using OpenBabel for {name}, skipping it")
                            utils.log_error(smiles, name)
                            continue

                    # Retrieve result from queue
                    elif not queue.empty():
                        bin_amsol_mol, bin_conf_rings, mol2_str, error = queue.get() 

                        if error:
                            logger.error(f"Error in generating initial conformation using RDKit for {name}, skipping it: {error}")
                            utils.log_error(smiles, name)
                            continue
                        confgen = ConformerGenerator.from_existing_data(smiles, 
                                                                        name, 
                                                                        bin_amsol_mol, 
                                                                        bin_conf_rings, 
                                                                        mol2_str, 
                                                                        request_alignment, 
                                                                        mode,
                                                                        tolerance,
                                                                        VERBOSE)
                    else:
                        logger.error(f"Unknown error in generating initial conformation for {name}, skipping it.")
                        utils.log_error(smiles, name)
                        continue
            except Exception as e:
                logger.error(f"Error in generating initial conformation for {name}, skipping it: {e}")
                utils.log_error(smiles, name)
                continue

            if args.timing: embed_time = time.time() # Time for embedding
            
            if 'pdbqt' in args.format: confgen.to_pdbqt()

            if any(format in args.format for format in ['sdf', 'mol2', 'db2', 'db2.tgz']):
                confgen.conf_sampling(numConfs=numConfs,
                                      energywindow=energywindow,
                                      ignoreTorlib=ignoreTorlib,
                                      AllowNonRing=False,
                                      request_alignment=request_alignment,
                                     )
            
            if args.timing: sampling_time = time.time() # Time for sampling
            if 'sdf' in args.format: confgen.to_sdf()
            if 'mol2' in args.format: confgen.to_mol2()
            if 'db2' in args.format:
                try: 
                    confgen.to_db2(longname = longname, env=env, cleanup=cleanup)
                except Exception as e:
                    logger.error(f"Error in converting {name} to DB2 format: {e}")
                    utils.log_error(smiles, name)
                    continue
            if 'db2.tgz' in args.format:
                try: 
                    confgen.to_db2(longname = longname, env=env, cleanup=cleanup, tarfile = output)
                except Exception as e:
                    logger.error(f"Error in converting {name} to DB2 format: {e}")
                    utils.log_error(smiles, name)
                    continue
                
            if args.timing: 
                db2_time = time.time()
                logging_time += f'{name},{embed_time-start},{confgen.amsol_time-sampling_time},{sampling_time-embed_time},{db2_time-confgen.amsol_time},{db2_time-start}\n'

    # Use this method to remove the tarball if it is empty. 
    # The "with open" method is better to handle unexpected error that lead to corrupted files 
    if ('db2.tgz' not in args.format) and os.path.exists(f"db2/{input_file}.db2.tgz"):
        os.remove(f"db2/{input_file}.db2.tgz")
    
    if not(utils.is_slurm_job()):
        folders_to_remove = ['3d', 'solv'] if cleanup else ['solv']
        for folder in folders_to_remove:
            if os.path.exists(folder) and os.path.isdir(folder):
                try:
                    os.rmdir(folder)
                except: pass

    if ('db2.tgz' not in args.format) and ('db2' not in args.format):
        try:
            if os.path.exists('db2') and os.path.isdir('db2') and not os.listdir('db2'):
                os.rmdir('db2')
        except: pass

    if args.timing:
        with open('eirvs_timing.csv', 'a') as f:
            f.write(logging_time)

class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def _format_action_invocation(self, action):
        """
        Override to customize the argument display in the help message.
        Suppress the metavar formatting like `$short $metavar, $long=$metavar`.
        """
        if not action.option_strings:
            return super()._format_action_invocation(action)

        parts = []
        for option_string in action.option_strings:
            parts.append(option_string)
        return ', '.join(parts)
    
def main():
    parser = argparse.ArgumentParser(description="Generate conformers for a given SMILES string/file.\nTwo-column files are required.",
                                     formatter_class=CustomHelpFormatter,
                                     add_help = False)  # Suppress default -h/--help)
    parser.add_argument('--input_files', '-i', type=str, default = None, help='Input file containing SMILES strings.')
    parser.add_argument('--smiles', '-s', type = str, default = None, help='Input SMILES string.')
    parser.add_argument('--prefix', '-p', type=str, default = 'db2', help='Prefix for the output files.')
    parser.add_argument('--format', '-f', type=str, nargs='+', default=['db2.tgz'], choices=['db2', 'db2.tgz', 'pdbqt', 'sdf', 'mol2'], help='Output format(s) (e.g., db2.tgz, pdbqt, sdf, mol2).')
    parser.add_argument('--mode', '-mode', type=str, default='vs', choices=['vs', 'extensive', 'ignoretorlib'], help='Mode for conformer generation (vs, extensive, ignoretorlib).')
    parser.add_argument('--tolerance', '-tol', type=float, default=30, help='Tolerance for dihedral angle sampling (default: 30).')
    parser.add_argument('--numconfs', '-nconfs', type=int, default=2000, help='Number of conformers to generate (default: 2000).')
    parser.add_argument('--nringconfs', '-nr', type=int, default=1, help='Number of ring conformers to generate (default: 1).')
    parser.add_argument('--debug', '-d', action='store_true', help='Enable verbose output for debugging.')
    parser.add_argument('--timeout', '-to',type=int, default=2, help='Timeout in minutes for RDKit-based conformation generation.')
    parser.add_argument('--energywindow', '-w',type=float, default=25.0, help='Energy window for conformer generation.')
    parser.add_argument('--numcores', '-j', type=int, default=4, help='Number of CPU cores to use (default: 4).')
    parser.add_argument('--method', '-m',type=str, choices=['rdkit', 'obabel', 'corina'], default='rdkit', help='Method for initial conformation generation.')
    parser.add_argument('--timing', action='store_true', help='Log timing information for each step.')
    parser.add_argument('--nocleanup', action='store_false', dest='cleanup', help='Do not remove intermediate files after processing.')
    parser.add_argument('--rigid', '-r', type=str, default=None, help='SMILES/SMARTS for conformers to be aligned to.')
    parser.add_argument("--help", "-h", action="help", help="Show this help message and exit")
    parser.add_argument('--randomSeed', '-rs',type=int, default=42, help=argparse.SUPPRESS)
    parser.add_argument('--test', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--synthon', action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    
    if args.input_files and args.smiles:
        parser.error('Please provide either input files or SMILES strings, not both.')

    if args.input_files is not None:
        for inFile in args.input_files:
            if not Path(inFile).is_file():
                parser.error(f'The input file: {inFile} does not exist.')
        args.input_files = [Path(inFile).resolve() for inFile in args.input_files]
        for inFile in args.input_files:
            gen_conf_chunk(pd.read_csv(inFile, sep=' ', header=None, names=['smiles', 'ids']), args, inFile.stem)
    else:
        if args.smiles is not None:
            smiles = args.smiles
            df = pd.DataFrame({'smiles': [smiles], 'ids': ['0']})
            gen_conf_chunk(df, args)
        else:
            parser.error('Please provide either input files or SMILES strings.')
    # Call the main function
    #process_files(args)
   
if __name__ == "__main__":
    main()