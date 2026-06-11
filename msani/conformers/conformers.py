"""Based partly on: 
    https://github.com/UnixJunkie/smi2sdf3d/blob/master/smi2sdf.py
    http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
    http://pubs.acs.org/doi/abs/10.1021/ci2004658
    https://greglandrum.github.io/rdkit-blog/posts/2024-02-11-more-multithreading.html

    Should try to sample all possible conformations based on dihedral angles sampling based on: https://github.com/dkoes/rdkit-scripts/blob/master/rdallconf.py
"""
# Author: Thua-Phong Lam, Jens Carlsson lab, Uppsala University
# Date: 2025-10-07

import logging
import os
import multiprocessing
import subprocess
import shutil
import random
import sys
import tarfile, io
import time
import argparse

from pandas import DataFrame, read_csv  # only what you use
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolAlign, PropertyPickleOptions

from msani.io.parsers import CustomHelpFormatter
from msani.conformers import utils, mol2writer, torsions
from msani.filtering import filters
from msani.db2 import solv
from msani.io.utils import log_error

# Check if Open Babel is installed
try:
    from openbabel.openbabel import OBMol
    OBABEL_AVAILABLE = True
    obabel_path = os.path.join(os.path.dirname(sys.executable), 'obabel') # Ensure that the exact obabel within the same conda environment is used
except:
    OBABEL_AVAILABLE = False
    pass

# Check if AMSOL is correctly installed
try:
    from msani.amsol import run_amsol
    if run_amsol.AMSOLEXE: AMSOL_AVAILABLE = True
    else: AMSOL_AVAILABLE = False
except ImportError:
    AMSOL_AVAILABLE = False
    pass

# Check if Meeko is installed
try:
    from meeko.preparation import MoleculePreparation
    from meeko.writer import PDBQTWriterLegacy
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False

# Check if the CPP-accelerated sampling module is available
try:
    import msani_confgen_cpp as cpp_sampler
    CPP_AVAILABLE = True
except:
    CPP_AVAILABLE = False

logger = logging.getLogger('msani')

SRLib = torsions.SmallRingLibrary()
Torlib = torsions.TorsionLibrary()

class ConformerGenerator:
    '''
    Class to generate conformers from SMILES strings.

    Parameters
    ----------
    smiles : str
        The SMILES string of the molecule.
    name : str, optional
        The name of the molecule, by default 'test'.
    forcefield : str, optional
        The force field to use for energy calculation, by default 'MMFF94s'.
    method : str, optional
        The method to use for embedding the molecule, by default 'rdkit'.
    pre_embed : bool, optional
        Whether the conformer has been pre-embed, just read-in again, by default False.
    randomSeed : int, optional
        The random seed for reproducibility, by default 42.
    num_ring_confs : int, optional
        The number of ring conformers to generate, by default 1.
    numcores : int, optional
        The number of CPU cores to use, by default 1.
    request_alignment : str, optional
        The ring alignment in SMILES or SMARTS to support constrained docking, by default None.
    ignoreTorlib : bool, optional
        Whether to ignore the torsion library, by default False.
    threshold : float, optional
        The threshold in Angstrom for non-bonded atom distance, by default 1.6.
    mode : str, optional
        The mode for conformer sampling, by default 'fixed'.
    tolerance : float, optional
        The tolerance for dihedral angle sampling, by default 30.
    torlib : TorsionLibrary, optional
        The torsion library to use, by default Torlib.
    VERBOSE : bool, optional    
        Whether to print verbose output, by default False.

    Examples
    ----------

    >>> from msani.conformers.conformers import ConformerGenerator
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
                 name = 'test', 
                 forcefield = 'MMFF94s',
                 method = 'rdkit',
                 pre_embed = False,
                 randomSeed = 42, 
                 num_ring_confs = 1, 
                 numcores = 1, 
                 request_alignment = None,
                 ignoreTorlib = False,
                 threshold = 1.6,
                 rmsd = 0.5,
                 mode = 'fixed',
                 tolerance = 30,
                 torlib = Torlib,

                 VERBOSE = False):
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
            mode (str): Mode for conformer sampling ('fixed', 'random', 'ignoretorlib)
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
        self.rmsd = rmsd
        self.mode = mode
        self.tolerance = tolerance
        self.torlib = torlib
        self.VERBOSE = VERBOSE

        # Initialize molecule
        self._initialize_molecule()
        self._substructure_perception()
        
        # Calculate net charge
        self.netcharge = sum(atom.GetFormalCharge() for atom in self.mol.GetAtoms())
        
        # Initialize output containers
        self.ring_confs = [] # List of RDKit Mol with different ring conformations
        
        # Setup force field
        self._setup_forcefield(forcefield)
        if pre_embed == False:
            if self.method == 'corina':
                self._embed_smiles_corina()
            elif self.method == 'obabel':
                self._embed_smiles_babel()
            elif self.method == 'rdkit':
                self._embed_smiles_rdkit()
            else:
                raise ValueError(f"Invalid embedding method: {self.method}. Supported methods are: rdkit, obabel, corina.")
    
    @classmethod
    def from_existing_data(cls,
                           smiles,
                           name,
                           amsol_mol,
                           ring_confs = None,
                           mol2_str = None,
                           request_alignment = None,
                           mode:str = 'vs',
                           tolerance = 30,
                           rmsd = 0.5,
                           VERBOSE=False):
        """Alternative constructor that initializes from existing data"""

        # Create a minimal instance
        instance = cls(smiles, name=name, pre_embed=True)
        
        # Override the instance attributes
        mol = Chem.Mol(amsol_mol)
        instance.amsol_mol = mol
        instance.mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(instance.amsol_mol, mmffVariant="MMFF94s")
        instance.method = 'rdkit'
        instance.ring_confs = [Chem.Mol(ring_conf) for ring_conf in ring_confs] if ring_confs else []
        instance.mol2_str = mol2_str
        instance.sulfo_matches = utils.find_sulfonamide_like_scaffolds(instance.amsol_mol)
        instance.request_alignment = request_alignment
        instance.VERBOSE = VERBOSE
        instance.mode = mode
        instance.tolerance = tolerance
        instance.rmsd = rmsd    
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
            self.mp.SetMMFFEleTerm(False)
        # elif forcefield == 'UFF': Not supported now...
        #    self.mp = rdForceFieldHelpers.UFFGetMoleculeProperties(self.mol_H)
        
        # Store forcefield type for later use
        self.forcefield = forcefield

    def _substructure_perception(self):
        """Identify important substructures in the molecule"""
        # Get ring information
        self.planar_rings, self.non_planar_rings = utils.get_flexible_ring(
            self.mol_H, 
            SRLib
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
        self.cycloheptatriene_like = utils.find_cycloheptatriene(self.mol_H)
        
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
                print('\tFound 5-membered nitrogen aromatic rings')
                for match in self.conjugated_substituted_nitrogen_5aro: 
                    print(f'\t {match}')
            if self.conjugated_substituted_nitrogen_6aro:
                print('\tFound 6-membered nitrogen aromatic rings')
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
            
        # Determine number of initial conformations needed
        self.num_initialConfs = 50 if (self.sulfo_matches or self.non_planar_rings or self.flippable_Ns) else 10

    #================= Generation of initial conformers =========================
    # There are three methods for the generation of the initial conformers:
    # 1. RDKit
    # 2. Open Babel
    # 3. Corina    
    def _embed_smiles_rdkit(self):
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
            if not(self.flippable_Cs) and not(self.flippable_Ns): 
                # If there are flippable C or N atoms, we need to generate more to filter out the favorable
                params.pruneRmsThresh = 0.35 
            params.randomSeed = self.randomSeed # For reproducibility
            #params.useRandomCoords = True
            conf_ring_descriptors_df = DataFrame()
            if CPP_AVAILABLE:
                result_mols = cpp_sampler.embed_multiple_confs(mol = self.mol_H, 
                                                               numConfs = self.num_initialConfs, 
                                                               params = params, 
                                                               constraints = self)
                
                # Collect conformer data efficiently
                conformer_data_list = []
                for conformer in result_mols.GetConformers():
                    energy = conformer.GetDoubleProp(f'MMFF_Energy')
                    conf_data = utils.classify_confs(conformer, 
                                                          energy, 
                                                          self.non_planar_rings, 
                                                          self.flippable_Ns, 
                                                          self.flippable_Cs,
                                                          self.sulfo_matches)
                    conformer_data_list.append(conf_data)
                
                # Create DataFrame once from all collected data
                conf_ring_descriptors_df = DataFrame(conformer_data_list)
            else:
                print('C++ extension not available, please install it by create the environment again with:\n' \
                'mamba create -f environment.yml\n' \
                'conda activate msani\n' \
                'pip install -e .')
                exit(1)
            
                
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

        if len(conf_ring_descriptors_df) == 0:
            print(f"ETKDGv3 also failed for {self.name}, using OpenBabel")
            logger.warning(f"ETKDGv3 also failed for {self.name}, using OpenBabel")
            try:
                self._embed_smiles_babel()
            except RuntimeError as e:
                logger.error(f"Error in generating initial conformation using OpenBabel for {self.name}, skipping it {e}")
                log_error(self.smiles, self.name)
                return
            return
        
        conf_ring_descriptors_df.sort_values(['equatorial_subs_Ns', 'equatorial_subs_Cs', 'Energy'],
                                            ascending=[False, False, True], inplace=True) 
        

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
                    conf_idx = temp_mol.AddConformer(row.iloc[0], assignId=True)
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
            
            if  self.cycloheptatriene_like: 
                self.num_ring_confs = max(2, self.num_ring_confs) # Cycloheptatriene has two puckering ring conformations
                print(f"Found cycloheptatriene in {self.name}, setting num_ring_confs to {self.num_ring_confs}")

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

    def _embed_smiles_babel(self, timeout: int = 30):
        '''
        Embed the SMILES string using Open Babel. CLI version is used as it is found more flexible 
        than the RDKit version.

        Parameters
        ----------
        timeout : int
            Maximum number of seconds to wait for OpenBabel before killing the process (default: 30).
        '''
                                            # -h: add hs; gen3d
        cmd = [str(obabel_path), f"-:{self.smiles}", "-h", "--gen3d", "-osdf"]

        # Execute the command and capture stdout
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.communicate()  # Drain buffers and prevent zombie process
            raise RuntimeError(
                f"(OpenBabel timed out)"
            )

        # Convert the SDF output from stdout to an RDKit molecule
        mol_rdkit = Chem.MolFromMolBlock(stdout, removeHs=False)
        
        # Set molecule properties
        mol_rdkit.SetProp("_Name", self.name)
        self.mol_H = Chem.Mol(mol_rdkit) 
        conjugated_substituted_nitrogen_5aro = utils.find_conjugated_substituted_nitrogen_5aro(mol_rdkit)
        conjugated_substituted_nitrogen_6aro = utils.find_conjugated_substituted_nitrogen_6aro(mol_rdkit)

        barbiturate_matches = utils.find_barbiturates(mol_rdkit)
        hydantoin_matches = utils.find_hydantoins(mol_rdkit)
        substituted_N_barbi_hydan_like = utils.find_substituted_N_barbi_hydan_like(mol_rdkit, barbiturate_matches, hydantoin_matches)
        planar_rings, _ = utils.get_flexible_ring(mol_rdkit, SRLib) 
        self.mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_rdkit, mmffVariant="MMFF94s")
        self.mp.SetMMFFEleTerm(False)
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_rdkit, self.mp, confId=0)
        if conjugated_substituted_nitrogen_5aro:
            for a, b, c, d, e, f in conjugated_substituted_nitrogen_5aro:
                ff.MMFFAddTorsionConstraint(a, b, c, d, False, 178, 182, 1)
                ff.MMFFAddTorsionConstraint(a, b, f, e, False, 178, 182, 1)
                ff.MMFFAddTorsionConstraint(b, c, d, e, False, -2, 2, 1)
                ff.MMFFAddTorsionConstraint(d, e, f, b, False, -2, 2, 1)
        if conjugated_substituted_nitrogen_6aro:
            for a, b, c, d, e, f, g in conjugated_substituted_nitrogen_6aro:
                ff.MMFFAddTorsionConstraint(a, b, c, d, False, 178, 182, 1)
                ff.MMFFAddTorsionConstraint(a, b, g, f, False, 178, 182, 1)
                ff.MMFFAddTorsionConstraint(b, c, d, e, False, -2, 2, 1)
                ff.MMFFAddTorsionConstraint(e, f, g, b, False, -2, 2, 1)
        if barbiturate_matches:
            for match in barbiturate_matches:
                n = len(match)
                for i in range(n):
                    a, b, c, d = [match[(i + j) % n] for j in range(4)]
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, -2, 2, 1) 
        if hydantoin_matches:
            for match in hydantoin_matches:
                n = len(match)
                for i in range(n):
                    a, b, c, d  = [match[(i + j) % n] for j in range(4)]
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, -2, 2, 1)  
        if substituted_N_barbi_hydan_like:
            for a, b, c, d in substituted_N_barbi_hydan_like:
                ff.MMFFAddTorsionConstraint(a, b, c, d, False, 178, 182, 1)
        if planar_rings:
            for ring in planar_rings:
                # For planar rings, we need to ensure that the ring is planar.
                # This is done by setting the dihedral angles to +-2.
                n = len(ring)
                for i in range(n-1):
                    a, b, c, d = ring[i], ring[i + 1], ring[(i + 2) % n], ring[(i + 3) % n ]
                    ff.MMFFAddTorsionConstraint(a, b, c, d, False, -2, 2, 1)
        ff.Minimize()
        self.ring_confs = [Chem.Mol(mol_rdkit)] # Replicate the output from embed_rdkit
        self.amsol_mol = Chem.Mol(mol_rdkit) # An RDKit Mol Object with upto 10 confs for AMSOL
        mol2_obj = mol2writer.Mol2Writer(mol_rdkit)
        self.mol2_str = mol2_obj.write_mol2()

    def _embed_smiles_corina(self):
        '''
        Embed the SMILES string using CORINA and return the mol, net_charge,
        rigid_scaffolds, and flexible_scaffolds
        '''
        self.mol2_str, self.ring_confs = utils.embed_smiles_corina(self.smiles, self.name, self.num_ring_confs, self.VERBOSE)
        self.mol_H = Chem.Mol(self.ring_confs[0])
        self.mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(self.ring_confs[0], mmffVariant="MMFF94s")
        self.amsol_mol = Chem.Mol(self.ring_confs[0]) # An RDKit Mol Object with upto 10 confs for AMSOL
        self.sulfo_matches = [] # No sulfonamide flipping in CORINA

    # ========== Torsional sampling =========================
    def stochastic_sampling(self,
                            mol,
                            tolerance_level,
                            match_torlib,
                            numConfs,
                            window = 25,
                            max_attempts = 50000,
                            eps = 1,
                            random_method = 'uniform',
                            timeout_conf = 1):
        """
        
        """
        if not(CPP_AVAILABLE): 
            print('C++ extension not available, please install it by create the environment again with:\n' \
            'mamba create -f environment.yml\n' \
            'conda activate msani\n' \
            'pip install -e .')
            exit(1)
        try:
            result = cpp_sampler.stochastic_sampling_continuous(
                    mol = mol,
                    match_torlib = match_torlib,
                    tolerance_level = tolerance_level,
                    numConfs = numConfs,
                    window = window,
                    max_attempts = max_attempts,
                    rmsd = self.rmsd,
                    clash_threshold = self.threshold,
                    hetero_H_bonds = self.hetero_H_bonds,
                    randomSeed = self.randomSeed,
                    random_method = random_method,
                    mmff_variant = self.forcefield,
                    eps = eps,
                    timeout_conf = int(timeout_conf * 60),
                    verbose = self.VERBOSE

                )
            if result is not None and hasattr(result, 'GetNumConformers'):
                return(result)
        except Exception as e:
            print(f'Error in C++ extension: {e}')

    def conf_sampling(self, 
                      numConfs=2000, 
                      energywindow = 25,
                      eps = 1, 
                      ignoreTorlib = False, 
                      AllowNonRing = False, 
                      timeout_conf = 1,
                      request_alignment=None):
        """
        Perceive the allowed dihedral angles and call stochastic sampling to generate conformers.

        Notes:
            For assymetric sulfonamides, there would be two versions of rigid scaffolds handled by msani.

        Args:
            numConfs (int): Number of conformers to generate.
            energywindow (float): Energy window for conformer generation.
            eps (float): DIelectric constant for electrostatic interactions.
            ignoreTorlib (bool): Whether to ignore the torsion library.
            AllowNonRing (bool): Whether to allow the full sampling of non-ring compounds.
            request_alignment (list): List of atom indices for alignment.
        
        """
        self.conf_sampled = True
        self.mp.SetMMFFEleTerm(True) #Turn on back otherwise it would produce unfeasible conformers
        self.mp.SetMMFFDielectricConstant(eps) #Set the dielectric constant for electrostatic interactions
        if self.request_alignment is None and request_alignment is not None:
            self.request_alignment = request_alignment
        self.rot_bonds = utils.getDihedralMatches(self.ring_confs[0])

        if self.mode == 'fixed' or self.mode == 'ignoretorlib':
            # This is still experimental, call conf_samplingv2, where Torlib is read differently
            # and the angle chosen would be deterministic. By default, peak +- 30 degrees, if tol2 >= 30.
            self.conf_samplingv2(numConfs = numConfs,
                                 energywindow = energywindow,
                                 AllowNonRing = AllowNonRing,
                                 timeout_conf = timeout_conf,
                                 request_alignment = request_alignment,
                                 eps = eps,
                                 ignoreTorlib = ignoreTorlib)
            return
        
        possible_numConfs, match_torlib, self.hetero_H_bonds = utils.count_confs_by_rotbonds(
                                                    mol = self.ring_confs[0],
                                                    rot_bonds = self.rot_bonds,
                                                    amide_bonds = self.amide_linkages,
                                                    ignoretorlib = ignoreTorlib,
                                                    torlib = self.torlib,
                                                    VERBOSE = self.VERBOSE
                                                )
        

        # Find the rigid part only once outside the loop to save processing time
        self.atom_maps, self.label_map = utils.find_rigid_part(self.ring_confs[0], request_alignment)
        # Molecules which don't have rings are not of interest --> only sample limitedly.
        if request_alignment and not self.atom_maps:
            log_error(self.smiles, self.name)
            return
        # For very flexible molecules, we need to sample more, then filter by energy later
        else:
            numConfs = int(numConfs/len(self.ring_confs))
        
        if (self.label_map) and not (AllowNonRing): numConfs = 30
    

        if self.VERBOSE:
            rigid_info = f'\tFound {self.atom_maps} ({self.label_map}) as a rigid part' if self.label_map else f'\tFound {self.atom_maps} (rings) as rigid parts'
            print(rigid_info)
            print("\tMaximum possible conformers based on rotatable bonds: ", possible_numConfs)


        if (len(match_torlib) == 0):
            print('No rotatable bonds found, returning the original conformation')
            return

        for idx, mol in enumerate(self.ring_confs):
            if self.VERBOSE: print(f"\tHandling ring/sulfonamide conformation {idx+1}/{len(self.ring_confs)}")
            original_mol = Chem.Mol(mol)
            processing_mol = Chem.Mol(mol)
              
            # Only remap the match_torlib when sulfo_matches is found
            if self.sulfo_matches: 
                _, match_torlib, _ = utils.count_confs_by_rotbonds(
                                                            mol = mol,
                                                            rot_bonds = self.rot_bonds,
                                                            amide_bonds = self.amide_linkages,
                                                            ignoretorlib = ignoreTorlib,
                                                            torlib = self.torlib,
                                                            VERBOSE = self.VERBOSE
                                                            )
            if self.VERBOSE: print('\tRunning stochastic torsional sampling')
            
            result = self.stochastic_sampling( mol = processing_mol,
                                                tolerance_level = 2,
                                                match_torlib = match_torlib,
                                                numConfs = numConfs,
                                                window = energywindow,
                                                max_attempts = 50000,
                                                eps = eps,
                                                timeout_conf = timeout_conf,
                                                random_method= 'uniform'
                                                )

            if result.GetNumConformers() == 1: # No conformers are generated, use initial conformation instead
                print(f'Failed for stochastic sampling, use the original conformation')
                continue
            
            largest_ring = max(self.atom_maps, key=len)
            for confId in range(result.GetNumConformers()):
                rdMolAlign.AlignMol(result, original_mol, confId, 0, atomMap=[(i, i) for i in largest_ring])
            self.ring_confs[idx] = Chem.Mol(result)


    # Deterministic sampling
    # This is still experimental, call conf_samplingv2, where Torlib is read differently
    # and the angle chosen would be deterministic. By default, peak +- 30 degrees, if tol2 >= 30.
    def stochastic_sampling_v2(self,
                               mol,
                               angle_map,
                               score_map,
                               numConfs,
                               possible_numConfs,
                               importance_order,
                               window = 25,
                               max_attempts=50_000,
                               hetero_H_bonds = [],
                               timeout_conf = 1,
                               eps = 1):
        """
        Call the C++ extension for stochastic sampling.
        Args:
            mol (rdkit.Chem.Mol): The molecule to sample.
            angle_map (list): List of allowed angles for each rotatable bond.
            score_map (list): List of scores for each angle.
            numConfs (int): Number of conformers to generate.
            possible_numConfs (int): Maximum possible number of conformers.
            importance_order (list): Order of importance for sampling rotatable bonds.
            window (float): Energy window for conformer generation.
            max_attempts (int): Maximum number of attempts for sampling.
            hetero_H_bonds (list): List of heteroatom-hydrogen bonds to consider.
            timeout_conf (int): Timeout for conformer generation in minutes.
            eps (float): Dielectric constant for electrostatic interactions.
        """
        if not(CPP_AVAILABLE): 
            print('C++ extension not available, please install it by create the environment again with:' \
            'mamba create -f environment.yml' \
            'conda activate msani' \
            'pip install -e .')
            exit(1)
        try:

            result = cpp_sampler.stochastic_sampling_discrete(
                    mol=mol,
                    angle_map=angle_map,
                    score_map=score_map,
                    importance_order=importance_order.tolist(),
                    hetero_H_bonds=hetero_H_bonds,
                    numConfs=numConfs,
                    possible_numConfs=possible_numConfs,
                    window=window,
                    max_attempts=max_attempts,
                    timeout_conf = int(timeout_conf * 60),  # Convert minutes to seconds (int)
                    rmsd=self.rmsd,
                    randomSeed=self.randomSeed,
                    clash_threshold = self.threshold,  
                    verbose = self.VERBOSE,
                    mmff_variant = self.forcefield,
                    eps = eps,

                )
            if result is not None and hasattr(result, 'GetNumConformers'):
                return(result)
        except Exception as e:
            print(f'Error in C++ extension: {e}')

    
    def conf_samplingv2(self,
                        numConfs=2000,
                        energywindow = 25,
                        ignoreTorlib = False, 
                        AllowNonRing=False,
                        timeout_conf = 1,
                        request_alignment=None,
                        eps = 1):
        """
        Perceive the allowed dihedral angles and call stochastic sampling to generate conformers.
        """
        possible_numConfs, angle_map, score_map, self.rot_bonds, self.hetero_H_bonds = utils.count_confs_by_rotbonds_v2(mol = self.ring_confs[0],
                                                                                   rot_bonds = self.rot_bonds,
                                                                                   amide_bonds = self.amide_linkages,
                                                                                   ignoretorlib = ignoreTorlib,
                                                                                   torlib = self.torlib,
                                                                                   VERBOSE=self.VERBOSE)
        importance_order = utils.get_importance_order(self.ring_confs[0], self.rot_bonds)

        if self.VERBOSE:
            print(f"\tTheory: {possible_numConfs} possible conformations")
            print(f"\tRotatable bonds: {self.rot_bonds}")
            print(f"\tPossible angles:")
            for idx in range(len(angle_map)):
                print(f"\t{angle_map[idx]} {score_map[idx]}")

        # Find the rigid part only once outside the loop to save processing time
        self.atom_maps, self.label_map = utils.find_rigid_part(self.ring_confs[0], request_alignment)

        if request_alignment and not self.atom_maps:
            log_error(self.smiles, self.name)
            logger.error('No substructure found for the requested alignment. No conformers are generated')
            return
        else:
            # In case of sulfonamides and cycloheptatrienes, we divine the numConfs by the number of ring conformers
            numConfs = int(numConfs/len(self.ring_confs))

        # Molecules which don't have rings are not of interest --> only sample limitedly.
        if (self.label_map) and not (AllowNonRing): numConfs = 30

        if (len(self.rot_bonds) == 0):
            print('No rotatable bonds found, returning the original conformation')
            return
        
        for idx, mol in enumerate(self.ring_confs):
            if self.VERBOSE: print(f"\tHandling ring/sulfonamide conformation {idx+1}/{len(self.ring_confs)}")
            original_mol = Chem.Mol(mol)
            processing_mol = Chem.Mol(mol)
              
            # Only remap the match_torlib when sulfo_matches is found
            if self.sulfo_matches: 
                possible_numConfs, angle_map, score_map, _, _ = utils.count_confs_by_rotbonds_v2(mol = mol,
                                                                                                rot_bonds = self.rot_bonds,
                                                                                                ignoretorlib = ignoreTorlib,
                                                                                                amide_bonds = self.amide_linkages,
                                                                                                torlib = self.torlib,
                                                                                                VERBOSE = self.VERBOSE)
            if self.VERBOSE: print('\tRunning stochastic torsional sampling')
            
            result = self.stochastic_sampling_v2(mol = processing_mol,
                                                  angle_map = angle_map,
                                                  score_map = score_map,
                                                  numConfs = numConfs,
                                                  possible_numConfs = possible_numConfs,
                                                  importance_order =  importance_order, 
                                                  window = energywindow, 
                                                  max_attempts = 50_000,
                                                  hetero_H_bonds = self.hetero_H_bonds, 
                                                  timeout_conf = timeout_conf,
                                                  eps = eps)
            if result.HasProp('Failed_sampling') and result.GetProp('Failed_sampling') == '1':
                if self.VERBOSE: print(f'Failed to find any confs for {self.name}, using the random dihedral angles approach as a fallback')
                possible_numConfs, match_torlib, self.hetero_H_bonds = utils.count_confs_by_rotbonds(mol = self.ring_confs[0],
                                                             rot_bonds = self.rot_bonds,
                                                             amide_bonds = self.amide_linkages,
                                                             ignoretorlib = ignoreTorlib,
                                                             torlib = self.torlib,
                                                             VERBOSE = self.VERBOSE)
                result = self.stochastic_sampling( mol = processing_mol,
                                                tolerance_level = 2,
                                                match_torlib = match_torlib,
                                                numConfs = numConfs,
                                                window = energywindow,
                                                max_attempts = 50000,
                                                eps = eps,
                                                timeout_conf = timeout_conf,
                                                random_method= 'uniform'
                                                )
                
                if result.HasProp('Failed_sampling') and result.GetProp('Failed_sampling') == '1':
                    if self.VERBOSE: print(f'Failed for stochastic sampling for {self.name}, use the original conformation')
                    continue
            
            
            if request_alignment: 
                largest_ring = self.atom_maps[0]
            else:
                largest_ring = max(self.atom_maps, key=len)
            for confId in range(result.GetNumConformers()):
                rdMolAlign.AlignMol(result, original_mol, confId, 0, atomMap=[(i, i) for i in largest_ring])
            self.ring_confs[idx] = Chem.Mol(result)


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

        # Remove the Failed_sampling property from the SDF file if it exists
        for mol in self.ring_confs:
            if mol.HasProp('Failed_sampling'): mol.ClearProp('Failed_sampling')

        # Write the conformers to SDF file(s)
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
        Write the conformers to a Mol2 file.

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
        if not MEEKO_AVAILABLE:
            raise ImportError('Please install the meeko package using "pip install meeko" to use this script.\nIn case you are using Python >= 3.12, install it from the Github repository.')

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
        if not AMSOL_AVAILABLE:
            raise ImportError("Please install the AMSOL to msani/amsol to use this function.")
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
            log_error(self.smiles, self.name)
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
            log_error(self.smiles, self.name)
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
    """Write content to a file."""
    with open(file, "w") as f:
        f.write(content)

def write_to_tarball(ball, data, name):
    """Write data to a tarball with the specified name."""
    tar = tarfile.TarInfo(name=name)
    tar.size = len(data)
    ball.addfile(tar, io.BytesIO(data))

def gen_conf_chunk(df: DataFrame, args, input_file='0'):
    """
    Generate conformers for a given DataFrame of SMILES strings and save them in different formats.
    
    Args:
        df (DataFrame): DataFrame containing SMILES strings and other relevant information.
        args (Namespace): Parsed arguments containing various configuration options.

            Keys may include:
                randomSeed (int): Seed for random number generation.\n
                numconfs (int): Number of conformations to generate.\n
                debug (bool): Verbose output for debugging.\n
                cleanup (bool): Whether to remove intermediate files after processing.\n
                energywindow (float): Energy window for conformer sampling.\n
                timeout (int): Timeout (in minutes) for RDKit-based conformation generation.\n
                ignoretorlib (bool): Whether to ignore torsion library constraints.\n
                rigid (str): SMILES/SMARTS for conformers to be aligned to.\n
                nringconfs (int): Number of ring conformers to generate.\n
                numcores (int): Number of CPU cores to use for parallel processing.\n
                method (str): Method for initial conformation generation ('rdkit', 'obabel', 'corina').\n
                timing (bool): If enabled, logs timing information for each step.\n
                smiles (bool): If True, skips restarting logic.\n
                format (list): Output formats to generate (e.g., 'pdbqt', 'sdf', 'mol2', 'db2').\n
        input_file (str): Name of the input file (default is '0').
    
    Returns:
        None
    """
    if df.empty:
        logger.warning("Empty DataFrame provided, skipping conformation generation.")
        return
    if 'mol' not in df.columns:
        df.loc[:,'mol'] = df['smiles'].apply(Chem.MolFromSmiles)
    df = filters.Filters.remove_exotic_chem_to_db2(df)
    randomSeed, numConfs, VERBOSE, cleanup, energywindow, timeout, request_alignment, nr, numcores, mode, tolerance, allowNonring = \
        args.randomSeed, args.numconfs, args.debug, args.cleanup, args.energywindow, args.timeout, args.rigid, args.nringconfs, args.numcores, args.mode, args.tolerance, args.allowNonring
    
    ignoreTorlib = (args.mode == 'ignoretorlib')
    request_alignment = Chem.MolFromSmarts(utils.canonicalize_if_smiles(request_alignment)) if request_alignment else None
    if not(ignoreTorlib):
        if args.torsion:
            Torlib.add_custom_rules_from_file(args.torsion, debug = VERBOSE)

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
        if not(os.path.exists('msani_timing.csv')): 
            if 'sdf' in args.format:
                with open('msani_timing.csv', 'w') as f: f.write('Name,Initial embedding,Torsional sampling,SDF,Total\n')
            else:
                with open('msani_timing.csv', 'w') as f: f.write('Name,Initial embedding,AMSOL,Torsional sampling,Mol2DB2,Total\n')
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
            longname = row['longname'] if args.synthon else None
            if name in processed_mols:
                print(f"Skipping {name} as it already exists")
                logger.info(f"Skipping {name} as it already exists")
                continue
            logger.info(f"Handling {name}")
            if VERBOSE: print(f"Handling {name}")
            if args.timing: start = time.time() 
            try:
                if args.method == 'corina':
                    confgen = ConformerGenerator(smiles, name, method='corina', mode = mode, tolerance=tolerance, rmsd=args.rmsd, VERBOSE=VERBOSE)
                elif args.method == 'obabel':
                    confgen = ConformerGenerator(smiles, name, method='obabel', mode = mode, tolerance=tolerance, rmsd=args.rmsd, VERBOSE=VERBOSE)
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
                            confgen = ConformerGenerator(smiles, name, num_ring_confs=nr, method='obabel', tolerance=tolerance, rmsd=args.rmsd, VERBOSE=VERBOSE)
                        except Exception as e:
                            logger.error(f"Error in generating initial conformation using OpenBabel for {name}, skipping it {e}")
                            log_error(smiles, name)
                            continue
                        if confgen.amsol_mol is None:
                            logger.error(f"Error in generating initial conformation using OpenBabel for {name}, skipping it")
                            log_error(smiles, name)
                            continue

                    # Retrieve result from queue
                    elif not queue.empty():
                        bin_amsol_mol, bin_conf_rings, mol2_str, error = queue.get() 

                        if error:
                            logger.error(f"Error in generating initial conformation using RDKit for {name}, skipping it {error}")
                            log_error(smiles, name)
                            continue
                        confgen = ConformerGenerator.from_existing_data(smiles=smiles, 
                                                                        name=name, 
                                                                        amsol_mol=bin_amsol_mol, 
                                                                        ring_confs=bin_conf_rings, 
                                                                        mol2_str=mol2_str, 
                                                                        request_alignment=request_alignment, 
                                                                        mode=mode,
                                                                        tolerance=tolerance,
                                                                        rmsd=args.rmsd,
                                                                        VERBOSE=VERBOSE)
                    else:
                        logger.error(f"Unknown error in generating initial conformation for {name}, skipping it.")
                        log_error(smiles, name)
                        continue
            except Exception as e:
                logger.error(f"Error in generating initial conformation for {name}, skipping it: {e}")
                log_error(smiles, name)
                continue

            if args.timing: embed_time = time.time() # Time for embedding
            
            if 'pdbqt' in args.format: confgen.to_pdbqt()

            if any(format in args.format for format in ['sdf', 'mol2', 'db2', 'db2.tgz']):
                try:
                    confgen.conf_sampling(numConfs=numConfs,
                                        energywindow=energywindow,
                                        ignoreTorlib=ignoreTorlib,
                                        AllowNonRing=allowNonring,
                                        eps = args.eps,
                                        timeout_conf=args.timeout_conf,
                                        request_alignment=request_alignment,
                                        )
                except Exception as e:
                    logger.error(f"Error in conformational sampling for {name}: {e}")
                    log_error(smiles, name)
                    continue
            if args.timing: sampling_time = time.time() # Time for sampling
            if 'sdf' in args.format: 
                confgen.to_sdf()
                if args.timing: 
                    sdf_time = time.time() # Time for sdf
                    logging_time += f'{name},{embed_time-start},{sampling_time-embed_time},{sdf_time-sampling_time},{sdf_time-start}\n'

            if 'mol2' in args.format: confgen.to_mol2()
            if 'db2' in args.format:
                try: 
                    confgen.to_db2(longname = longname, env=env, cleanup=cleanup)
                except Exception as e:
                    logger.error(f"Error in converting {name} to DB2 format: {e}")
                    log_error(smiles, name)
                    continue
            if 'db2.tgz' in args.format:
                try: 
                    confgen.to_db2(longname = longname, env=env, cleanup=cleanup, tarfile = output)
                except Exception as e:
                    logger.error(f"Error in converting {name} to DB2 format: {e}")
                    log_error(smiles, name)
                    continue
                
            if args.timing and ('db2' in args.format): 
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
        with open('msani_timing.csv', 'a') as f:
            f.write(logging_time)

def main():
    parser = argparse.ArgumentParser(description="Generate conformers for a given SMILES string/file.\nTwo-column files are required.",
                                     formatter_class=CustomHelpFormatter,
                                     add_help = False)  # Suppress default -h/--help)
    parser.add_argument('--input_files', '-i', type=str, default = None, help='Input file containing SMILES strings.')
    parser.add_argument('--smiles', '-s', type = str, default = None, help='Input SMILES string.')
    parser.add_argument('--prefix', '-p', type=str, default = 'db2', help='Prefix for the output files.')
    parser.add_argument('--format', '-f', type=str, nargs='+', default=['db2.tgz'], choices=['db2', 'db2.tgz', 'pdbqt', 'sdf', 'mol2'], help='Output format(s) (e.g., db2.tgz, pdbqt, sdf, mol2).')
    parser.add_argument('--mode', '-mode', type=str, default='fixed', choices=['fixed', 'random', 'ignoretorlib'], help='Mode for conformer generation (fixed, random, ignoretorlib).')
    parser.add_argument('--tolerance', '-tol', type=float, default=30, help='Tolerance for dihedral angle sampling (default: 30).')
    parser.add_argument('--allowNonring', '-anr', action='store_true', help='Allow full sampling of non-ring compounds.')
    parser.add_argument('--eps', type=float, default=1, help='The dielectric constant for electrostatic calculations (default: 1 - vacuum).')
    parser.add_argument('--numconfs', '-nconfs', type=int, default=2000, help='Number of conformers to generate (default: 2000).')
    parser.add_argument('--nringconfs', '-nr', type=int, default=1, help='Number of ring conformers to generate (default: 1).')
    parser.add_argument('--debug', '-d', action='store_true', help='Enable verbose output for debugging.')
    parser.add_argument('--timeout', '-to',type=int, default=2, help='Timeout in minutes for RDKit-based conformation generation.')
    parser.add_argument('--timeout_conf', '-toc', type=int, default=2, help='Timeout in minutes for conformational sampling (default: 2 minutes).')
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
    parser.add_argument('--torsion', '-t', type=str, default=None, help='File containing custom torsion rules.')
    parser.add_argument('--rmsd', '-rmsd', type=float, default=0.5, help=argparse.SUPPRESS)

    args = parser.parse_args()

    
    if args.input_files and args.smiles:
        parser.error('Please provide either input files or SMILES strings, not both.')

    if args.input_files is not None:
        for inFile in args.input_files:
            if not Path(inFile).is_file():
                parser.error(f'The input file: {inFile} does not exist.')
        args.input_files = [Path(inFile).resolve() for inFile in args.input_files]
        for inFile in args.input_files:
            gen_conf_chunk(read_csv(inFile, sep=' ', header=None, names=['smiles', 'ids']), args, inFile.stem)
    else:
        if args.smiles is not None:
            smiles = args.smiles
            df = DataFrame({'smiles': [smiles], 'ids': ['0']})
            gen_conf_chunk(df, args)
        else:
            parser.error('Please provide either input files or SMILES strings.')
    # Call the main function
    #process_files(args)
   
if __name__ == "__main__":
    main()