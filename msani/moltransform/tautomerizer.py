import argparse
import multiprocessing as mp
import platform
from concurrent.futures import ProcessPoolExecutor
import logging
import shutil
from itertools import tee

from numpy import array_split, arange
from pandas import DataFrame, read_csv
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize

from msani.moltransform.neutralizer import Neutralizer
from msani.io.parsers import CustomHelpFormatter

# forkserver avoids inheriting parent locks (RDKit allocator, logging) on Linux.
# Windows only supports 'spawn'; macOS prefers 'spawn' too (fork is deprecated).
_MP_START_METHOD = 'forkserver' if platform.system() == 'Linux' else 'spawn'

RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('msani')

TAUTOMER_RULES_PATH = Path(__file__).parent.parent / 'Data' / 'tautomers_v3.txt'

TAUTOMER_PARAMS = rdMolStandardize.CleanupParameters()
TAUTOMER_PARAMS.tautomerRemoveSp3Stereo = False
TAUTOMER_PARAMS.tautomerRemoveBondStereo = False
TAUTOMER_PARAMS.tautomerRemoveIsotopicHs = False
TAUTOMER_PARAMS.maxTransforms = 1000
TAUTOMER_PARAMS.maxTautomers = 1000
TAUTOMER_PARAMS.tautomerReassignStereo = False

TE = rdMolStandardize.TautomerEnumerator(TAUTOMER_PARAMS) 

vinyl = Chem.MolFromSmarts(
    '[CX3;!$(C-[!#6&!H0])]=[CX3;!$(C-[!#6&!H0])]'
    )

integrity_substructs = {
    'amide_like': Chem.MolFromSmarts('[#6^2;$([#6](=,:[!#6])~[!#6]),$([#6](~[!#6])(~[!#6])=,:*)]'),
    'ketone_like': Chem.MolFromSmarts('[#6^2;$([#6](~[#1,#6])(=,:[#6])-[OH]),$([#6](~[#1,#6])(~[#6])=O)]'),
    'carboxylic_acid': Chem.MolFromSmarts('C(=O)[OH]'),
    'sulfoximine_like': Chem.MolFromSmarts('S=N'),
    'aliphatic_alcohol': Chem.MolFromSmarts('[OH,SH&X2]-[C;$([CH^3](-[C;!$(C=[!C])])-[#6;!$([C&z1]=[!C])]),$([CH^3]([OH])(-[O,N,S])-[#6]),$([CH2]),$([C!H0^3]-[CX3;!$(C-[!#6&!H0])]=[CX3;!$(C-[!#6&!H0])])]'),
    'phenol': Chem.MolFromSmarts('[OH]-c1ccccc1'),
    'ethylene_like': Chem.MolFromSmarts('[C;$([CX4H2;!$(C-[!#6])!$(C-C=[!#6])!$(C-C=C-C=[!#6])]-[CX4;!$(C-[!#6])]),$([CX4H2;!$(C-[!#6])]-[CX4;!$(C-[!#6])!$(C-C=[!#6])!$(C-C=C-C=[!#6])])]'),
    'methyl' : Chem.MolFromSmarts('[CH3]')

}

try:
    substructure_terms = rdMolStandardize.GetDefaultTautomerScoreSubstructs()
    del substructure_terms[8] #Methyl rule. We don't want to penalize terminal alkenes.
    del substructure_terms[2] #C=,:O rule. We introduce our own rules for carbonyls.
    del substructure_terms[0] #benzoquinone rule. We introduce our own rules for benzoquinones.
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("aro-5", "a1aaaa1", 100)
        )
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("aro-6", "a1aaaaa1", 100)
        )   
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("aro-7", "a1aaaaaa1", 100)
        )
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("thiophene", "[#6]:[#16]", 2)
        )
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("carbonyl", "[#6]=[#8,#16]", 2)
        )
    # Fix aromatic benzoquinone
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("aromatic benzoquinone", "[c;$(c1(=A)c(:a)cccc1),$(c1(=A)ccc(:a)cc1)]1[c;!$(c12c(=A)c3c(cccc3)cc1ccc(:a)c2)!$(c12c(=A)c3c(cccc3)cc1cc(:a)cc2)!$(c12c(=A)c3c(cc(:a)cc3)cc1cccc2)!$(c12c(=A)c3c(ccc(:a)c3)cc1cccc2)!$(c12cc3c(cc(:a)cc3)cc1cccc2=A)!$(c12cc3c(ccc(:a)c3)cc1cccc2=A)!$(c1c(=A)c2c(cc3c(ccc(:a)c3)c2)cc1)!$(c1c(=A)c2c(cc3c(cc(:a)cc3)c2)cc1)]cccc1", -100)
        )
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("p-benzoquinone", "[#6]1(-[#6]=,:[#6]-[#6](-[#6]=,:[#6]-1)=,:[N,S,O])=,:[N,S,O]", 94)
        )
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("o-benzoquinone", "[#6]1(-[#6](=,:[N,S,O])-[#6]=,:[#6](-[#6]=,:[#6]-1))=,:[N,S,O]", 94)
        )
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("2-or-3-OH furane and di-OH-pyrrols", "[a;$(c1([OH])[c!$(c~[OX1,OH,SX1,SX2H])][o,s][c!$(c~[OX1,OH,SX1,SX2H])][c!$(c~[OX1,OH,SX1,SX2H])]1),$([c!$(c~[OX1,OH,SX1,SX2H])]1c([OH])[o,s][c!$(c~[OX1,OH,SX1,SX2H])][c!$(c~[OX1,OH,SX1,SX2H])]1),$(c1([OH])[nX3]c([OH])[c!$(c~[OX1,OH,SX1,SX2H])][c!$(c~[OX1,OH,SX1,SX2H])]1),$(c1([OH])[nX3][c!$(c~[OX1,OH,SX1,SX2H])]c([OH])[c!$(c~[OX1,OH,SX1,SX2H])]1)]1aaaa1", -99)
        )
except AttributeError as e:
    from rdkit import rdBase
    rdkit_version = rdBase.rdkitVersion
    print(f'The new tautomerizer requires RDKit version >= 2024.9.3, your current version: {rdkit_version}')
    print(f'Use `pip install rdkit>=2024.9.3`')
    exit(1)


def score_func(mol):
    """Customized scoring function for tautomerizer: from
    https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/MolStandardize/Wrap/testMolStandardize.py"""

    return (rdMolStandardize.ScoreHeteroHs(mol) + rdMolStandardize.ScoreSubstructs(mol, substructure_terms))
            

def pairwise(iterable):
    """Utility function to iterate in a pairwise fashion."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def check_configurations(matches, mol):
    ''' A helper function to check the configurations of the vinyl bonds in the molecule.'''
    config = []
    for match in matches:
        for bond_idx in (pairwise(match)):
            bond = mol.GetBondBetweenAtoms(bond_idx[0], bond_idx[1])
            config.append(bond.GetBondType())
    return config

logo =r""" _____           _                            _              
|_   _|         | |                          (_)             
  | | __ _ _   _| |_ ___  _ __ ___   ___ _ __ _ _______ _ __ 
  | |/ _` | | | | __/ _ \| '_ ` _ \ / _ \ '__| |_  / _ \ '__|
  | | (_| | |_| | || (_) | | | | | |  __/ |  | |/ /  __/ |   
  \_/\__,_|\__,_|\__\___/|_| |_| |_|\___|_|  |_/___\___|_|   
                                        From MolSanitizer
                                        """

class Tautomerizer:
    """
    A class for tautomerization of molecules.
    
    Parameters
    ----------

    neutralize: bool, default True.
        Whether to neutralize the input molecule before tautomer standardization and enumeration.

    taurdkit: bool, default True.
        Whether to use RDKit tautomer enumerator scoring function to canonicalize the input tautomer.

    numcores: int, default 1.
        Number of cores to use for multiprocessing. If 1, no multiprocessing is used.
    
    debug: bool, default False.
        Print debug message

    Examples
    ----------
    
    >>> from msani.moltransform.tautomerizer import Tautomerizer
    >>> tautomerizer = Tautomerizer(numcores= 4, neutralize= False)
    
    Tautomerize a molecule from a SMILES

    >>> tautomers = tautomerizer.tautomerize(smiles='c1ccccc1O')

    Tautomerize a molecule from an RDKit Mol object

    >>> mol = Chem.MolFromSmiles('c1ccccc1O')
    >>> tautomers = tautomerizer.tautomerize(mol = mol)

    Tautomerize a DataFrame of molecules with SMILES strings:

    >>> tautomers_df = tautomerize_df(df, smiles_column = 'smiles', name_column = 'ids')"""

    def __init__(self,
                 smartsFile = TAUTOMER_RULES_PATH,
                 taurdkit=True,
                 neutralize=True,
                 numcores = 1,
                 debug = False):

        self.debug = debug
        self.taurdkit = taurdkit
        self.neutralize = neutralize
        self.numcores = numcores
        if smartsFile is None: smartsFile = TAUTOMER_RULES_PATH
        self.reactions = self.load_reactions(smartsFile)
        self.integrity_substructs = integrity_substructs
        self.standardizing_reactions = [r for r in self.reactions if not r[1]]
        self.enumerating_reactions = [r for r in self.reactions if r[1]]
        self.debug = debug

    def __repr__(self):
        cls_name = self.__class__.__name__
        attrs = ', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())
        return f'{cls_name}({attrs})'
    
    def load_reactions(self, file_path: str):
        """Load the reactions from a file containing SMARTS strings.

        Args:
            file_path (str): Path to the file containing the reactions in SMARTS strings.

        Returns:
            list: A list containing the reactions in the form of [rdkit.Chem.rdChemReactions object, name].
        """
        reactions = []
        with open(file_path, 'r') as file:
            for line in file:
                # Skip the comments and silent rules
                if line.startswith('#'): 
                    continue
                smarts = line.strip().split()
                if smarts:
                    try:
                        # Name, Enumerate, Reaction SMARTS
                        reactions.append((smarts[0], bool(int(smarts[1])), AllChem.ReactionFromSmarts(smarts[2])))
                    except: 
                        logger.error(f"Error loading reaction: {line}")
        if self.debug:
            logger.info(f"Loaded {len(reactions)} reactions from {file_path}")
        return reactions
    
    def count_defined_stereo_doublebonds(self, mol):
        num_stereo = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                stereo = bond.GetStereo()
                if stereo != rdchem.BondStereo.STEREONONE:
                    num_stereo += 1
        return num_stereo
    
    @staticmethod
    def _get_substruct_configurations(mol, pattern, ref_match=None):
        """Returns substructure match and its configuration if present."""
        matches = mol.GetSubstructMatches(pattern)
        if ref_match is not None:
            return (matches, check_configurations(ref_match, mol)) if matches else ([], None)
        return (matches, check_configurations(matches, mol)) if matches else ([], None)
    
    def tautomer_canonicalize_rdkit(self, mol: Chem.Mol):
        """Tautomerize the input molecule using the RDKit TautomerEnumerator class.
        This is the first step to make the input result in a canonical tautomer for later fixes."""

        try:
            #canonical_tautomer = te.Canonicalize(mol, score_func)
            canonical_tautomer = None
            tautomers = [] # The tautomers would be list of (mol, score)
            max_score = -9999
            original_mol_substructs = dict()
            for name, substruct in self.integrity_substructs.items():
                original_mol_substructs[name] = set(mol.GetSubstructMatches(substruct))
            
            initial_chiral_centers = len(Chem.FindMolChiralCenters(mol))
            intial_defined_double_bonds = self.count_defined_stereo_doublebonds(mol) 

            if self.debug: 
                for name, substruct in original_mol_substructs.items():
                    print(f"\tInitial {name} substructures: {len(substruct)}")
                print(f"\tInitial defined double bonds: {intial_defined_double_bonds}")
                print(f"\tInitial chiral centers: {initial_chiral_centers}")

            # Enumerate the tautomers
            for tau in TE.Enumerate(mol):
                # Check chiral centers
                tau_chiral_centers = len(Chem.FindMolChiralCenters(tau))
                if tau_chiral_centers < initial_chiral_centers:
                    if self.debug: print(f"\t{Chem.MolToSmiles(tau)} has fewer chiral centers than the input, skipping")
                    continue
                
                # Check stereo double bonds
                tau_defined_double_bonds = self.count_defined_stereo_doublebonds(tau)
                if tau_defined_double_bonds < intial_defined_double_bonds:
                    if self.debug: print(f"\t{Chem.MolToSmiles(tau)} has fewer defined double bonds than the input, skipping")
                    continue
                
                broken = False
                for name, substruct in self.integrity_substructs.items():
                    if len(original_mol_substructs[name]) > 0:
                        tau_substruct = set(tau.GetSubstructMatches(substruct))
                        if not(original_mol_substructs[name].issubset(tau_substruct)):
                            if self.debug: print(f"\t{Chem.MolToSmiles(tau)} has lost {len(original_mol_substructs[name]-tau_substruct)} {name} substructures, skipping")
                            broken = True
                            break
                        
                if broken: continue #Any of the integrity test failed, skip this tautomer

                # All passed, add to the list
                score = score_func(tau)
                tautomers.append((tau, score))

                if score > max_score: max_score = score

            # Sort the tautomers by score
            tautomers.sort(key=lambda x: (-x[1], Chem.MolToSmiles(x[0])))

            if self.debug:
                print(f"\tInitial score: {score_func(mol)}") # To avoid the warning of not having a score function
                print(f"\tFound {len(tautomers)} tautomers")
                for t in tautomers:
                    print(f"\t\t{Chem.MolToSmiles(t[0])} {t[1]}")
                print(f"\tMax score: {max_score}")
            
            # If the canonical tautomer is the same SCORE as the input,
            # we believe more in the input than the output.
            # Return the input molecule
            initial_score = score_func(mol)
            if max_score == initial_score:
                if self.debug: print(f"Same score, use input molecule")
                return mol
            
            # Emulate the Canonicalize function
            # Pick the one that has the same "configuration" of the double bonds as the input molecule
            # Lexicographically min first
            equal_tautomers = [(t[0], t[1], Chem.MolToSmiles(t[0])) for t in tautomers if t[1] >= max_score - 4] 
            
            
            if len(equal_tautomers) > 1:
                if self.debug: 
                    print(f"\tFound {len(equal_tautomers)} tautomers with the nearly similar score:")
                    for t in equal_tautomers:
                        print(f'\t  {t[2]} {t[1]}')
                ref_vinyl, ref_conf_vinyl = Tautomerizer._get_substruct_configurations(mol, vinyl)
                if ref_vinyl:
                    if self.debug: 
                        print(f"\tReference vinyl configurations: {ref_conf_vinyl}")
                    potential_pool = []
                    for tautomer, score, smiles in (equal_tautomers): 
                        match_vinyl, config_vinyl = Tautomerizer._get_substruct_configurations(tautomer, vinyl, ref_vinyl)
                        # print('\t'+smiles)
                        # print(f"\tAcrylic: {match_acrylic} {ref_acrylic}")
                        # print(f"\tVinyl: {match_vinyl} {ref_vinyl}")
                        if len(match_vinyl) != len(ref_vinyl): continue
                        # Check for the same bonds as in the reference, whether the bond configurations changed
                        if config_vinyl == ref_conf_vinyl:
                            # if self.debug: print(f'\tChanged to {smiles}')
                            if self.debug: print(f'\tAdding {smiles} to the potential pool')
                            potential_pool.append((tautomer, score ,smiles))
                            continue
                            # canonical_tautomer = tautomer
                            # break
                    if potential_pool:
                        if any(Chem.MolToSmiles(mol) == smiles for _, _, smiles in potential_pool):
                            if initial_score == potential_pool[0][1]:
                                if self.debug: 
                                    print(f"\tInput has the highest score, take it.")
                                canonical_tautomer = mol
                            else:
                                if self.debug:
                                    print(f"\tInput is in the potential pool, but less scored, take first one.")
                                canonical_tautomer = potential_pool[0][0]
                        else:
                            if self.debug: 
                                print(f'\tInput not found, use the first one in the pool.')
                            canonical_tautomer = potential_pool[0][0]
                    if canonical_tautomer is None:
                        if self.debug: print(f"\tNone of the tautomers have the same configuration as the input molecule.")
                        canonical_tautomer = equal_tautomers[0][0]    
                else: 
                    # No vinyl bonds found, prioritize the tautomers also without vinyl bonds
                    canonical_tautomer = None
                    if self.debug: print(f"\tNo vinyl bonds found, picking the first one that also has no vinyl bond.")
                    for tautomer, _, smiles in (equal_tautomers):
                        match_vinyl, _ = Tautomerizer._get_substruct_configurations(tautomer, vinyl)
                        if not match_vinyl:
                            canonical_tautomer = tautomer
                            if self.debug: print(f'\tChanged to {smiles}')
                            break
                    # A fallback if no tautomer without vinyl bonds is found
                    if canonical_tautomer == None: 
                        if self.debug: print(f"\tNo tautomer without vinyl bonds found, picking the first one.")
                        canonical_tautomer = equal_tautomers[0][0]
            else:
                # No equal tautomers found, just pick the first one
                canonical_tautomer = equal_tautomers[0][0]

        except Exception as e:
            if self.debug: print(f"Error canonicalizing molecule: {Chem.MolToSmiles(mol)} {e}")
            logger.info(f"Error canonicalizing molecule: {Chem.MolToSmiles(mol)}")
            return mol
        return canonical_tautomer

    def standardize(self, mol:Chem.Mol):
        """
        Standardize the input molecule using the standardizing reactions.

        Args:
            mol (rdkit.Chem.rdchem.Mol): The reactant molecule.

        Returns:
            None. The unique product SMILES strings are added to self.collection.
        """
        for name, _, rxn in self.standardizing_reactions:
            rxn.Initialize()
            runs = 0
            # 10 is likely enough to avoid infinite loops. Minimum rule needs to match 3 atoms. 
            # If the molecule needs to apply more than 10 times the same rule, it is unlikely a leadlike molecule.
            while mol.HasSubstructMatch(rxn.GetReactantTemplate(0)) and runs < 10:
                if self.debug: 
                    logger.info(f"\tApplying {name} to {Chem.MolToSmiles(mol)}")
                    print(f"\tApplying {name} to {Chem.MolToSmiles(mol)}")
                new_mol = rxn.RunReactants((mol,))[0][0]
                error = Chem.SanitizeMol(new_mol, catchErrors=True)
                if error == 0:
                    mol = new_mol
                else:
                    if self.debug: 
                        print(f"\tError sanitizing molecule: {Chem.MolToSmiles(mol)}")
                    break
                runs += 1

        return mol

    def enumerate(self, mol: Chem.Mol, max_tautomers: int = 10):
        """
        Enumerate different combination of different tautomer substructures of a molecule.

        Args:
            mol (rdkit.Chem.rdchem.Mol): The reactant molecule.
            max_tautomers (int): Hard limit on maximum tautomers to prevent combinatorial explosion.
        Returns:
            list: A list of unique SMILES strings of the tautomer substructures.
        """
        unique_smiles = [Chem.MolToSmiles(mol)]
        for name, _, rxn in self.enumerating_reactions:
            #rxn.Initialize()
            i = 0
            while i < len(unique_smiles):
                if len(unique_smiles) >= max_tautomers:
                    if self.debug: print(f"\tReached max_tautomers ({max_tautomers}) limit, stopping enumeration.")
                    break
                    
                current_smiles = unique_smiles[i]
                current_mol = Chem.MolFromSmiles(current_smiles)
                products = rxn.RunReactants((current_mol,))
                if products:
                    if self.debug: print(f"\tApplying {name} to {current_smiles}")
                    for product in products:
                        new_mol = product[0]
                        error = Chem.SanitizeMol(new_mol, catchErrors=True)
                        if error == 0:
                            new_smiles = Chem.MolToSmiles(new_mol)
                            if new_smiles not in unique_smiles:
                                unique_smiles.append(new_smiles)
                                if len(unique_smiles) >= max_tautomers:
                                    break # Instantly break inner loop if cap hit
                i += 1
            if len(unique_smiles) >= max_tautomers:
                break # Break outer reaction loop as well if cap hit
                
        return unique_smiles


    def tautomerize(self, 
                    smiles:str = None, 
                    mol: Chem.Mol = None,
                    name: str = None) -> list:
        """
        Tautomerize the input molecule. (Either SMILES or RDKit molecule object must be provided)

        Args: 
            
            smiles (str): SMILES string of the molecule.
            mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.
            name (str): Name of the molecule. - mainly for debugging purposes.

        Returns:
            Returns a list SMILES strings of the tautomers.
        """
        
        if (mol is None) and (smiles is None):
            raise ValueError("Either SMILES or RDKit molecule object must be provided.")
        
        if mol is None and smiles:
            mol = Chem.MolFromSmiles(smiles)
        else:
            smiles = Chem.MolToSmiles(mol)
        if self.debug:
            print(f"Tautomerizing {name}...")
            logger.info(f"Tautomerizing {name}...")
        if self.neutralize:
            mol = Neutralizer.neutralize_mol(mol)
        if self.taurdkit:
            mol = self.tautomer_canonicalize_rdkit(mol)
            
        # Step 2: Standardize the molecule
        standardized_mol = self.standardize(mol)
        if self.debug:
            print(f"\tStandardized to {Chem.MolToSmiles(standardized_mol)}")
        unique_tautomers = self.enumerate(standardized_mol)
        # for _, tautomer in enumerate(unique_tautomers):
        #     unique_tautomers[_] = tautomer.replace('[C]', 'C').replace('[H]', 'H').replace('[CH]', 'C')\
        #         .replace('[N]', 'N').replace('[O]', 'O').replace('[S]', 'S').replace('[cH]', 'c').replace('[n]', 'n')
        return unique_tautomers
    
    def tautomerize_df_mp(self, 
                          df: DataFrame,
                          smiles_column: str = 'smiles',
                          mol_column: str = 'mol',
                          name_column: str = 'ids') -> DataFrame:
        """
        Tautomerize a dataframe of molecules using ProcessPoolExecutor with
        forkserver start method and an initializer to avoid:
          1. Pipe-buffer deadlocks from imap_unordered with large results.
          2. RDKit/logging lock inheritance via os.fork() on Linux.
          3. Re-pickling the Tautomerizer object for every chunk.
        """
        if mol_column not in df.columns:
            df = df.copy()  # Make a single shallow copy to safely add mol_column
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        num_cores = min(self.numcores, len(df))  # Avoid using more cores than data chunks

        # Use numpy array_split to produce contiguous index slices (no data copies)
        chunk_indices = array_split(arange(len(df)), num_cores)
        chunks = [df.iloc[indices] for indices in chunk_indices]

        # _MP_START_METHOD: forkserver on Linux (no inherited locks), spawn elsewhere.
        ctx = mp.get_context(_MP_START_METHOD)
        with ProcessPoolExecutor(
                max_workers=num_cores,
                mp_context=ctx,
                initializer=_init_tautomerizer_worker,
                initargs=(self,)) as executor:
            # list() drains eagerly via an internal thread — safe from pipe-block.
            batch_results = list(executor.map(_tautomerize_chunk_worker, chunks))

        results = [row for batch in batch_results if batch for row in batch]
        return DataFrame(results)

    def tautomerize_df(self, 
                       df: DataFrame,
                       smiles_column: str = 'smiles',
                       mol_column: str = 'mol',
                       name_column: str = 'ids') -> DataFrame:
        """
        Tautomerize a dataframe of molecules using multiprocessing or single-core based on `num_cores`.

        Args:
            df (DataFrame): The input DataFrame containing SMILES strings.
            smiles_column (str): The column name containing SMILES strings (default: 'smiles').
            mol_column (str): The column name for RDKit molecule objects (default: 'mol').
            name_column (str): The column name for molecule identifiers (default: 'ids').
        Returns:
            DataFrame: A DataFrame with tautomerized molecules, including SMILES and identifiers.
        """
        if df.empty:
            return df
        if mol_column not in df.columns:
            df = df.copy() # Make a shallow copy to safely add mol_column
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)
            
        if self.numcores > 1:
            return self.tautomerize_df_mp(df, smiles_column, mol_column, name_column)
        else:
            return DataFrame(_process_tautomer_rows(df, self, smiles_column, mol_column, name_column))

                       

# ---------------------------------------------------------------------------
# Module-level worker helpers for ProcessPoolExecutor
# These must be at module level (not nested) to be picklable.
# ---------------------------------------------------------------------------

_tautomerizer_worker = None  # per-worker singleton set by the initializer


def _init_tautomerizer_worker(tautomerizer):
    """Initializer run once per worker process."""
    global _tautomerizer_worker
    _tautomerizer_worker = tautomerizer


def _tautomerize_chunk_worker(chunk):
    """Top-level worker function dispatched by ProcessPoolExecutor."""
    return _process_tautomer_rows(chunk, _tautomerizer_worker,
                                  'smiles', 'mol', 'ids')


def _process_tautomer_rows(df, tautomerizer, smiles_column, mol_column, name_column):
    """
    Common function to process tautomerization for both single-core and multiprocessing.
    """
    results = []
    for _, row in df.iterrows():
        mol = row[mol_column]
        longname = row.get('longname', None)
        original_idx = row.get('original_idx', None)  # For debugging purposes
        tautomers_smiles = tautomerizer.tautomerize(mol=mol, name = row[name_column])

        for _, tautomer in enumerate(tautomers_smiles):
            results.append({name_column: row[name_column],
                            mol_column: Chem.MolFromSmiles(tautomer),
                            smiles_column: tautomer,
                            'longname': longname,
                            'original_idx': original_idx
                            })

    return results


def main():
    print(logo)
    parser = argparse.ArgumentParser(description='Tautomerize a molecule', formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input', default=None, type=str, help='Input file containing SMILES strings and names')
    parser.add_argument('-s', '--smiles', default=None, type=str, help='SMILES string of the molecule')
    parser.add_argument('-l', '--library', default=TAUTOMER_RULES_PATH, type=str, help='Tautomer rules file (default: msani/Data/tautomers_v3.txt)')
    parser.add_argument('-t', '--template', action='store_true', help='Create a template for tautomer rules')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('-o', '--output', default='tautomers_output.smi', type=str, help='Output file')
    parser.add_argument('-j', '--numcores', default=4, type=int, help='Number of cores to use for multiprocessing')
    args = parser.parse_args()
    if args.smiles and args.input:
        raise ValueError("Either SMILES or input file must be provided.")
    if args.template:
        # Copy the tautomer rules file to the current directory
        output_filename = f"tautomer_template_msani.txt"
        shutil.copy(TAUTOMER_RULES_PATH, output_filename)
        print(f"Successfully create a template for tautomerization to {output_filename}")
        return
    tautomerizer = Tautomerizer(smartsFile=args.library,
                                numcores=args.numcores,
                                debug=args.debug)
    if args.smiles:
        tautomers = tautomerizer.tautomerize(smiles=args.smiles)
        print(f"Tautomerized SMILES:")
        for tautomer in tautomers: print(tautomer)
    elif args.input:
        df = read_csv(args.input, names = ['smiles', 'ids'], sep = r'\s+', header=None)
        tautomers = tautomerizer.tautomerize_df(df)
        tautomers[['smiles', 'ids']].to_csv(args.output, header=False, sep = ' ', index=False)



if __name__=="__main__":
    main()