import argparse
import multiprocessing as mp
import logging
import shutil
from functools import partial
from itertools import tee

from pandas import DataFrame, read_csv
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize

from msani.moltransform.neutralizer import Neutralizer
from msani.io.parsers import CustomHelpFormatter

RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('msani')

TAUTOMER_RULES_PATH = Path(__file__).parent.parent / 'Data' / 'tautomers_v3.txt'

TAUTOMER_PARAMS = rdMolStandardize.CleanupParameters()
TAUTOMER_PARAMS.tautomerRemoveSp3Stereo = False
TAUTOMER_PARAMS.tautomerRemoveBondStereo = False
TAUTOMER_PARAMS.tautomerRemoveIsotopicHs = False
TAUTOMER_PARAMS.maxTransforms = 1000
TAUTOMER_PARAMS.maxTautomers = 1000

TE = rdMolStandardize.TautomerEnumerator(TAUTOMER_PARAMS) 

allylic = Chem.MolFromSmarts(
    '[CX4&!H0;!$(C-[!#6&!H0])]-[CX3;!$(C-[!#6&!H0])]=!@[CX3;!$(C-[!#6&!H0])]'
    )

try:
    substructure_terms = rdMolStandardize.GetDefaultTautomerScoreSubstructs()
    del substructure_terms[8] #Methyl rule. We don't want to penalize terminal alkenes.
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("amide", "[NH1,NH2]-C=O", 1)
        )
    substructure_terms.append(
        rdMolStandardize.SubstructTerm("corr_rdkit_feature1", "a1:a:a2:a:a:a:a:a-2:a:1", 199)
        )
    #substructure_terms.append(rdMolStandardize.SubstructTerm("aromatic methylidene", "c=C", -1))
except AttributeError as e:
    from rdkit import rdBase
    rdkit_version = rdBase.rdkitVersion
    print(f'The new tautomerizer requires RDKit version >= 2024.9.3, your current version: {rdkit_version}')
    print(f'Use `pip install rdkit>=2024.9.3`')
    exit(1)


def score_func(mol):
    """Customized scoring function for tautomerizer: from
    https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/MolStandardize/Wrap/testMolStandardize.py"""

    return (rdMolStandardize.ScoreRings(mol) + rdMolStandardize.ScoreHeteroHs(mol) +
            rdMolStandardize.ScoreSubstructs(mol, substructure_terms))

def pairwise(iterable):
    """Utility function to iterate in a pairwise fashion."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def check_configurations(matches, mol):
    ''' A helper function to check the configurations of the allylic bonds in the molecule.'''
    config = []
    for match in matches:
        for bond_idx in (pairwise(match)):
            bond = mol.GetBondBetweenAtoms(bond_idx[0], bond_idx[1])
            config.append(bond.GetBondType())
    return config

logo = """ _____           _                            _              
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

    Example use:
    ----------
    
    >>> from msani.moltransform.tautomerizer import Tautomerizer\n
    >>> tautomerizer = Tautomerizer(numcores= 4, neutralize= False)\n
    Tautomerize a molecule from a SMILES
    >>> tautomers = tautomerizer.tautomerize(smiles='c1ccccc1O')\n

    Tautomerize a molecule from an RDKit Mol object
    >>> mol = Chem.MolFromSmiles('c1ccccc1O')\n
    >>> tautomers = tautomerizer.tautomerize(mol = mol)\n

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
        
    def tautomer_canonicalize_rdkit(self, mol: Chem.Mol):
        """Tautomerize the input molecule using the RDKit TautomerEnumerator class.
        This is the first step to make the input result in a canonical tautomer for later fixes."""

        try:
            #canonical_tautomer = te.Canonicalize(mol, score_func)
            canonical_tautomer = None
            tautomers = [] # The tautomers would be list of (mol, score)
            max_score = -9999
            
            # Enumerate the tautomers
            for tau in TE.Enumerate(mol):
                score = score_func(tau)
                tautomers.append((tau, score))
                if score > max_score: max_score = score

            if self.debug:
                print(f"\tInitial score: {score_func(mol)}") # To avoid the warning of not having a score function
                print(f"\tFound {len(tautomers)} tautomers")
                for t in tautomers:
                    print(f"\t\t{Chem.MolToSmiles(t[0])} {t[1]}")
                print(f"\tMax score: {max_score}")
            
            # If the canonical tautomer is the same SCORE as the input,
            # we believe more in the input than the output.
            # Return the input molecule
            if max_score == score_func(mol):
                if self.debug: print(f"Same score, use input molecule")
                return mol
            
            # Emulate the Canonicalize function
            # Pick the one that has the same "configuration" of the double bonds as the input molecule
            # Lexicographically min first
            equal_tautomers = sorted(
                [(t[0], t[1], Chem.MolToSmiles(t[0])) for t in tautomers if t[1] == max_score],
                key=lambda x: x[2]
                ) 
            
            if len(equal_tautomers) > 1:
                if self.debug: print(f"\tFound {len(equal_tautomers)} tautomers with the same score: {[t[2] for t in equal_tautomers]}")
                matches = mol.GetSubstructMatches(allylic)
                if matches:
                    reference_configuration = check_configurations(matches, mol)
                    if self.debug: print(f"\tReference configuration: {reference_configuration}")
                    for tautomer in (equal_tautomers): 
                        # Pick the lexicographically min one with the same configuration as reference
                        tautomer_configuration = check_configurations(matches, tautomer[0])
                        if tautomer_configuration == reference_configuration:
                            if self.debug: print(f'\tChanged to {tautomer[2]}')
                            canonical_tautomer = tautomer[0]
                            break
                    if canonical_tautomer is None:
                        if self.debug: print(f"\tNone of the tautomers have the same configuration as the input molecule.")
                        canonical_tautomer = equal_tautomers[0][0]    
                else: 
                    # No allylic bonds found, prioritize the tautomers also without allylic bonds
                    canonical_tautomer = None
                    if self.debug: print(f"\tNo allylic bonds found, picking the first one that also has no allylic bond.")
                    for tautomer in (equal_tautomers):
                        matches = tautomer[0].GetSubstructMatches(allylic)
                        if not matches:
                            canonical_tautomer = tautomer[0]
                            if self.debug: print(f'\tChanged to {tautomer[2]}')
                            break
                    # A fallback if no tautomer without allylic bonds is found
                    if canonical_tautomer == None: 
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
                    if self.debug: print(f"\tError sanitizing molecule: {Chem.MolToSmiles(mol)}")
                    logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)}")
                    break
                runs += 1

        return mol

    def enumerate(self, mol: Chem.Mol):
        """
        Enumerate different combination of different tautomer substructures of a molecule.
        Args:
            mol (rdkit.Chem.rdchem.Mol): The reactant molecule.
        Returns:
            list: A list of unique SMILES strings of the tautomer substructures.
        """
        unique_smiles = [Chem.MolToSmiles(mol)]
        for name, _, rxn in self.enumerating_reactions:
            #rxn.Initialize()
            i = 0
            while i < len(unique_smiles):
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
                i += 1
        return unique_smiles


    def tautomerize(self, 
                    smiles:str = None, 
                    mol: Chem.Mol = None,
                    name: str = None) -> list:
        """
        Tautomerize the input molecule.
        Args: (Either SMILES or RDKit molecule object must be provided)
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
            initial_mol = Chem.Mol(mol)
            # Step 1: Use RDKit TautomerEnumerator to canonicalize the input molecule
            rdkit_canonical = self.tautomer_canonicalize_rdkit(mol)
            initial_chiral_centers = len(Chem.FindMolChiralCenters(initial_mol))
            rdkit_chiral_centers = len(Chem.FindMolChiralCenters(rdkit_canonical))
            # Canonicalization loses the stereochemistry,
            # so we keep the original molecule if it has more chiral centers
            if self.debug: 
                print(f"\tInitial chiral centers: {initial_chiral_centers}, RDKit canonical chiral centers: {rdkit_chiral_centers}")
            
            # Number of defined double bonds
            intial_defined_double_bonds = self.count_defined_stereo_doublebonds(initial_mol) 
            # Number of defined double bonds
            rdkit_defined_double_bonds = self.count_defined_stereo_doublebonds(rdkit_canonical) 
            if self.debug: 
                print(f"\tInitial defined double bonds: {intial_defined_double_bonds}, RDKit canonical defined double bonds: {rdkit_defined_double_bonds}")
            # If the canonical tautomer has less defined double bonds, we keep the original molecule
            if rdkit_defined_double_bonds < intial_defined_double_bonds or\
                rdkit_chiral_centers < initial_chiral_centers:
                if self.debug: print("\tViolation of stereochemistry, using the original molecule")
                mol = initial_mol 
            else: mol = rdkit_canonical
        
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
        Tautomerize a dataframe of molecules using multiprocessing with chunked processing.
        """
        if mol_column not in df.columns:
            df.loc[:,mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        num_cores = min(self.numcores, len(df))  # Avoid using more cores than data chunks
        chunks = [df.iloc[i::num_cores] for i in range(num_cores)]

        process_func = partial(_process_tautomer_rows, tautomerizer=self,
                               smiles_column=smiles_column,
                               mol_column=mol_column,
                               name_column=name_column)

        results = []
        with mp.Pool(processes=num_cores) as pool:
            async_results = [pool.apply_async(process_func, (chunk,)) for chunk in chunks]

            for async_result in async_results:
                try:
                    results.extend(async_result.get())  # Timeout for safety
                except Exception as e:
                    logger.error(f"Error processing a tautomer batch: {str(e)}")

        return DataFrame(results)

    def tautomerize_df(self, 
                       df: DataFrame,
                       smiles_column: str = 'smiles',
                       mol_column: str = 'mol',
                       name_column: str = 'ids') -> DataFrame:
        """
        Tautomerize a dataframe of molecules using multiprocessing or single-core based on `num_cores`.
        """
        if df.empty:
            return df
        if mol_column not in df.columns:
            df.loc[:,mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)
        if self.numcores > 1:
            return self.tautomerize_df_mp(df, smiles_column, mol_column, name_column)
        else:
            return DataFrame(_process_tautomer_rows(df, self, smiles_column, mol_column, name_column))

                       

def _process_tautomer_rows(df, tautomerizer, smiles_column, mol_column, name_column):
    """
    Common function to process tautomerization for both single-core and multiprocessing.
    """
    results = []
    for _, row in df.iterrows():
        mol = row[mol_column]
        highlights = row.get('highlights', None)
        tautomers_smiles = tautomerizer.tautomerize(mol=mol, name = row[name_column])

        if len(tautomers_smiles) == 1:
            results.append({name_column: row[name_column],
                            mol_column: Chem.MolFromSmiles(tautomers_smiles[0]),
                            smiles_column: tautomers_smiles[0],
                            'highlights': highlights})
        else:
            two_digits = len(tautomers_smiles) >= 10
            for i, tautomer in enumerate(tautomers_smiles):
                results.append({name_column: f"{row[name_column]}_{i+1:02}" if two_digits else f"{row[name_column]}_{i+1}",
                                mol_column: Chem.MolFromSmiles(tautomer),
                                smiles_column: tautomer,
                                'highlights': highlights})

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
