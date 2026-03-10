import argparse
import logging
import shutil
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from io import StringIO

from numpy import array_split, arange
from pandas import DataFrame, read_csv, concat
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from msani.moltransform.neutralizer import Neutralizer
from msani.io.parsers import CustomHelpFormatter

RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('msani')
PROTONATION_RULES_PATH = Path(__file__).parent.parent / 'Data' / 'ionizations_v3.txt'

logo=r""" _____            _              
|_   _|          (_)             
  | |  ___  _ __  _ _______ _ __ 
  | | / _ \| '_ \| |_  / _ \ '__|
 _| || (_) | | | | |/ /  __/ |   
 \___/\___/|_| |_|_/___\___|_|   
                                 
                From MolSanitizer
"""

    
class Ionizer:
    """
        A class to protonate the input molecules using a set of predefined reactions.

        Parameters
        ----------

        pH: int, default 7.
            The pH value to use as a reference for protonation.

        pH_range: int, default 0.
            The range of pH values to consider for protonation. If 0, only the pH value is considered, else a range of pH values are considered.

        numcores: int, default 1.
            Number of cores to use for multiprocessing. If 1, no multiprocessing is used.
        
        neutralize: bool, default False.
            Neutralize the input molecules before protonation.

        debug: bool, default False.
            Print debug message

        Examples
        ----------
        
        >>> from msani.moltransform.ionizer import Ionizer
        >>> ionizer = Ionizer(pH = 7, pH_range = 2)

        ionize a single molecule from SMILES string:

        >>> results = ionizer.ionize(smiles = 'CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1')

        ionize a DataFrame of molecules with SMILES strings:

        >>> df = ionizer.ionize_df(mol_df)
        """
    def __init__(self,
                 smartsFile = PROTONATION_RULES_PATH,
                 pH: int = 7,
                 pH_range: int = 0,
                 numcores: int = 1,
                 neutralize: bool = True,
                 debug = False):
        
        self.pH = pH
        self.pH_range = pH_range
        self.numcores = numcores
        self.debug = debug
        if smartsFile is None: smartsFile = PROTONATION_RULES_PATH
        self.rules = self.load_protonation_rules(smartsFile)
        self.enumerating_rules = self.rules[self.rules['Enumerate'] == 1]['FUNCTIONAL_GROUP'].to_list()
        self.neutralize = neutralize
        if self.debug: 
            print(f'Loaded {len(self.rules)} protonation rules. Of these, {len(self.enumerating_rules)} rules are enumerating')


        self.pH_values = (self.pH-0.1,
                          self.pH,
                          self.pH+0.1) \
        if self.pH_range == 0\
            else (max(0, self.pH - self.pH_range - 0.1),
                  max(0, self.pH - self.pH_range),
                  max(0, self.pH - self.pH_range + 0.1),
                  self.pH - 0.1,
                  self.pH,
                  self.pH + 0.1,
                  min(14, self.pH + self.pH_range - 0.1),
                  min(14, self.pH + self.pH_range), 
                  min(14, self.pH + self.pH_range + 0.1))

        self.rules_across_pH = {}
        for pH in self.pH_values:
            rules = self.extract_necessary_rules(pH)
            self.rules_across_pH[pH] = rules

    def __repr__(self):
        cls_name = self.__class__.__name__
        attrs = ', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())
        return f'{cls_name}({attrs})'
        
    @staticmethod
    def load_protonation_rules(file_path: str):
        '''
        Load the protonation rules from a file containing SMARTS strings.
        Expected format from the text file: 

        =================  =====  ====  ==========  ====================  ===
        FUNCTIONAL_GROUP   pKa    TYPE  Enumerate   REACTION              REF
        =================  =====  ====  ==========  ====================  ===
        amine              10     BASE  1           [reagent]>>[product]
        =================  =====  ====  ==========  ====================  ===
        '''
        with open(file_path, 'r') as f:
            uncommented = [line for line in f if not (line.startswith('#')) and line.strip()]

        rules = read_csv(StringIO(''.join(uncommented)), sep=r"\s+", header=None,
                            names=['FUNCTIONAL_GROUP', 'pKa', 'TYPE', 'Enumerate', 'REACTION', 'REF'])
        
        rules.loc[:, 'Mol'] = rules['REACTION'].apply(lambda x: AllChem.ReactionFromSmarts(x))
        acid_rules = rules[rules['TYPE'] == 'ACID'].copy()
        base_rules = rules[rules['TYPE'] == 'BASE'].copy()
        acid_rules.sort_values(by=['pKa', 'Enumerate'], ascending=[True, True], inplace=True)
        base_rules.sort_values(by=['pKa', 'Enumerate'], ascending=[False, True], inplace=True)
        rules = concat([acid_rules, base_rules], ignore_index=True)
        return rules
    
    def extract_necessary_rules(self, pH: int):
        """
        Processes ionization rules from a file and filters them based on pH, pKa, and type.
        
        Parameters:
            pH (int): The pH value to filter the rules.
            
        Returns:
            DataFrame: Filtered rules based on the pH.
        """

        # Filter based on pH for ACID and BASE rules
        filtered_rules = self.rules[
            ((self.rules['TYPE'] == 'ACID') & (self.rules['pKa'] <= pH)) |
            ((self.rules['TYPE'] == 'BASE') & (self.rules['pKa'] >= pH))
        ]

        reaction_list = []

        for i, row in filtered_rules.iterrows():
            reaction_list.append((row['Mol'], row['FUNCTIONAL_GROUP'], f"({row['TYPE']} - est. pKa {row['pKa']})"))
        if self.debug:
            print(f'Parsed {len(reaction_list)} rules for pH {round(pH, 1)}')    
        return reaction_list

    def recursive_reaction(self, mol, rxn, collection, visited = None):
        """
        Recursively apply a reaction to a molecule and its products, collecting unique products.

        Args:
            mol (rdkit.Chem.rdchem.Mol): The reactant molecule.
            rxn (rdkit.Chem.rdChem.Reaction): The reaction to apply.
            collection (set): A set of unique products (in SMILES) generated from the reaction.
            visited (set): A set of SMILES strings for molecules already processed to avoid redundancy.
            
        Returns:
            set: A set of unique products (in SMILES) generated from the reaction.
        """
        if visited is None:
            visited = set()

        mol_smiles = Chem.MolToSmiles(mol)
        last_successful_smiles = mol_smiles  # Keep track of the last valid molecule
        
        # Check if the molecule has already been visited
        if mol_smiles in visited:
            return collection  # Skip redundant processing

        visited.add(mol_smiles)  # Mark the molecule as visited
        reactive = False
        outcomes = rxn.RunReactants((mol,))
        if outcomes:  # Check if there are any outcomes            
            reactive = True
            product = outcomes[0][0]
            try:
                error = Chem.SanitizeMol(product, catchErrors=True)
                if error == 0:
                    last_successful_smiles = Chem.MolToSmiles(product)
                    self.recursive_reaction(product, rxn, collection, visited)
                else:
                    if self.debug: print(f"Sanitization error for molecule: {Chem.MolToSmiles(product)}")
                    reactive = False
            except Exception as e:
                logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}. Exception: {e}")
        
        if not reactive: 
            collection.add(last_successful_smiles)  # Add the last valid molecule if it is not reactive
        return collection  # Return the collection of products
    
    def _process_reaction_outcomes(self, outcomes, rxn, name, original_smi):
        """
        Process reaction outcomes and handle sanitization errors.
        
        Returns:
            set: Set of valid product SMILES.
        """
        products = set()
        
        if name in self.enumerating_rules:
            # Process all outcomes for enumerating rules
            for outcome in outcomes:
                product_smiles = self._sanitize_and_process_product(
                    outcome[0], rxn, original_smi
                )
                products.update(product_smiles)
        else:
            # Process only first outcome for non-enumerating rules
            if outcomes:
                product_smiles = self._sanitize_and_process_product(
                    outcomes[0][0], rxn, original_smi
                )
                products.update(product_smiles)
        
        return products

    def _sanitize_and_process_product(self, product, rxn, fallback_smi):
        """
        Sanitize a product molecule and apply recursive reactions.
        
        Returns:
            set: Set containing either the processed products or fallback SMILES.
        """
        try:
            error = Chem.SanitizeMol(product, catchErrors=True)
            if error == 0:
                return self.recursive_reaction(product, rxn, set())
            else:
                if self.debug: 
                    print(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}")
                return {fallback_smi}
        except Exception as e:
            logger.error(f"Exception during sanitization: {e}")
            if self.debug:
                print(f"Exception during sanitization: {e}")
            return {fallback_smi}
    
    def apply_reactions(self, input_mol, reactions):
        """
        Apply a list of reactions to a list of molecules and return the resulting products.
        
        Parameters:
            mol_list (list): A list of RDKit molecule objects.
            reactions (list): A list of reactions in the form [(rdkit.Chem.rdChem.Reaction object, name, additional information), ...].
            
        Returns:
            set: A set of unique products (in SMILES) generated from the reaction.
        """
        current_smiles = {Chem.MolToSmiles(input_mol)} # Start with the initial molecule
        try:
            for rxn, name, additional_info in reactions:
                next_smiles = set()  # Reset output for the next reaction
                for smi in current_smiles:
                    mol = Chem.MolFromSmiles(smi)
                    if mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                        outcomes = rxn.RunReactants((mol,))
                        if self.debug: print(f"\tApplying reaction {name} {additional_info} to {smi}")
                        products = self._process_reaction_outcomes(outcomes, rxn, name, smi)
                        next_smiles.update(products)
                    else:
                        next_smiles.add(smi)  # No match, keep original
                current_smiles = next_smiles
        except Exception as e:
            logger.error(f"Error applying reactions: {str(e)}")
            if self.debug: print(f"Error applying reactions: {str(e)}")
            return [Chem.MolToSmiles(input_mol)]
        return current_smiles
    
    def check_duplicated_rule_combs(self, mol: Chem.Mol, pH: float) -> tuple:
        """
        Check if the molecule matches any of the rules for a given pH
        and return the rule combinations.
        """
        rule_combinations = []
        for rule in self.rules_across_pH[pH]:
            if mol.HasSubstructMatch(rule[0].GetReactantTemplate(0)):
                rule_combinations.append(rule[1])
                if rule[1] == 'heteroacid' and pH > 6.0 and \
                    mol.HasSubstructMatch(Chem.MolFromSmarts('[OH1&+0;$(O-P(=O)(-[OH])-[#6&+0,#8&+0])]')): # Special case for phosphates
                        rule_combinations.append('phosphate-2stage')
        return tuple(rule_combinations)
    
    def ionize(self, smiles: str = None, mol: Chem.Mol = None) -> list:
        """
        Protonate the input molecule using a set of predefined reactions.
        
        Parameters:

            smiles (str): The SMILES string of the molecule.
            
            mol (Chem.Mol): The RDKit molecule object.
            
        Returns:

            list: A list of protonated SMILES strings.
        """
        if smiles is None and mol is None:
            raise ValueError("Either 'smiles' or 'mol' must be provided.")
        if mol is None:
            mol = Chem.MolFromSmiles(smiles)

        if self.neutralize:
            mol = Neutralizer.neutralize_mol(mol)

        variation_sets = set()
        applied_rule_combinations = set()
        for pH in self.pH_values:
            rule_combinations = self.check_duplicated_rule_combs(mol, pH)
            if rule_combinations in applied_rule_combinations:
                if self.debug: print(f'Skipping pH {round(pH, 1)} due to the duplicated rule applied.')
                continue
            applied_rule_combinations.add(rule_combinations)
            if self.debug: print('Processing pH:', round(pH, 1))
            # variations = list(self.recursive_reaction(mol, self.rules_across_pH[pH], set()))
            variations = list(self.apply_reactions(mol, self.rules_across_pH[pH]))
            variation_sets.update(variations)
        return list(variation_sets)



 
    def ionize_df_mp(self,
                     df: DataFrame,
                     smiles_column: str = 'smiles',
                     name_column: str = 'ids',
                     mol_column: str = 'mol') -> DataFrame:
        """
        Protonate the input molecules using ProcessPoolExecutor with
        forkserver start method and an initializer to avoid:
          1. Pipe-buffer deadlocks from imap_unordered with large results.
          2. RDKit/logging lock inheritance via os.fork() on Linux.
          3. Re-pickling the Ionizer object for every chunk.

        Parameters:
            df (DataFrame): The input DataFrame containing SMILES strings.
            smiles_column (str): The column name containing SMILES strings (default: smiles).
            name_column (str): The column name for molecule identifiers (default: ids).
            mol_column (str): The column name for RDKit molecule objects (default: mol).

        Returns:
            DataFrame: A DataFrame with protonated molecules, including SMILES and identifiers.
        """
        # Ensure mol_column exists
        if mol_column not in df.columns:
            df = df.copy()  # Make a single shallow copy to safely add mol_column
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        # Determine number of cores
        numcores = min(self.numcores, len(df))  # Prevent using more cores than data chunks

        # Handle contiguous chunks using numpy array_split
        chunk_indices = array_split(arange(len(df)), numcores)
        chunks = [df.iloc[indices] for indices in chunk_indices]

        # forkserver: worker starts in a clean process (no inherited locks from parent).
        # initializer: Ionizer is pickled once at worker startup, not once per chunk.
        ctx = mp.get_context('forkserver')
        with ProcessPoolExecutor(
                max_workers=numcores,
                mp_context=ctx,
                initializer=_init_ionizer_worker,
                initargs=(self,)) as executor:
            # list() drains eagerly via an internal thread — safe from pipe-block.
            batch_results = list(executor.map(_ionize_chunk_worker, chunks))

        results = [row for batch in batch_results if batch for row in batch]
        return DataFrame(results)
    
    def ionize_df(self, 
                  df: DataFrame,
                  smiles_column: str = 'smiles',
                  name_column: str = 'ids',
                  mol_column: str = 'mol') -> DataFrame: 
        """
        Protonate the input molecules using multiprocessing or single core based on `numcores`.
        
        Parameters:
            df (DataFrame): The input DataFrame containing SMILES strings.
            smiles_column (str): The column name containing SMILES strings (default: smiles).
            name_column (str): The column name for molecule identifiers (default: ids).
            mol_column (str): The column name for RDKit molecule objects (default: mol).

        Returns:
            DataFrame: A DataFrame with protonated molecules, including SMILES and identifiers.
        """
        if len(df) == 0:
            return df
        if mol_column not in df.columns:
            df = df.copy() # Make a shallow copy to safely add mol_column
            df[mol_column] = df[smiles_column].apply(lambda x: Chem.MolFromSmiles(x))
            
        if self.numcores > 1:
            return self.ionize_df_mp(df, smiles_column, name_column, mol_column)
        else:
            return DataFrame(_process_ionization_rows(df, self, smiles_column, mol_column, name_column))


# ---------------------------------------------------------------------------
# Module-level worker helpers for ProcessPoolExecutor
# These must be at module level (not nested) to be picklable.
# ---------------------------------------------------------------------------

_ionizer_worker = None  # per-worker singleton set by the initializer


def _init_ionizer_worker(ionizer):
    """Initializer run once per worker process."""
    global _ionizer_worker
    _ionizer_worker = ionizer


def _ionize_chunk_worker(chunk):
    """Top-level worker function dispatched by ProcessPoolExecutor."""
    return _process_ionization_rows(chunk, _ionizer_worker, 'smiles', 'mol', 'ids')


def _process_ionization_rows(df, ionizer, smiles_column, mol_column, name_column):
    """
    Common function to process ionization for both single-core and multiprocessing.
    """
    results = []
    for _, row in df.iterrows():
        longname = row.get('longname', None)
        original_idx = row.get('original_idx', None)  # For debugging purposes
        protonated_smiles = ionizer.ionize(mol=row[mol_column])

        for _, smiles in enumerate(protonated_smiles):
            results.append({name_column: row[name_column],
                            mol_column: Chem.MolFromSmiles(smiles),
                            smiles_column: smiles,
                            'longname': longname,
                            'original_idx': original_idx})

    return results

def main():
    print(logo)
    parser = argparse.ArgumentParser(description='Protonate a molecule or a file of SMILES', formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input', default=None, type=str, help='Input file containing SMILES strings and names')
    parser.add_argument('-l', '--library', default=PROTONATION_RULES_PATH, type=str, help='Protonation rules file (default: msani/Data/ionizations_v3.txt)')
    parser.add_argument('-t', '--template', action='store_true', help='Create a template for protonation rules')
    parser.add_argument('-s', '--smiles', default=None, type=str, help='SMILES string of the molecule')
    parser.add_argument('-p', '--pH', default=7, type=int, help='pH value to use for ionization')
    parser.add_argument('-r', '--pH_range', default=0, type=int, help='Range of pH values to consider for ionization')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('-j', '--numcores', default=4, type=int, help='Number of cores to use for multiprocessing')
    parser.add_argument('-o', '--output', default='protonated_molecules.smi', type=str, help='Output file to save protonated molecules')

    args = parser.parse_args()
    if args.smiles and args.input:
        raise ValueError("Either SMILES or input file must be provided.")
    
    if args.template:
        # Copy the ionization rules file to the current directory
        output_filename = f"ionization_template_msani.txt"
        shutil.copy(PROTONATION_RULES_PATH, output_filename)
        print(f"Successfully create a template for protonation to {output_filename}")
        return  # Exit after creating the library
    
    ionizer = Ionizer(pH = args.pH, 
                      smartsFile = args.library,
                      pH_range = args.pH_range, 
                      numcores = args.numcores,
                      debug = args.debug)
    if args.smiles:
        results = ionizer.ionize(smiles = args.smiles)
        print(f"Protonated SMILES:")
        for result in results:
            print(result)
    else:
        df = read_csv(args.input, sep =r'\s+', header=None, names=['smiles', 'ids'])
        ionized_df = ionizer.ionize_df(df)
        ionized_df[['smiles', 'ids']].to_csv(args.output, index=False, header=False, sep =' ')

if __name__ == '__main__':
    main()
