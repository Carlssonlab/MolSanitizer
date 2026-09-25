from contextlib import nullcontext
import argparse
import logging
import platform
import shutil
import multiprocessing as mp
from io import StringIO

from pandas import DataFrame, read_csv, concat
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from msani.moltransform.neutralizer import Neutralizer
from msani.io.parsers import CustomHelpFormatter

# forkserver avoids inheriting parent locks (RDKit allocator, logging) on Linux.
# Windows only supports 'spawn'; macOS prefers 'spawn' too (fork is deprecated).
_MP_START_METHOD = 'forkserver' if platform.system() == 'Linux' else 'spawn'

RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('msani')
PROTONATION_RULES_PATH = Path(__file__).parent.parent / 'Data' / 'ionizations_v3.txt'
PHOSPHATE_2STAGE_PATTERN = Chem.MolFromSmarts(
    '[OH1&+0;$(O-P(=O)(-[OH])-[#6&+0,#8&+0])]'
)

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

    def _recursive_reaction_mols(self, mol, rxn, collection, visited=None):
        """Recursively apply one rule, retaining terminal product molecules."""
        if visited is None:
            visited = set()

        mol_smiles = Chem.MolToSmiles(mol)
        if mol_smiles in visited:
            return collection

        visited.add(mol_smiles)
        reactive = False
        outcomes = rxn.RunReactants((mol,))
        if outcomes:
            reactive = True
            product = outcomes[0][0]
            try:
                error = Chem.SanitizeMol(product, catchErrors=True)
                if error == 0:
                    self._recursive_reaction_mols(
                        product, rxn, collection, visited
                    )
                else:
                    if self.debug:
                        print(
                            "Sanitization error for molecule: "
                            f"{Chem.MolToSmiles(product)}"
                        )
                    reactive = False
            except Exception as e:
                logger.info(
                    f"Error sanitizing molecule: {Chem.MolToSmiles(product)}. "
                    f"Exception: {e}"
                )

        if not reactive:
            collection[mol_smiles] = mol
        return collection

    def recursive_reaction(self, mol, rxn, collection, visited=None):
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
        products = self._recursive_reaction_mols(mol, rxn, {}, visited)
        collection.update(products)
        return collection

    def _process_reaction_outcome_mols(
            self, outcomes, rxn, name, original_smi, original_mol):
        """Sanitize reaction outcomes and retain their product molecules."""
        products = {}
        selected_outcomes = (
            outcomes if name in self.enumerating_rules else outcomes[:1]
        )
        for outcome in selected_outcomes:
            products.update(self._sanitize_and_process_product_mols(
                outcome[0], rxn, original_smi, original_mol
            ))
        return products

    def _sanitize_and_process_product_mols(
            self, product, rxn, fallback_smi, fallback_mol):
        """Return terminal product molecules or the unchanged fallback molecule."""
        try:
            error = Chem.SanitizeMol(product, catchErrors=True)
            if error == 0:
                return self._recursive_reaction_mols(product, rxn, {})
            if self.debug:
                print(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}")
            return {fallback_smi: fallback_mol}
        except Exception as e:
            logger.error(f"Exception during sanitization: {e}")
            if self.debug:
                print(f"Exception during sanitization: {e}")
            return {fallback_smi: fallback_mol}
    
    def _process_reaction_outcomes(self, outcomes, rxn, name, original_smi):
        """
        Process reaction outcomes and handle sanitization errors.
        
        Returns:
            set: Set of valid product SMILES.
        """
        original_mol = Chem.MolFromSmiles(original_smi)
        return set(self._process_reaction_outcome_mols(
            outcomes, rxn, name, original_smi, original_mol
        ))

    def _sanitize_and_process_product(self, product, rxn, fallback_smi):
        """
        Sanitize a product molecule and apply recursive reactions.
        
        Returns:
            set: Set containing either the processed products or fallback SMILES.
        """
        fallback_mol = Chem.MolFromSmiles(fallback_smi)
        return set(self._sanitize_and_process_product_mols(
            product, rxn, fallback_smi, fallback_mol
        ))

    def apply_reactions_mols(self, input_mol, reactions):
        """Apply ordered rules using canonical keys mapped to product molecules."""
        initial_mol = Chem.Mol(input_mol)
        initial_smiles = Chem.MolToSmiles(initial_mol)
        current_mols = {initial_smiles: initial_mol}
        try:
            for rxn, name, additional_info in reactions:
                next_mols = {}
                for smi, mol in current_mols.items():
                    if mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                        outcomes = rxn.RunReactants((mol,))
                        if self.debug:
                            print(
                                f"\tApplying reaction {name} "
                                f"{additional_info} to {smi}"
                            )
                        next_mols.update(self._process_reaction_outcome_mols(
                            outcomes, rxn, name, smi, mol
                        ))
                    else:
                        next_mols[smi] = mol
                current_mols = next_mols
        except Exception as e:
            logger.error(f"Error applying reactions: {str(e)}")
            if self.debug:
                print(f"Error applying reactions: {str(e)}")
            return {initial_smiles: initial_mol}
        return current_mols

    def apply_reactions(self, input_mol, reactions):
        """
        Apply a list of reactions to a list of molecules and return the resulting products.
        
        Parameters:
            mol_list (list): A list of RDKit molecule objects.
            reactions (list): A list of reactions in the form [(rdkit.Chem.rdChem.Reaction object, name, additional information), ...].
            
        Returns:
            set: A set of unique products (in SMILES) generated from the reaction.
        """
        return set(self.apply_reactions_mols(input_mol, reactions))
    
    def _cache_rule_matches(self, mol: Chem.Mol) -> tuple[dict[int, bool], bool]:
        """Match each unique pH rule once against the original molecule."""
        rule_matches = {}
        matched_heteroacid_above_ph6 = False

        for pH in self.pH_values:
            for reaction, name, _ in self.rules_across_pH[pH]:
                reaction_key = id(reaction)
                if reaction_key not in rule_matches:
                    rule_matches[reaction_key] = mol.HasSubstructMatch(
                        reaction.GetReactantTemplate(0)
                    )
                if (name == 'heteroacid' and pH > 6.0 and
                        rule_matches[reaction_key]):
                    matched_heteroacid_above_ph6 = True

        phosphate_match = (
            matched_heteroacid_above_ph6 and
            mol.HasSubstructMatch(PHOSPHATE_2STAGE_PATTERN)
        )
        return rule_matches, phosphate_match

    def check_duplicated_rule_combs(
            self,
            mol: Chem.Mol,
            pH: float,
            rule_matches: dict[int, bool] | None = None,
            phosphate_match: bool | None = None) -> tuple:
        """
        Check if the molecule matches any of the rules for a given pH
        and return the rule combinations.
        """
        rule_combinations = []
        for rule in self.rules_across_pH[pH]:
            reaction_key = id(rule[0])
            matches = (
                rule_matches[reaction_key]
                if rule_matches is not None
                else mol.HasSubstructMatch(rule[0].GetReactantTemplate(0))
            )
            if matches:
                rule_combinations.append(rule[1])
                if rule[1] == 'heteroacid' and pH > 6.0 and (
                    phosphate_match
                    if phosphate_match is not None
                    else mol.HasSubstructMatch(PHOSPHATE_2STAGE_PATTERN)
                ): # Special case for phosphates
                    rule_combinations.append('phosphate-2stage')
        return tuple(rule_combinations)
    
    def _ionize_mols(
            self,
            smiles: str = None,
            mol: Chem.Mol = None) -> dict[str, Chem.Mol]:
        """Ionize into canonical SMILES keys mapped to product molecules."""
        if smiles is None and mol is None:
            raise ValueError("Either 'smiles' or 'mol' must be provided.")
        if mol is None:
            mol = Chem.MolFromSmiles(smiles)
        else:
            # Reactions sanitize their products in place, so isolate all work
            # from a molecule owned by the caller or another DataFrame row.
            mol = Chem.Mol(mol)

        if self.neutralize:
            mol = Neutralizer.neutralize_mol(mol)

        variation_mols = {}
        applied_rule_combinations = set()
        rule_matches, phosphate_match = self._cache_rule_matches(mol)
        for pH in self.pH_values:
            rule_combinations = self.check_duplicated_rule_combs(
                mol, pH, rule_matches, phosphate_match
            )
            if rule_combinations in applied_rule_combinations:
                if self.debug: print(f'Skipping pH {round(pH, 1)} due to the duplicated rule applied.')
                continue
            applied_rule_combinations.add(rule_combinations)
            if self.debug: print('Processing pH:', round(pH, 1))
            variation_mols.update(self.apply_reactions_mols(
                mol, self.rules_across_pH[pH]
            ))
        return variation_mols

    def ionize(self, smiles: str = None, mol: Chem.Mol = None) -> list:
        """Protonate the input molecule and return canonical product SMILES."""
        return list(self._ionize_mols(smiles=smiles, mol=mol))



 
    def ionize_df_mp(self,
                     df: DataFrame,
                     smiles_column: str = 'smiles',
                     name_column: str = 'ids',
                     mol_column: str = 'mol',
                     *, pool=None) -> DataFrame:
        """
        Protonate molecules using an existing initialized pool or a local pool.
        Local pools use forkserver on Linux and spawn elsewhere. The initializer
        installs the rules once per worker rather than once per input row.

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

        numcores = min(self.numcores, len(df))
        ctx = mp.get_context(_MP_START_METHOD)
        # Package the row data efficiently using python dicts to avoid memory bloat of large Pandas Series
        keys = df.columns.tolist()
        rows = ((dict(zip(keys, r)), smiles_column, name_column, mol_column) for r in df.itertuples(index=False, name=None))
        
        results = []
        pool_context = nullcontext(pool) if pool is not None else ctx.Pool(
            processes=numcores,
            initializer=_init_ionizer_worker,
            initargs=(self,))
        with pool_context as pool:
            # Ordered results keep the output and its id numbering reproducible.
            for res in pool.imap(_process_single_ionization_row, rows, chunksize=10):
                if res:
                    results.extend(res)
                    
        return DataFrame(results)
    
    def ionize_df(self, 
                  df: DataFrame,
                  smiles_column: str = 'smiles',
                  name_column: str = 'ids',
                  mol_column: str = 'mol',
                  *, pool=None, parallel=True) -> DataFrame:
        """
        Protonate the input molecules using multiprocessing or single core based on `numcores`.
        """
        if len(df) == 0:
            return df
        if mol_column not in df.columns:
            df = df.copy() # Make a shallow copy to safely add mol_column
            df[mol_column] = df[smiles_column].apply(lambda x: Chem.MolFromSmiles(x))
            
        if parallel and self.numcores > 1:
            return self.ionize_df_mp(df, smiles_column, name_column, mol_column, pool=pool)
        else:
            keys = df.columns.tolist()
            rows = ((dict(zip(keys, r)), smiles_column, name_column, mol_column) for r in df.itertuples(index=False, name=None))
            results = []

            for row_args in rows:
                res = _process_single_ionization_row(row_args, worker=self)
                if res:
                    results.extend(res)
            return DataFrame(results)

# ---------------------------------------------------------------------------
# Module-level worker helpers for multiprocessing
# ---------------------------------------------------------------------------

_ionizer_worker = None

def _init_ionizer_worker(ionizer):
    """Initializer run once per worker process to prevent memory spikes from pickling."""
    global _ionizer_worker
    _ionizer_worker = ionizer

def _process_single_ionization_row(args, worker=None):
    """
    Worker function to process ionization for a single row.
    Caught exceptions simply log an error and return an empty list so the job continues.
    """
    row, smiles_column, name_column, mol_column = args
    worker = _ionizer_worker if worker is None else worker
    results = []
    
    try:
        longname = row.get('longname', None)
        original_idx = row.get('original_idx', None)
        
        protonated_mols = worker._ionize_mols(mol=row[mol_column])

        for smiles, protonated_mol in protonated_mols.items():
            results.append({
                name_column: row[name_column],
                mol_column: protonated_mol,
                smiles_column: smiles,
                'longname': longname,
                'original_idx': original_idx
            })
    except Exception as e:
        logger.error(f"Error during ionization for {row[name_column]}: {str(e)}")
        
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
