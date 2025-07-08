import argparse
import pandas as pd
import logging
import multiprocessing as mp

from io import StringIO
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from functools import partial

from .neutralizer import Neutralizer

RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('msani')
PROTONATION_RULES_PATH = Path(__file__).parent.parent / 'Data' / 'ionizations_v3.txt'

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

        Example use:
        ----------
        
        >>> from msani.ionizer import Ionizer\n
        >>> ionizer = Ionizer(pH = 7, pH_range = 2)\n
        >>> results = ionizer.ionize(smiles = 'CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1')\n
        >>> df = ionizer.ionize_df(mol_df, pH=7, pH_range=0, num_cores=1, debug=False)\n
        """
    def __init__(self,
                 smartsFile = PROTONATION_RULES_PATH,
                 pH: int = 7,
                 pH_range: int = 0,
                 num_cores: int = 1,
                 neutralize: bool = True,
                 debug = False):
        self.pH = pH
        self.pH_range = pH_range
        self.num_cores = num_cores
        self.debug = debug
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
            FUNCTIONAL_GROUP	pKa	    TYPE	Enumerate   REACTION                REF
            amine	            10	    BASE	1           [reagent]>>[product]
        '''
        with open(file_path, 'r') as f:
            uncommented = [line for line in f if not (line.startswith('#')) and line.strip()]

        rules = pd.read_csv(StringIO(''.join(uncommented)), sep=r"\s+", header=None,
                            names=['FUNCTIONAL_GROUP', 'pKa', 'TYPE', 'Enumerate', 'REACTION', 'REF'])
        
        rules['Mol'] = rules['REACTION'].apply(lambda x: AllChem.ReactionFromSmarts(x))
        acid_rules = rules[rules['TYPE'] == 'ACID'].copy()
        base_rules = rules[rules['TYPE'] == 'BASE'].copy()
        acid_rules.sort_values(by=['pKa', 'Enumerate'], ascending=[True, False], inplace=True)
        base_rules.sort_values(by=['pKa', 'Enumerate'], ascending=[False, False], inplace=True)
        rules = pd.concat([acid_rules, base_rules], ignore_index=True)
        return rules
    
    def extract_necessary_rules(self, pH: int):
        """
        Processes ionization rules from a file and filters them based on pH, pKa, and type.
        
        Parameters:
            pH (int): The pH value to filter the rules.
            
        Returns:
            pd.DataFrame: Filtered rules based on the pH.
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

    def apply_reactions(self, input_mol, reactions):
        """
        Apply a list of reactions to a list of molecules and return the resulting products.
        
        Parameters:
            mol_list (list): A list of RDKit molecule objects.
            reactions (list): A list of reactions in the form [(rdkit.Chem.rdChem.Reaction object, name, additional information), ...].
            
        Returns:
            set: A set of unique products (in SMILES) generated from the reaction.
        """
        output = [Chem.MolToSmiles(input_mol)]  # Start with the initial molecule
        for rxn, name, additional_information in reactions:
            input = list(output).copy()
            output = set()  # Reset output for the next reaction
            for smi in input:
                mol = Chem.MolFromSmiles(smi)
                if mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                    outcomes = rxn.RunReactants((mol,))
                    if self.debug:
                        print(f"\tApplying reaction {name} {additional_information} to {Chem.MolToSmiles(mol)}")
                    if name in self.enumerating_rules:
                        for outcome in outcomes:
                            product = outcome[0]
                            try:
                                error = Chem.SanitizeMol(product, catchErrors=True)
                                if error == 0: output.update(self.recursive_reaction(product, rxn, set()))
                            except Exception as e:
                                logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}. Exception: {e}")
                    else:
                        product = outcomes[0][0]
                        try:
                            error = Chem.SanitizeMol(product, catchErrors=True)
                            if error == 0: output.update(self.recursive_reaction(product, rxn, set()))
                        except Exception as e:
                            logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}. Exception: {e}")
                else:
                    output = input.copy()
        # print(output)
        return output
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
                     df: pd.DataFrame,
                     smiles_column: str = 'smiles',
                     name_column: str = 'ids',
                     mol_column: str = 'mol') -> pd.DataFrame:
        """
        Protonate the input molecules using multiprocessing with chunked DataFrame processing.
        """
        # Ensure mol_column exists
        if mol_column not in df.columns:
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        # Determine number of cores
        num_cores = min(self.num_cores, len(df))  # Prevent using more cores than data chunks

        # Split DataFrame into chunks
        chunks = [df.iloc[i::num_cores] for i in range(num_cores)]

        # Create a partial function for multiprocessing
        process_func = partial(_process_ionization_rows, ionizer=self,
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
                    logger.error(f"Error processing a molecule batch: {str(e)}")

        return pd.DataFrame(results)
    
    def ionize_df(self, 
                  df: pd.DataFrame,
                  smiles_column: str = 'smiles',
                  name_column: str = 'ids',
                  mol_column: str = 'mol') -> pd.DataFrame:
        """
        Protonate the input molecules using multiprocessing or single core based on `num_cores`.
        """
        if len(df) == 0:
            return df
        if mol_column not in df.columns:
            df[mol_column] = df[smiles_column].apply(lambda x: Chem.MolFromSmiles(x))
        if self.num_cores > 1:
            return self.ionize_df_mp(df, smiles_column, name_column, mol_column)
        else:
            return pd.DataFrame(_process_ionization_rows(df, self, smiles_column, mol_column, name_column))


def _process_ionization_rows(df, ionizer, smiles_column, mol_column, name_column):
    """
    Common function to process ionization for both single-core and multiprocessing.
    """
    results = []
    for _, row in df.iterrows():
        highlights = row.get('highlights', None)
        protonated_smiles = ionizer.ionize(mol=row[mol_column])

        if len(protonated_smiles) == 1:
            results.append({name_column: row[name_column],
                            mol_column: Chem.MolFromSmiles(protonated_smiles[0]),
                            smiles_column: protonated_smiles[0],
                            'highlights': highlights})
        else:
            two_digits = len(protonated_smiles) >= 10
            for i, smile in enumerate(protonated_smiles):
                results.append({name_column: f"{row[name_column]}_{i+1:02}" if two_digits else f"{row[name_column]}_{i+1}",
                                mol_column: Chem.MolFromSmiles(smile),
                                smiles_column: smile,
                                'highlights': highlights})
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Protonate a molecule or a file of SMILES')
    parser.add_argument('-i', '--input', default=None, type=str, help='Input file containing SMILES strings and names')
    parser.add_argument('-s', '--smiles', default=None, type=str, help='SMILES string of the molecule')
    parser.add_argument('-p', '--pH', default=7, type=int, help='pH value to use for ionization')
    parser.add_argument('-r', '--pH_range', default=0, type=int, help='Range of pH values to consider for ionization')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('-j', '--num_cores', default=4, type=int, help='Number of cores to use for multiprocessing')
    parser.add_argument('-o', '--output', default='protonated_molecules.smi', type=str, help='Output file to save protonated molecules')
    args = parser.parse_args()
    if args.smiles and args.input:
        raise ValueError("Either SMILES or input file must be provided.")

    ionizer = Ionizer(pH = args.pH, 
                      pH_range = args.pH_range, 
                      num_cores = args.num_cores,
                      debug = args.debug)
    if args.smiles:
        results = ionizer.ionize(smiles = args.smiles)
        for result in results:
            print(result)
    else:
        df = pd.read_csv(args.input, sep =r'\s+', header=None, names=['smiles', 'ids'])
        ionized_df = ionizer.ionize_df(df)
        ionized_df[['smiles', 'ids']].to_csv(args.output, index=False, header=False, sep =' ')

if __name__ == '__main__':
    main()
