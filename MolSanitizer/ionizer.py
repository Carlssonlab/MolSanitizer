import argparse
from pathlib import Path
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import logging
import sys
import multiprocessing as mp
from functools import partial


RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('molsani')
PROTONATION_RULES_PATH = Path(__file__).parent / 'Data' / 'ionizations_v3.txt'

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
        
        >>> from MolSanitizer.ionizer import Ionizer\n
        >>> ionizer = Ionizer(pH = 7, pH_range = 2)\n
        >>> results = ionizer.ionize(smiles = 'CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1')\n
        >>> df = ionizer.ionize_df(mol_df, pH=7, pH_range=0, num_cores=1, debug=False)\n
        """
    def __init__(self,
                 pH: int = 7,
                 pH_range: int = 0,
                 num_cores: int = 1,
                 neutralize: bool = True,
                 debug = False):
        self.pH = pH
        self.pH_range = pH_range
        self.num_cores = num_cores
        self.debug = debug
        self.rules = self.load_protonation_rules(PROTONATION_RULES_PATH)
        self.enumerating_rules = self.rules[self.rules['Enumerate'] == 1]['FUNCTIONAL_GROUP'].to_list()
        self.neutralize = neutralize
        if self.debug: 
            print(f'Loaded {len(self.rules)} protonation rules')
            print(f'In these, {len(self.enumerating_rules)} rules are enumerating')


        self.pH_values = (self.pH,) if self.pH_range == 0\
            else (max(0, self.pH - self.pH_range), self.pH, min(14, self.pH + self.pH_range))

        self.rules_across_pH = {}
        for pH in self.pH_values:
            rules = self.extract_necessary_rules(pH)
            self.rules_across_pH[pH] = rules
    
    @staticmethod
    def load_protonation_rules(file_path: str):
        '''
        Load the protonation rules from a file containing SMARTS strings.
        Expected format from the text file: 
            FUNCTIONAL_GROUP	pKa	    TYPE	Enumerate   REACTION
            amine	            10	    BASE	1           [reagent]>>[product]
        '''

        rules = pd.read_csv(file_path, sep=r"\s+")
        rules['Mol'] = rules['REACTION'].apply(lambda x: AllChem.ReactionFromSmarts(x))
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
            ((self.rules['TYPE'] == 'ACID') & (self.rules['pKa'] < pH)) |
            ((self.rules['TYPE'] == 'BASE') & (self.rules['pKa'] > pH))
        ]

        reaction_list = []

        for i, row in filtered_rules.iterrows():
            reaction_list.append((row['Mol'], row['FUNCTIONAL_GROUP']))
        if self.debug:
            print(f'Parsed {len(reaction_list)} rules for pH {pH}')    
        return reaction_list

    @staticmethod
    def neutralize_mol(mol: Chem.Mol) -> Chem.Mol:
        """Neutralize the input molecule by balancing the charges on atoms.
        Adapted from RDKit Cookbook: https://rdkit.org/docs/Cookbook.html
        Args:

            mol (rdkit.Chem.rdchem.Mol): The input molecule.
        Returns:

            rdkit.Chem.rdchem.Mol: The neutralized molecule.
        """

        pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
            smiles = Chem.MolToSmiles(mol)
            return Chem.MolFromSmiles(smiles)
        else: return mol

    def recursive_reaction(self, mol, reactions, collection, visited=None, debug=False):
        """
        Recursively apply reactions to the molecule until no new products are generated.

        Args:
            mol (rdkit.Chem.rdchem.Mol): The reactant molecule.
            reactions (list): A list of reactions in the form [(rdkit.Chem.rdChem.Reaction object, name), ...].
            collection (set): A set of unique products (in SMILES) generated from the reaction.
            visited (set): A set of SMILES strings for molecules already processed to avoid redundancy.
            debug (bool): Enable debug logging.

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
        for rxn, name in reactions:
            outcomes = rxn.RunReactants((mol,))
            if outcomes:  # Check if there are any outcomes            
                reactive = True
                if debug:
                    print(f"\tApplying reaction: {name} to {Chem.MolToSmiles(mol)}")
                if name in self.enumerating_rules:
                    for outcome in outcomes:
                        product = outcome[0]
                        try:
                            error = Chem.SanitizeMol(product, catchErrors=True)
                            if error == 0:
                                last_successful_smiles = Chem.MolToSmiles(product)
                                self.recursive_reaction(product, reactions, collection, visited, debug)
                            else:
                                if debug: print(f"Sanitization error for molecule: {Chem.MolToSmiles(product)}")
                                reactive = False
                                continue
                        except Exception as e:
                            logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}. Exception: {e}")
                else:
                    product = outcomes[0][0]
                    try:
                        error = Chem.SanitizeMol(product, catchErrors=True)
                        if error == 0:
                            last_successful_smiles = Chem.MolToSmiles(product)
                            self.recursive_reaction(product, reactions, collection, visited, debug)
                        else:
                            if debug: print(f"Sanitization error for molecule: {Chem.MolToSmiles(product)}")
                            reactive = False
                            continue
                    except Exception as e:
                        logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}. Exception: {e}")
        
        if not reactive: 
            collection.add(last_successful_smiles)  # Add the last valid molecule if it is not reactive
        return collection  # Return the collection of products

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
            mol = Ionizer.neutralize_mol(mol)


        variation_sets = set()
        for pH in self.pH_values:
            variations = list(self.recursive_reaction(mol, self.rules_across_pH[pH], set()))
            variation_sets.update(variations)
        return list(variation_sets)



    def ionize_df(self,
                df: pd.DataFrame,
                smiles_column: str = 'smiles',
                name_column: str = 'ids',
                mol_column: str = 'mol') -> pd.DataFrame:
        """
        Protonate the input molecules using multiprocessing.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            smiles_column (str): The name of the column containing SMILES strings.
            mol_column (str): The name of the column containing RDKit molecule objects.
            name_column (str): The name of the column containing molecule names.

        Returns:
            pd.DataFrame: The expanded DataFrame with protonated molecules.
        """
        # Ensure mol_column exists
        if mol_column not in df.columns:
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        # Create a partial function for multiprocessing
        process_func = partial(_process_ionization, ionizer=self,
                            smiles_column=smiles_column,
                            mol_column=mol_column,
                            name_column=name_column)

        results = []

        # Use multiprocessing Pool
        num_cores = self.num_cores if self.num_cores > 1 else 1  # Ensure at least 1 core
        with mp.Pool(processes=num_cores) as pool:
            # Submit tasks asynchronously
            async_results = [pool.apply_async(process_func, (row,)) for _, row in df.iterrows()]

            # Collect results
            for async_result in async_results:
                try:
                    results.extend(async_result.get(timeout=60))  # Get result with a timeout
                except mp.TimeoutError:
                    logger.warning("Timeout occurred while ionizing a molecule. Skipping.")
                except Exception as e:
                    logger.error(f"Error processing a molecule: {str(e)}")

        return pd.DataFrame(results)



def _process_ionization(row, ionizer, smiles_column, mol_column, name_column):
    """ A helper to pickle and run ionization in multiprocessing """
    
    mol = row[mol_column]
    highlights = row.get('highlights', None)
    protonated_smiles = ionizer.ionize(mol=mol)

    results = []
    if len(protonated_smiles) == 1:
        results.append({name_column: row[name_column],
                        mol_column: Chem.MolFromSmiles(protonated_smiles[0]),
                        smiles_column: protonated_smiles[0],
                        'highlights': highlights})
    else:
        for i, smile in enumerate(protonated_smiles):
            results.append({name_column: f"{row[name_column]}_{i+1}",
                            mol_column: Chem.MolFromSmiles(smile),
                            smiles_column: smile,
                            'highlights': highlights})

    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tautomerize a molecule')
    parser.add_argument('-i', '--input', default=None, type=str, help='Input file containing SMILES strings and names')
    parser.add_argument('-s', '--smiles', default=None, type=str, help='SMILES string of the molecule')
    parser.add_argument('-p', '--pH', default=7, type=int, help='pH value to use for ionization')
    parser.add_argument('-r', '--pH_range', default=2, type=int, help='Range of pH values to consider for ionization')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode')
    args = parser.parse_args()
    if args.smiles and args.input:
        raise ValueError("Either SMILES or input file must be provided.")

    ionizer = Ionizer(pH = args.pH, pH_range = args.pH_range, debug = args.debug)
    if args.smiles:
        print(ionizer.ionize(smiles = args.smiles))
    else:
        df = pd.read_csv(args.input, sep =r'\s+', header=None, names=['smiles', 'ids'])
        ionized_df = ionizer.ionize_df(df)
        ionized_df[['smiles', 'ids']].to_csv('protonated_molecules.smi', index=False, header=False, sep =' ')
    #print(ionizer.ionize(smiles = sys.argv[1]))
    #print(ionizer.rules_across_pH)