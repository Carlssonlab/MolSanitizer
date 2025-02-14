import argparse
import multiprocessing as mp
from functools import partial
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import pandas as pd
import logging
from pathlib import Path
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit
logger = logging.getLogger('molsani')

TAUTOMER_RULES_PATH = Path(__file__).parent / 'Data' / 'tautomers_v2.txt'

# These can be reused along multiprocessing and needs to be outside the class to resolve pickling problem.
TAUTOMER_PARAMS = rdMolStandardize.CleanupParameters()
TAUTOMER_PARAMS.tautomerRemoveSp3Stereo = False
TAUTOMER_PARAMS.tautomerRemoveBondStereo = False
TAUTOMER_PARAMS.tautomerRemoveIsotopicHs = False
TAUTOMER_PARAMS.maxTransforms = 1000
TAUTOMER_PARAMS.maxTautomers = 1000
te = rdMolStandardize.TautomerEnumerator(TAUTOMER_PARAMS) 


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
    
    >>> from MolSanitizer.tautomerizer import Tautomerizer\n
    >>> tautomerizer = Tautomerizer(numcores= 4, neutralize= False)\n
    >>> tautomers = tautomerizer.tautomerize(smiles='c1ccccc1O')\n
    >>> tautomers = tautomerizer.tautomerize(mol = RDKit Mol object)\n
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
        self.reactions = self.load_reactions(smartsFile)
        self.standardizing_reactions = [r for r in self.reactions if not r[1]]
        self.enumerating_reactions = [r for r in self.reactions if r[1]]
        self.debug = debug

    def load_reactions(self, file_path: str):
        """Load the reactions from a file containing SMARTS strings.

        Args:
            file_path (str): Path to the file containing the reactions in SMARTS strings.

        Returns:
            list: A list containing the reactions in the form of [rdkit.Chem.rdChemReactions object, name].
        """
        reactions = []
        with open(file_path, 'r') as file:
            next(file)  # Skip first line
            for line in file:
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
    
    @staticmethod
    def neutralize_mol(mol: Chem.Mol):
        """Neutralize the input molecule by balancing the charges on atoms.
        Adapted from RDKit Cookbook: https://rdkit.org/docs/Cookbook.html"""

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
        
    def tautomer_canonicalize_rdkit(self, mol: Chem.Mol):
        """Tautomerize the input molecule using the RDKit TautomerEnumerator class.
        This is the first step to make the input result in a canonical tautomer for later fixes."""


        try:
            canonical_tautomer = te.Canonicalize(mol)
            # If the canonical tautomer is the same SCORE as the input,
            # we believe more in the input than the output.
            # Return the input molecule
            if te.ScoreTautomer(mol) == te.ScoreTautomer(canonical_tautomer):
                return mol
        except Exception:
            logger.info(f"Error tautomerizing molecule: {Chem.MolToSmiles(mol)}")
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
            while mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                if self.debug: print(f"Applying {name} to {Chem.MolToSmiles(mol)}")
                new_mol = rxn.RunReactants((mol,))[0][0]
                error = Chem.SanitizeMol(new_mol, catchErrors=True)
                if error == 0:
                    mol = new_mol
                else:
                    logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)}")
                    break

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
                    mol: Chem.Mol = None) -> list:
        """
        Tautomerize the input molecule.
        Args: (Either SMILES or RDKit molecule object must be provided)
            smiles (str): SMILES string of the molecule.
            mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.
        
        Returns:
            Returns a list SMILES strings of the tautomers.
        """
        
        if (mol is None) and (smiles is None):
            raise ValueError("Either SMILES or RDKit molecule object must be provided.")
        
        if mol is None and smiles:
            mol = Chem.MolFromSmiles(smiles)
        
        if self.neutralize:
            mol = Tautomerizer.neutralize_mol(mol)
        if self.taurdkit:
            # Step 1: Use RDKit TautomerEnumerator to canonicalize the input molecule
            rdkit_canonical = self.tautomer_canonicalize_rdkit(mol)
            initial_chiral_centers = len(Chem.FindMolChiralCenters(mol))
            rdkit_chiral_centers = len(Chem.FindMolChiralCenters(rdkit_canonical))
            # Canonicalization loses the stereochemistry,
            # so we keep the original molecule if it has more chiral centers
            mol = mol if rdkit_chiral_centers < initial_chiral_centers else rdkit_canonical

            if self.debug:
                print(f"Using RDKit tautomerizer for {smiles}...\n\tTurned to {Chem.MolToSmiles(mol)}")
        
        # Step 2: Standardize the molecule
        standardized_mol = self.standardize(mol)
        if self.debug:
            print(f"Standardized to {Chem.MolToSmiles(standardized_mol)}")
        unique_tautomers = self.enumerate(standardized_mol)
        return unique_tautomers
    
    def tautomerize_df(self, df: pd.DataFrame,
                    smiles_column: str = 'smiles',
                    mol_column: str = 'mol',
                    name_column: str = 'ids') -> pd.DataFrame:
        """
        Tautomerize a dataframe of molecules using multiprocessing.

        Args:
            df (pd.DataFrame): The input dataframe containing the molecules.
            smiles_column (str): The name of the column containing SMILES strings.
            mol_column (str): The name of the column containing RDKit molecule objects.
            name_column (str): The name of the column containing molecule names.

        Returns:
            pd.DataFrame: The expanded dataframe with tautomers.
        """

        # Ensure mol_column is populated
        if mol_column not in df.columns:
            df[mol_column] = df[smiles_column].apply(Chem.MolFromSmiles)

        # Create a partial function to pass self and other parameters
        process_func = partial(_process_tautomer, tautomerizer=self,
                            smiles_column=smiles_column,
                            mol_column=mol_column,
                            name_column=name_column)

        results = []
        
        # Use multiprocessing Pool
        num_cores = self.numcores if self.numcores > 1 else 1  # Ensure at least 1 core
        with mp.Pool(processes=num_cores) as pool:
            # Submit tasks asynchronously
            async_results = [pool.apply_async(process_func, (row,)) for _, row in df.iterrows()]

            # Collect results
            for async_result in async_results:
                try:
                    results.extend(async_result.get(timeout=60))  # Get result with a timeout
                except mp.TimeoutError:
                    logger.warning("Timeout occurred while tautomerizing a molecule. Skipping.")
                except Exception as e:
                    logger.error(f"Error processing a molecule: {str(e)}")

        return pd.DataFrame(results)


def _process_tautomer(row, tautomerizer, smiles_column, mol_column, name_column):
        """ A helper to pickle and run tautomer enumeration in multiprocessing """
        mol = row[mol_column]
        highlights = row.get('highlights', None)
        tautomers_smiles = tautomerizer.tautomerize(mol=mol)

        results = []
        if len(tautomers_smiles) == 1:
            results.append({name_column: row[name_column],
                            mol_column: Chem.MolFromSmiles(tautomers_smiles[0]),
                            smiles_column: tautomers_smiles[0],
                            'highlights': highlights})
        else:
            two_digits = len(tautomers_smiles) >= 10
            for i, tautomer in enumerate(tautomers_smiles):
                results.append({name_column: row[name_column] + '_' + (f"{i+1:02}" if two_digits else f"{i+1}"),
                                mol_column: Chem.MolFromSmiles(tautomer),
                                smiles_column: tautomer,
                                'highlights': highlights})
        return results


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Tautomerize a molecule')
    parser.add_argument('-i', '--input', default=None, type=str, help='Input file containing SMILES strings and names')
    parser.add_argument('-s', '--smiles', default=None, type=str, help='SMILES string of the molecule')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode')
    args = parser.parse_args()
    if args.smiles and args.input:
        raise ValueError("Either SMILES or input file must be provided.")
    tautomerizer = Tautomerizer(debug=args.debug)
    if args.smiles:
        tautomers = tautomerizer.tautomerize(smiles=args.smiles)
        for tautomer in tautomers: print(tautomer)
    elif args.input:
        df = pd.read_csv(args.input, names = ['smiles', 'ids'], sep = r'\s+', header=None)
        tautomers = tautomerizer.tautomerize_df(df)
        tautomers[['smiles', 'ids']].to_csv('tautomers_output.smi', header=False, sep = ' ', index=False)
