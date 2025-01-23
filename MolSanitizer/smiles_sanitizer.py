from pathlib import Path

import pandas as pd

from rdkit import Chem, RDLogger

from rdkit.Chem import AllChem, SaltRemover
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import multiprocessing as mp
from functools import partial
from .filters import Filters
import logging
logger = logging.getLogger('molsani')
RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit

class SmilesSanitizer:
    """
    A class to store the filter options and conduct chemical modifications for molSanitizer. 
    Initialize the class with the desired filter options and apply the filters to the input DataFrame.
    
    Example use:

        sanitizer = SmilesSanitizer(
                    removesalts=True,
                    ha='>10',
                    logp='<5',
                    tautomers=True, 
                    protonation = True, 
                    debug=True)

        santized_df = sanitizer.Sanitize(df)

    """
    tautomer_params = None  # Define as a class variable
    enumerating_reactions = [
    'amine', 'vinylog_acid', 'squaric_acid', 'vinyl_acid', 'vanillin_like',
    'hydrazine', 'p-phenyldiamine', 'biguanides', 'cycloguanil-like'
    ]

    def __init__(self, 
                removesalts=False,
                custom = None,
                unwanted = None,
                pains = False,
                ha = None,
                logp = None,
                hba = None,
                hbd = None,
                mw = None,
                tautomers = False,
                neutralize = True,
                taurdkit = True,
                stereoisomers = False,
                max_stereoisomers = 8,
                protonation = False,
                pH = 7,
                pH_range = 0,
                numcores = 1,
                conformal = False,
                debug = False):

        self.removesalts = removesalts
        self.ha = ha
        self.logp = logp
        self.hba = hba
        self.hbd = hbd
        self.mw = mw
        self.custom = custom
        self.unwanted = unwanted
        self.pains = pains
        
        self.tautomers = tautomers
        self.neutralize = neutralize
        self.taurdkit = taurdkit
        self.stereoisomers = stereoisomers
        self.max_stereoisomers = max_stereoisomers
        self.protonation = protonation
        self.pH = pH
        self.pH_range = pH_range
        self.conformal = conformal    
        self.debug = debug
        self.numcores = numcores
        if SmilesSanitizer.tautomer_params is None: SmilesSanitizer.tautomer_params = self.get_tautomer_params()

    
    @staticmethod
    def get_tautomer_params():
        params = rdMolStandardize.CleanupParameters()
        params.tautomerRemoveSp3Stereo = False
        params.tautomerRemoveBondStereo = False
        params.tautomerRemoveIsotopicHs = False
        params.maxTransforms = 1000
        params.maxTautomers = 1000
        return params
    
    @staticmethod
    def load_reactions(file_path):
        """Load the reactions from a file containing SMARTS strings.

        Args:
            file_path (str): Path to the file containing the reactions in SMARTS strings.

        Returns:
            list: A list containing the reactions in the form of [rdkit.Chem.rdChemReactions object, name].
        """
        reactions = []
        with open(file_path, 'r') as file:
            for line in file:
                smarts = line.strip().split()
                if smarts:
                    reactions.append([AllChem.ReactionFromSmarts(smarts[0]), smarts[1]])
        return reactions

    @staticmethod
    def neutralization(df: pd.DataFrame, debug = False) -> pd.DataFrame:
        '''Neutralize the input molecules using the neutralize_atoms() function. 
        Turn on by default if the user trigger the tautomers or protonation flag.'''
        
        def neutralize_atoms(mol):
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
            return mol
        
        if debug: print('Neutralizing molecules...')
        logger.info('Neutralizing molecules...')
        for i, row in df.iterrows():
            try:
                row['mol'] = neutralize_atoms(row['mol'])
                row['smiles'] = Chem.MolToSmiles(row['mol'])
            except:
                logger.error(f"Error neutralizing molecule: {Chem.MolToSmiles(row['mol'])}")
                pass
        return df
    
    @staticmethod
    def _tautomerize_step1(mol, params):
        """Tautomerize the input molecule using the RDKit TautomerEnumerator class."""
        te = rdMolStandardize.TautomerEnumerator(params)

        try:
            canonical_tautomer = te.Canonicalize(mol)
            if te.ScoreTautomer(mol) == te.ScoreTautomer(canonical_tautomer):
                return mol
        except Exception:
            logger.info(f"Error tautomerizing molecule: {Chem.MolToSmiles(mol)}")
            return mol
        return canonical_tautomer

    @staticmethod
    def _process_molecule_tautomer( row, reactions, taurdkit, debug=False):
        """Process individual molecule: tautomerize and clean with SMARTS reactions."""
        if taurdkit:
            # Step 1: Use RDKit TautomerEnumerator to canonicalize the input molecule
            rdkit_canonical = SmilesSanitizer._tautomerize_step1(row['mol'], SmilesSanitizer.tautomer_params)
            initial_chiral_centers = len(Chem.FindMolChiralCenters(row['mol']))
            rdkit_chiral_centers = len(Chem.FindMolChiralCenters(rdkit_canonical))
            mol = row['mol'] if rdkit_chiral_centers < initial_chiral_centers else rdkit_canonical

            if debug:
                print(f"Using RDKit tautomerizer for {row['ids']}...\n\tTurned to {Chem.MolToSmiles(mol)}")
        else:
            mol = row['mol']

        if debug:
            logger.info(f"Processing tautomer: {row['ids']}, {Chem.MolToSmiles(mol)}")

        # Step 2: Apply corrections
        updates = list(SmilesSanitizer.recursive_reaction(mol, reactions, set(), set(), debug))
        updated_rows = []
        highlights = row.get('highlights', None)

        if len(updates) == 1:
            updated_rows.append({'smiles': updates[0],
                                 'ids': row['ids'],
                                 'mol': Chem.MolFromSmiles(updates[0]),
                                 'highlights': highlights})
        else:
            two_digits = len(updates) >= 10
            for i, update in enumerate(updates):
                updated_rows.append({
                    'smiles': update,
                    'ids': row['ids'] + '_' + (f"{i+1:02}" if two_digits else f"{i+1}"),
                    'mol': Chem.MolFromSmiles(update),
                    'highlights': highlights
                })

        return updated_rows

    @staticmethod
    def tautomerization(df: pd.DataFrame, taurdkit=True, num_cores=4, debug=False) -> pd.DataFrame:
        """
        High-level method for tautomers enumeration. Encapsulates all tautomerization logic.
        """
        smartsFile = Path(__file__).parent / 'Data' / 'tautomers.txt'
        reactions = SmilesSanitizer.load_reactions(smartsFile)

        # Parallel processing setup
        process_func = partial(SmilesSanitizer._process_molecule_tautomer, reactions=reactions, taurdkit=taurdkit, debug=debug)
        results = []

        with mp.Pool(processes=num_cores) as pool:
            # Submit all tasks to the pool in parallel
            async_results = [(row, pool.apply_async(process_func, (row,))) for _, row in df.iterrows()]

            # Collect results as they complete
            for row_data, async_result in async_results:
                try:
                    chunk_result = async_result.get(timeout=60)
                    results.extend(chunk_result)
                except (mp.TimeoutError, Exception) as e:
                    # Handle exceptions during processing
                    if isinstance(e, mp.TimeoutError):
                        logger.warning(f"Timeout occurred for compound {row_data['ids']}. Using original molecule.")
                    else:
                        logger.error(f"Error processing compound {row_data['ids']}: {str(e)}. Using original molecule.")

                    highlights = row_data.get('highlights', None)
                    results.append({
                        'smiles': row_data['smiles'],
                        'ids': row_data['ids'],
                        'mol': row_data['mol'],
                        'highlights': highlights
                    })

        return pd.DataFrame(results)

    
    
    @staticmethod
    def load_protonation_rules(file_path):
        '''
        Load the protonation rules from a file containing SMARTS strings.'''

        rules = pd.read_csv(file_path, sep=r"\s+")
        rules['Mol'] = rules['REACTION'].apply(lambda x: AllChem.ReactionFromSmarts(x))
        return rules

    @staticmethod
    def extract_necessary_rules(rule_library, pH):
        """
        Processes ionization rules from a file and filters them based on pH.
        
        Parameters:
            file_path (str): Path to the ionization rules file.
            pH (float): The pH value to filter the rules.
            
        Returns:
            pd.DataFrame: Filtered rules based on the pH.
        """


        # Filter based on pH for ACID and BASE rules
        filtered_rules = rule_library[
            ((rule_library['TYPE'] == 'ACID') & (rule_library['pKa'] < pH)) |
            ((rule_library['TYPE'] == 'BASE') & (rule_library['pKa'] > pH))
        ]

        reaction_list = []

        for i, row in filtered_rules.iterrows():
            reaction_list.append((row['Mol'], row['FUNCTIONAL_GROUP']))
        logger.info(f'Parsed {len(reaction_list)} rules for pH {pH}')    
        return reaction_list

    @staticmethod
    def recursive_reaction(mol, reactions, collection, visited=None, debug=False):
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
                if name in SmilesSanitizer.enumerating_reactions:
                    for outcome in outcomes:
                        product = outcome[0]
                        try:
                            error = Chem.SanitizeMol(product, catchErrors=True)
                            if error == 0:
                                last_successful_smiles = Chem.MolToSmiles(product)
                                SmilesSanitizer.recursive_reaction(product, reactions, collection, visited, debug)
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
                            SmilesSanitizer.recursive_reaction(product, reactions, collection, visited, debug)
                        else:
                            if debug: print(f"Sanitization error for molecule: {Chem.MolToSmiles(product)}")
                            reactive = False
                            continue
                    except Exception as e:
                        logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(product)}. Exception: {e}")
        
        if not reactive: 
            collection.add(last_successful_smiles)  # Add the last valid molecule if it is not reactive
        return collection  # Return the collection of products

    def ionization(df: pd.DataFrame, pH: int = 7, pH_range: int = 0, debug = False) -> pd.DataFrame:
        """Protonate the input molecules using a set of predefined reactions.
        The reactions are stored in a file with the following format:
        SMARTS reactants >> products [name]
        Default file: MolSanitizer/Data/ionizations_v2.txt

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            pH (int, optional): The pH value to use for protonation. Defaults to 7.
            pH_range (int, optional): The range of pH values to consider. Defaults to 0.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with protonated molecules.
        """
        # Get the absolute path to the template SMARTS file using pathlib
        smartsFile = Path(__file__).parent / 'Data' / 'ionizations_v2.txt'

        # Load reactions from file
        rule_library = SmilesSanitizer.load_protonation_rules(smartsFile)
        pH_values = (pH,) if pH_range == 0\
                    else (max(0, pH - pH_range), pH, min(14, pH + pH_range))

        all_rules = {}
        for pH in pH_values:
            rules = SmilesSanitizer.extract_necessary_rules(rule_library, pH)
            all_rules[pH] = rules
        # Apply reactions to each SMILES in the DataFrame
        charged_df = []
        for _, row in df.iterrows():
            if debug: print(f"Protonating: {row['ids']}, {Chem.MolToSmiles(row['mol'])}")
            variation_sets = set() 
            highlights = row.get('highlights', None)   
            for pH in pH_values:
                variations = list(SmilesSanitizer.recursive_reaction(row['mol'], all_rules[pH], set()))
                variation_sets.update(variations)
            if (len(variation_sets) == 1): 
                charged_df.append(
                    {'smiles': list(variation_sets)[0], 
                    'ids': row['ids'],
                    'mol': Chem.MolFromSmiles(list(variation_sets)[0]),
                    'highlights': highlights})
            else:
                for i, variation in enumerate(variation_sets):
                    charged_df.append({
                        'smiles': variation,
                        'ids': row['ids']+'_'+str(i+1),
                        'mol': Chem.MolFromSmiles(variation),
                        'highlights': highlights
                    })
        return pd.DataFrame(charged_df)
    
    @staticmethod
    def _generate_stereoisomers(mol, max_isomers):
        """
        Generate stereoisomers for a given molecule.
        
        Args:
        mol (rdkit.Chem.Mol): RDKit molecule object.
        max_isomers (int): Maximum number of stereoisomers to generate.
        
        Returns:
        list: List of stereoisomer molecules.
        """
        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True, maxIsomers=max_isomers)
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        return isomers

    @staticmethod
    def _process_molecule_stereoisomer(row_data, max_isomers=8):
        """
        A wrapper to generate stereoisomers for a given molecule using multiprocessing.
        Args:
        row_data (pd.Series): A row from the input DataFrame containing ['smiles', 'ids', 'mol'] and optionally 'highlights'.
        max_isomers (int): Maximum number of stereoisomers to generate.

        Returns:
        results (list): List of expanded stereoisomers dictionaries with 'smiles', 'ids', 'mol', and optionally 'highlights'.
        """
        mol = Chem.MolFromSmiles(row_data['smiles'])
        try:
            isomers = SmilesSanitizer._generate_stereoisomers(mol, max_isomers=max_isomers)
        except Exception as e:
            logger.error(f"Error generating stereoisomers for compound {row_data['ids']}: {row_data['smiles']}")
            isomers = [mol]
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        unassigned = [idx for idx, tag in centers if tag == '?']
        num_possible_isomers = 2 ** len(unassigned)
        result = []
        # Get highlights if available
        highlights = row_data.get('highlights', None)
        if len(isomers) == 1:
            result.append({'smiles': Chem.MolToSmiles(isomers[0]),
                            'ids': row_data['ids'],
                            'mol': isomers[0],
                            'highlights': highlights})
        else:
            if max_isomers > 0 and num_possible_isomers > max_isomers:
                logger.warning(f"{row_data['ids']}: Not all the stereoisomers are written out (capped at {max_isomers}/{num_possible_isomers}).")
                isomers = isomers[:max_isomers]
            two_digits = len(isomers) >= 10
            for i, isomer in enumerate(isomers):
                result.append({
                    'smiles': Chem.MolToSmiles(isomer, isomericSmiles=True),
                    'ids': row_data['ids'] + '.' + (f"{i+1:02}" if two_digits else f"{i+1}"),
                    'mol': isomer,
                    'highlights': highlights
                })
        
        return result

    @staticmethod
    def enum_stereoisomers(df: pd.DataFrame, max_isomers=8, numcores=4, debug=False) -> pd.DataFrame:
        """
        Generate stereoisomers for molecules in the 'smiles' column and expand the DataFrame using multiprocessing.
        
        Args:
        df (pd.DataFrame): DataFrame with 'smiles' and 'ids' columns.
        max_isomers (int): Maximum number of stereoisomers to generate for each molecule.
        numcores (int): Number of processes to use. Default is 4.
        debug (bool): Enable debug messages.
        
        Returns:
        pd.DataFrame: Expanded DataFrame with each stereoisomer as a separate row.
        """
        # Partial function to fix max_isomers as an argument
        process_func = partial(SmilesSanitizer._process_molecule_stereoisomer, max_isomers=max_isomers)
        results = []

        with mp.Pool(processes=numcores) as pool:
            # Submit all tasks to the pool in parallel, keeping track of the rows for error handling
            async_results = [(row, pool.apply_async(process_func, (row,))) for _, row in df.iterrows()]

            # Collect results as they complete
            for row_data, async_result in async_results:
                try:
                    chunk_result = async_result.get(timeout=60)
                    results.extend(chunk_result)
                except (mp.TimeoutError, Exception) as e:
                    # Differentiate the logging based on the type of exception if desired
                    if isinstance(e, mp.TimeoutError):
                        logger.warning(f"Timeout occurred for compound {row_data['ids']}. Using original molecule.")
                    else:
                        logger.error(f"Error processing compound {row_data['ids']}: {str(e)}. Using original molecule.")
                    
                    highlights = row_data.get('highlights', None)
                    results.append({
                        'smiles': row_data['smiles'],
                        'ids': row_data['ids'],
                        'mol': row_data['mol'],
                        'highlights': highlights
                    })

        return pd.DataFrame(results)

    def Sanitize(self, df: pd.DataFrame, rejected_file) -> pd.DataFrame:
        """
        Sanitize the input DataFrame using the specified filters and rule-based chemical modifications.
        
        Args:
            df (pd.DataFrame): Input DataFrame with ['smiles', 'ids'] columns. 
            
        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules that passed the filters.
        """
        df['mol'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x))
        df = Filters.remove_invalid_SMILES(df)
        if self.removesalts: df = Filters.saltstripping(df, debug=self.debug)
        if self.ha is not None: df = Filters.filter_by_ha(df, self.ha, rejectedFile=rejected_file, debug=self.debug)
        if self.logp is not None: df = Filters.filter_by_logp(df, self.logp, rejectedFile=rejected_file, debug=self.debug)
        if self.hba is not None: df = Filters.filter_by_hba(df, self.hba, rejectedFile=rejected_file, debug=self.debug)
        if self.hbd is not None: df = Filters.filter_by_hbd(df, self.hbd, rejectedFile=rejected_file, debug=self.debug)
        if self.mw is not None: df = Filters.filter_by_mw(df, self.mw, rejectedFile=rejected_file, debug=self.debug)
        if self.conformal: 
            df = Filters.standarizeFilters(df)
            return df
        
        if self.tautomers or self.protonation:
            if self.neutralize: df = SmilesSanitizer.neutralization(df)
        if self.tautomers: df = SmilesSanitizer.tautomerization(df, taurdkit=self.taurdkit, num_cores=self.numcores, debug=self.debug)
        if self.pains: df = Filters.painsFilter(df, rejectedFile=rejected_file, debug=self.debug)
        if self.unwanted is not None: df = Filters.unwantedFilter(df, rejectedFile=rejected_file, unwanted_option=self.unwanted, debug=self.debug)
        if self.custom is not None: df = Filters.customFilter(df, rejectedFile=rejected_file, smartsFile = self.custom,  debug=self.debug)
        if self.protonation: df = SmilesSanitizer.ionization(df, pH=self.pH, pH_range=self.pH_range, debug=self.debug)
        if self.stereoisomers: df = SmilesSanitizer.enum_stereoisomers(df, max_isomers=self.max_stereoisomers, debug=self.debug)
        return df
    

def loadSMARTSdata(smartsFile: str, unwanted_option=None) -> pd.DataFrame:
        """Load SMARTS patterns from a file and convert them to RDKit molecule objects.

        Args:
            smartsFile (str): Path to the file containing the SMARTS patterns.
            unwanted_option (list): The mode input by the user. Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the SMARTS patterns and their corresponding RDKit molecule objects.
        """
        if unwanted_option is not None: 
            # Using default substructure file 
            smarts_df = pd.read_csv(smartsFile, sep='\s+', header=0, names=['smarts','label', 'reason', 'mode', 'ref'])
            smarts_df = smarts_df[smarts_df["mode"].isin(unwanted_option)]
        else:
            # Using customized substructure file
            has_header = Filters.check_header(smartsFile)
            if has_header:
                logger.info(f'Found header in {smartsFile}')
                smarts_df = pd.read_csv(smartsFile, sep='\s+', header=0, usecols=[0,1], names=['smarts','label'])
            else:
                smarts_df = pd.read_csv(smartsFile, sep='\s+', header=None, usecols=[0,1], names=['smarts','label'])

        smarts_df['mol'] = smarts_df['smarts'].apply(lambda x: Chem.MolFromSmarts(x)) #do we need mergeHs here?
        
        return smarts_df

def load_reactions(file_path):
        """Load the reactions from a file containing SMARTS strings.

        Args:
            file_path (str): Path to the file containing the reactions in SMARTS strings.

        Returns:
            list: A list containing the reactions in the form of [rdkit.Chem.rdChemReactions object, name].
        """
        reactions = []
        with open(file_path, 'r') as file:
            for line in file:
                smarts = line.strip().split()
                if smarts:
                    reactions.append([AllChem.ReactionFromSmarts(smarts[0]), smarts[1]])
        return reactions
