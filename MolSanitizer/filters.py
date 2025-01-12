from pathlib import Path

import pandas as pd

from rdkit import Chem, RDLogger

from rdkit.Chem import AllChem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import multiprocessing as mp
from functools import partial
import logging
logger = logging.getLogger('molsani')
RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit

def get_tautomer_params():
    """Get standardized RDKit tautomer parameters."""
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    params.maxTransforms = 1000
    params.maxTautomers = 1000
    
    return params

# Initialize parameters once at module level
TAUTOMER_PARAMS = get_tautomer_params()
#uncharger = rdMolStandardize.Uncharger()
enumerating_reactions = ['amine', 'vinylog_acid', 'squaric_acid', 'vinyl_acid', 'vanillin_like',
                         'hydrazine', 'p-phenyldiamine', 'biguanides', 'cycloguanil-like']
        
def remove_invalid_SMILES(df:pd.DataFrame) -> pd.DataFrame:
    """Remove rows with invalid SMILES from the input DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame with 'smiles' column containing SMILES strings.

    Returns:
        pd.DataFrame: A new DataFrame chunk with valid SMILES strings.
    """
    # Log rows where 'mol' is None before dropping
    invalid_rows = df[df['mol'].isna()]
    for index, row in invalid_rows.iterrows():
        logger.warning(f"INVALID SMILES: Removing row at index {index}: {row.to_dict()}")

    # Remove rows where 'mol' is None
    df_cleaned = df.dropna(subset=['mol'])

    return df_cleaned


def stripSMILESsalt(mol, molRemover, debug = False):
    """Strip salts from the input molecule using the RDKit SaltRemover class.

    Args:
        mol (rdkit mol object): The input molecule.
        molRemover: The RDKit SaltRemover object.
        debug (bool, optional): Debug mode. Defaults to False.

    Returns:
        rdkit mol object: The stripped molecule.
    """
    res, deleted = molRemover.StripMolWithDeleted(mol)

    if debug and len(deleted) > 0: logger.info(f"Stripped salt {Chem.MolToSmiles(mol)}:  Salts:{' '.join([Chem.MolToSmiles(m) for m in deleted])}")
    
    if len(Chem.rdmolops.GetMolFrags(res)) > 1:
        # If still contains more than one fragment, retains the largest one
        rdMolStandardize.FragmentParentInPlace(res)
    return res

def removesalts(df: pd.DataFrame, debug = False) -> pd.DataFrame:
    """Remove salts from the input molecules using the RDKit SaltRemover class.

    Args:
        df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
        debug (bool, optional): Debug mode. Defaults to False.

    Returns:
        pd.DataFrame: A new DataFrame chunk with salt-stripped molecules.
    """
    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'salt_stripping.txt'
    if debug: logger.info(f'Parsing salts SMARTS file: {smartsFile.resolve()}')

    remover = SaltRemover.SaltRemover(defnFilename=smartsFile)

    # Remove salts in the list
    df['mol'] = df['mol'].apply(lambda x: stripSMILESsalt(x, remover))

    df['smiles'] = df['mol'].apply(lambda x: Chem.MolToSmiles(x))
    df=df[df['smiles']!=''] #Remove purely salt molecules
    return df

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

def neutralize(df: pd.DataFrame, debug = False) -> pd.DataFrame:
    '''Neutralize the input molecules using the neutralize_atoms() function. 
    Turn on by default if the user trigger the tautomers or protonation flag.'''

    if debug: print('Neutralizing molecules...')
    logger.info('Neutralizing molecules...')
    # for atom in df.iloc['mol'][0].GetAtoms():
    #     print(atom.GetFormalCharge())
    for i, row in df.iterrows():
        try:
            row['mol'] = neutralize_atoms(row['mol'])
            row['smiles'] = Chem.MolToSmiles(row['mol'])
        except:
            logger.error(f"Error neutralizing molecule: {Chem.MolToSmiles(row['mol'])}")
            pass
    return df

def tautomerize_step1(mol, params=TAUTOMER_PARAMS):
    """Tautomerize the input molecule using the RDKit TautomerEnumerator class.

    Args:
        mol (RO Mol): The input molecule.
        params (_type_): The RDKit CleanupParameters object.

    Returns:
        RO Mol: Canonical tautomeric form of the input molecule.
    """
    te = rdMolStandardize.TautomerEnumerator(params) 

    try:
        canonical_tautomer = te.Canonicalize(mol)
    except: # If tautomerization fails, return the input molecule
        logger.info(f"Error tautomerizing molecule: {Chem.MolToSmiles(mol)}")
        canonical_tautomer = mol
    return canonical_tautomer

def process_molecule_tautomer(row, reactions, taurdkit, debug=False):
    """Process individual molecule: tautomerize and clean with SMARTS reactions."""
    if taurdkit:
    # Step 1: Use RDkit TautomerEnumerator to canonicalize the input molecule
        rdkit_canonical = tautomerize_step1(row['mol'])
        initial_chiral_centers  = len(Chem.FindMolChiralCenters(row['mol']))
        rdkit_chiral_centers = len(Chem.FindMolChiralCenters(rdkit_canonical))
        # If canonical tauatomer has less chiral centers than the initial molecule, keep the initial molecule
        if rdkit_chiral_centers < initial_chiral_centers: mol = row['mol']
        else: mol = rdkit_canonical
    else: mol = row['mol']
    
    if debug:
        logger.info(f"Processing tautomer: {row['ids']}, {Chem.MolToSmiles(mol)}")
    
    # Step 2: Apply corrections
    updates = list(recursive_reaction(mol, reactions, set()))
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

def tautomers(df: pd.DataFrame, taurdkit=True, num_cores=4, debug=False) -> pd.DataFrame:
    """Parallelized tautomers enumeration using multiprocessing."""

    smartsFile = Path(__file__).parent / 'Data' / 'tautomers.txt'
    reactions = load_reactions(smartsFile)

    # Parallel processing setup
    process_func = partial(process_molecule_tautomer, reactions=reactions, taurdkit=taurdkit, debug=debug)
    results = []

    with mp.Pool(processes=num_cores) as pool:
        # Submit all tasks to the pool in parallel, keeping track of rows for error handling
        async_results = [(row, pool.apply_async(process_func, (row,))) for _, row in df.iterrows()]

        # Collect results as they complete
        for row_data, async_result in async_results:
            try:
                chunk_result = async_result.get(timeout=60)
                results.extend(chunk_result)  # Append all processed rows
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

def detectPAINS(mol, catalog):
    """Detect PAINS functional groups in a molecule.

    Args:
        mol (Rdkit Mol object): The molecule to be checked.
        catalog (FilterCatalog): The PAINS catalog.

    Returns:
        str: Return 'OK' if no PAINS functional groups are detected, otherwise return the description of the first PAINS functional group detected.
    """
    entry = catalog.GetFirstMatch(mol)  # Get the first matching PAINS
    if entry is not None:
        return ('PAINS violation: ' + str(entry.GetDescription().capitalize()))  # Indicate a PAINS violation
    else:
        return 'OK'

def pains(df: pd.DataFrame, rejectedFile, debug = False) -> pd.DataFrame:
    """Filter out PAINS functional groups using the RDKit PAINS catalog.

    Args:
        df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
        rejectedFile (_type_): Path to the file to save rejected molecules.
        debug (bool, optional): Debug mode. Defaults to False.

    Returns:
        pd.DataFrame: A new DataFrame chunk with molecules that passed the PAINS filter.
    """
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    df['reason'] = df['mol'].apply(lambda x: detectPAINS(x, catalog))
    rejected_df=df[df['reason']!='OK']
    rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','reason'], sep = ' ', header=False)
    df=df[df['reason']=='OK']
    return df

def check_header(file_path):
    """
    Check if the file has a header containing 'SMARTS' in the first column.
    Args:
        file_path (Path): Path to the CSV file.
        
    Returns:
        bool: True if the file has a header, False otherwise.
    """
    with open(file_path, 'r') as file:
        first_line = file.readline().strip().upper()
        return 'SMARTS' in first_line
    
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
        has_header = check_header(smartsFile)
        if has_header:
            logger.info(f'Found header in {smartsFile}')
            smarts_df = pd.read_csv(smartsFile, sep='\s+', header=0, usecols=[0,1], names=['smarts','label'])
        else:
            smarts_df = pd.read_csv(smartsFile, sep='\s+', header=None, usecols=[0,1], names=['smarts','label'])

    smarts_df['mol'] = smarts_df['smarts'].apply(lambda x: Chem.MolFromSmarts(x)) #do we need mergeHs here?
    
    return smarts_df

def filterbysmarts(mol, smarts_df: pd.DataFrame) -> str:
    for _, substructure in smarts_df.iterrows():
        if mol.HasSubstructMatch(substructure.mol):
            return substructure.label
    return 'OK'

def unwanted(df: pd.DataFrame, rejectedFile, unwanted_option, debug = False) -> pd.DataFrame:
    """Filter out unwanted substructures using a default list of SMARTS patterns.

        Args:
            df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
            rejectedFile (str): Path to the file to save rejected molecules.
            unwanted_option (list): The mode input by the user.
            debug (bool, optional): Debug mode. Defaults to False.

        Returns:
            pd.DataFrame: A new DataFrame chunk with molecules that passed the filter.
    """
    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'filter_out.csv'

    # Load smarts to clean  from file
    unwanted_df = loadSMARTSdata(smartsFile.resolve(), unwanted_option)

    # Apply reactions to each SMILES in the DataFrame
    df['reason'] = df['mol'].apply(lambda x: filterbysmarts(x, unwanted_df))
    rejected_df=df[df['reason']!='OK']
    rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','reason'], sep = ' ', header=False)
    df=df[df['reason']=='OK']
    return df

def custom(df: pd.DataFrame, rejectedFile, smartsFile, debug = False) -> pd.DataFrame:
    """Filter out unwanted substructures using a customized list of SMARTS patterns.

    Args:
        df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
        rejectedFile (_type_): Path to the file to save rejected molecules.
        smartsFile (_type_): Path to the file containing the customized list of SMARTS patterns.
        debug (bool, optional): Debug mode. Defaults to False.

    Returns:
        pd.DataFrame: A new DataFrame chunk with molecules that passed the filter.
    """
    # Load smarts to clean  from file
    unwanted_df = loadSMARTSdata(smartsFile)

    # Apply reactions to each SMILES in the DataFrame
    df['reason'] = df['mol'].apply(lambda x: filterbysmarts(x, unwanted_df))
    rejected_df=df[df['reason']!='OK']
    rejected_df.to_csv(rejectedFile, index=False, mode='a', columns=['smiles','ids','reason'], sep = ' ', header=False)
    df=df[df['reason']=='OK']
    return df

def generate_stereoisomers(mol, max_isomers=0):
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


def process_molecule_stereoisomer(row_data, max_isomers=8, is_synthon=False):
    mol = Chem.MolFromSmiles(row_data['smiles'])
    try:
        isomers = generate_stereoisomers(mol, max_isomers=max_isomers)
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

def stereoisomers(df: pd.DataFrame, max_isomers=32, numcores=4, debug=False) -> pd.DataFrame:
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
    process_func = partial(process_molecule_stereoisomer, max_isomers=max_isomers)
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
   

# Function to read SMARTS reactions from a file
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

def load_protonation_rules(file_path):
    '''
    Load the protonation rules from a file containing SMARTS strings.'''

    rules = pd.read_csv(file_path, sep=r"\s+")
    rules['Mol'] = rules['REACTION'].apply(lambda x: AllChem.ReactionFromSmarts(x))
    return rules

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

# Function to apply reactions to a molecule
def recursive_reaction(mol, reactions, collection, visited=None):
    """
    Recursively apply reactions to the molecule until no new products are generated.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The reactant molecule.
        reactions (list): A list of reactions in the form [(rdkit.Chem.rdChem.Reaction object, name), ...].
        collection (set): A set of unique products (in SMILES) generated from the reaction.
        visited (set): A set of SMILES strings for molecules already processed to avoid redundancy.

    Returns:
        set: A set of unique products (in SMILES) generated from the reaction.
    """
    if visited is None:
        visited = set()
    
    mol_smiles = Chem.MolToSmiles(mol)
    
    # Check if the molecule has already been visited
    if mol_smiles in visited:
        return collection  # Skip redundant processing

    visited.add(mol_smiles)  # Mark the molecule as visited

    reactive = False
    for rxn, name in reactions:
        outcomes = rxn.RunReactants((mol,))
        if outcomes:  # Check if there are any outcomes            
            reactive = True
            if name in enumerating_reactions:
                for outcome in outcomes:
                    product = outcome[0]
                    try:
                        Chem.SanitizeMol(product)
                        recursive_reaction(product, reactions, collection, visited)  # Recurse with the new product
                    except:
                        print(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)} to {Chem.MolToSmiles(product)}")
                        logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)} to {Chem.MolToSmiles(product)}")
                        pass
            else:
                product = outcomes[0][0]
                try:
                        Chem.SanitizeMol(product)
                        recursive_reaction(product, reactions, collection, visited)  # Recurse with the new product
                except:
                    print(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)} to {Chem.MolToSmiles(product)}")
                    logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)} to {Chem.MolToSmiles(product)}")
                    pass
    if not reactive: 
        collection.add(mol_smiles)  # Add the initial molecule if it is not reactive
    return collection  # Return the collection of products

def protonation(df: pd.DataFrame, pH: int = 7, pH_range: int = 0, debug = False) -> pd.DataFrame:
    """Protonate the input molecules using a set of predefined reactions.
    The reactions are stored in a file with the following format:
    SMARTS reactants >> products [name]
    Default file: MolSanitizer/Data/ionizations.txt

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
    rule_library = load_protonation_rules(smartsFile)
    pH_values = (pH,) if pH_range == 0\
                else (max(0, pH - pH_range), pH, pH + pH_range)

    all_rules = {}
    for pH in pH_values:
        rules = extract_necessary_rules(rule_library, pH)
        all_rules[pH] = rules
    # Apply reactions to each SMILES in the DataFrame
    charged_df = []
    for _, row in df.iterrows():
        variation_sets = set() 
        highlights = row.get('highlights', None)   
        for pH in pH_values:
            variations = list(recursive_reaction(row['mol'], all_rules[pH], set()))
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



def applyStandarizeFilters(mol, params):

    taut_uncharged_parent_clean_mol = None

    try:

        clean_mol = rdMolStandardize.Cleanup(mol, params) 

        # if many fragments, get the "parent"
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)

        # try to neutralize molecule
        uncharger = rdMolStandardize.Uncharger()
        uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
        
        # tautomer enumerator
        te = rdMolStandardize.TautomerEnumerator(params) 
        taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)

    except:

        logger.info(f'Molecule NOT processed: {Chem.MolToSmiles(mol)}')

    return taut_uncharged_parent_clean_mol


def standarizeFilters(df: pd.DataFrame) -> pd.DataFrame:

    # Flavio filters to prepara SMILES databases

    # follows the steps in
    # https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb
    # https://www.youtube.com/watch?v=eWTApNX8dJQ
    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule

    #Sc, Y, In, Sn, W, Ac are not included for now (need to check how they can hit other structures).
    organometallics={'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Ga', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',\
    'Cd', 'La', 'Hf ', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Ce', 'Pr', 'Nd',\
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf ', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Ge', 'Sb'}


    # Filter out rows where 'smiles' contains any organometallics
    filtered_df = df[~df['smiles'].apply(lambda x: any(om in x for om in organometallics))]

    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False

    filtered_df['mol'] = filtered_df['mol'].apply(lambda x:  applyStandarizeFilters(x, params))
    filtered_df['smiles'] = filtered_df['mol'].apply(lambda x: Chem.MolToSmiles(x))

    return filtered_df
