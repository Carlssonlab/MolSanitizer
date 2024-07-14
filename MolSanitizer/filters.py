from pathlib import Path

import pandas as pd

from rdkit import Chem, RDConfig, RDLogger

from rdkit.Chem import AllChem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

import logging
logger = logging.getLogger('molsani')
RDLogger.DisableLog('rdApp.*') # To disable error messages with kekulizing tautomers from RDKit

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

    df['mol'] = df['mol'].apply(lambda x: stripSMILESsalt(x, remover))
    df['smiles'] = df['mol'].apply(lambda x: Chem.MolToSmiles(x))
    df=df[df['smiles']!=''] #Remove purely salt molecules
    return df


def tautomerize_step1(mol, params):
    """Tautomerize the input molecule using the RDKit TautomerEnumerator class.

    Args:
        mol (RO Mol): The input molecule.
        params (_type_): The RDKit CleanupParameters object.

    Returns:
        RO Mol: Canonical tautomeric form of the input molecule.
    """
    te = rdMolStandardize.TautomerEnumerator(params) 
    canonical_tautomer = te.Canonicalize(mol)
    return canonical_tautomer
    
def tautomers(df: pd.DataFrame, debug = False) -> pd.DataFrame:
    """Tautomers enumeration for the input molecules using the RDKit TautomerEnumerator class and cleaning using an in-house SMARTS reaction list.

    Args:
        df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
        debug (bool, optional): Debug mode. Defaults to False.

    Returns:
        pd.DataFrame: A new DataFrame chunk with tautomerized molecules.
    """
    # Step 1: Tautomerize the input molecules using the RDKit TautomerEnumerator class
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    params.maxTransforms = 10000
    params.maxTautomers = 10000

    df['mol'] = df['mol'].apply(lambda x:  tautomerize_step1(x, params))

    # Step 2: Clean the tautomerized molecules using an in-house SMARTS reaction list
    smartsFile = Path(__file__).parent / 'Data' / 'tautomers.txt'
    reactions = load_reactions(smartsFile)

    # Apply reactions to each SMILES in the DataFrame
    updated_df = []
    for _, row in df.iterrows():
        if debug: logger.info(f"Processing tautomer: {row['ids']}, {Chem.MolToSmiles(row['mol'])}")
        updates = list(recursive_reaction(row['mol'], reactions, set()))
        if (len(updates) == 1): 
            updated_df.append(
                {'smiles': updates[0], 'ids': row['ids'],'mol': Chem.MolFromSmiles(updates[0])})
        else:
            for i, update in enumerate(updates):
                print(update)
                updated_df.append({
                    'smiles': update,
                    'ids': row['ids']+'_'+str(i+1),
                    'mol': Chem.MolFromSmiles(update)
                })
    return pd.DataFrame(updated_df)

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


def stereoisomers(df: pd.DataFrame, max_isomers = 0, debug = False) -> pd.DataFrame:
    """
    Generate stereoisomers for molecules in the 'mol' column and expand the DataFrame.
    
    Ref: https://www.rdkit.org/new_docs/source/rdkit.Chem.EnumerateStereoisomers.html

    Args:
    df (pd.DataFrame): DataFrame with a 'mol' column containing RDKit molecule objects.
    max_isomers (int): Maximum number of stereoisomers to generate for each molecule.
    
    Returns:
    pd.DataFrame: Expanded DataFrame with each stereoisomer as a separate row.
    """
    # TODO: Log of how many stereoisomers for each compounds?
    product_df = []
    for _, row in df.iterrows():
        try:
            isomers = generate_stereoisomers(row['mol'], max_isomers)
        except Exception as e:
            logger.info(f"Error generating stereoisomers for compound {row['ids']}: {Chem.MolToSmiles(row['mol'])}")
            isomers = []
            continue
        if (len(isomers) == 1): 
            product_df.append(
                {'smiles': row['smiles'], 'ids': row['ids'],'mol': row['mol']})
        else:
            for i, isomer in enumerate(isomers):
                product_df.append({
                    'smiles': Chem.MolToSmiles(isomer, isomericSmiles=True),
                    'ids': row['ids']+'_'+str(i+1),
                    'mol': isomer
                })
    return pd.DataFrame(product_df)
   

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

# Function to apply reactions to a molecule
def recursive_reaction(mol, reactions, collection):
    """Recursively conduct the reaction until no more products are generated.

    Args:
        mol (rdkit.Chem.rdchem.Mol object): the reactants for the reaction
        reactions (list):   A list of reactions in the form of SMARTS strings
                            format: [rdkit.Chem.rdchem.Mol object, name]
        collection (set):   A set of unique products (in SMILES) generated from the reaction

    Returns:
        collection (set):   A set of unique products (in SMILES) generated from the reaction
    """
    reactive = False
    for rxn, name in reactions:
        outcomes = rxn.RunReactants((mol,))
        if outcomes:  # Check if there are any outcomes
            reactive = True
            if name == 'amine':
                for outcome in outcomes:
                    product = outcome[0]
                    try:
                        Chem.SanitizeMol(product)
                        recursive_reaction(product, reactions, collection)  # Recurse with the new product
                    except:
                        logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)} to {Chem.MolToSmiles(product)}")
                        pass
            else:
                product = outcomes[0][0]
                try:
                        Chem.SanitizeMol(product)
                        recursive_reaction(product, reactions, collection)  # Recurse with the new product
                except:
                    logger.info(f"Error sanitizing molecule: {Chem.MolToSmiles(mol)} to {Chem.MolToSmiles(product)}")
                    pass
    if not reactive: 
        collection.add(Chem.MolToSmiles(mol))  # Add the initial molecule if it is not reactive
    return collection  # Return the collection of products

def protonation(df: pd.DataFrame, debug = False) -> pd.DataFrame:
    """Protonate the input molecules using a set of predefined reactions.
    The reactions are stored in a file with the following format:
    SMARTS reactants >> products [name]
    Default file: MolSanitizer/Data/ionizations.txt

    Args:
        df (pd.DataFrame): Input DataFrame with 'mol' column containing RDKit molecule objects.
        debug (bool, optional): Debug mode. Defaults to False.

    Returns:
        pd.DataFrame: A new DataFrame chunk with protonated molecules.
    """
    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'ionizations.txt'

    # Load reactions from file
    reactions = load_reactions(smartsFile)

    # Apply reactions to each SMILES in the DataFrame
    charged_df = []
    for _, row in df.iterrows():
        variations = list(recursive_reaction(row['mol'], reactions, set()))
        if (len(variations) == 1): 
            charged_df.append(
                {'smiles': variations[0], 'ids': row['ids'],'mol': Chem.MolFromSmiles(variations[0])})
        else:
            for i, variation in enumerate(variations):
                charged_df.append({
                    'smiles': variation,
                    'ids': row['ids']+'_'+str(i+1),
                    'mol': Chem.MolFromSmiles(variation)
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


    return filtered_df
