from pathlib import Path
import csv

import pandas as pd

from rdkit import Chem, RDConfig

from rdkit.Chem import AllChem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

import logging
logger = logging.getLogger('molsani')

def remove_invalid_SMILES(df):

    # Log rows where 'mol' is None before dropping
    invalid_rows = df[df['mol'].isna()]
    for index, row in invalid_rows.iterrows():
        logger.warning(f"INVALID SMILES: Removing row at index {index}: {row.to_dict()}")

    # Remove rows where 'mol' is None
    df_cleaned = df.dropna(subset=['mol'])

    return df_cleaned


def stripSMILESsalt(mol, molRemover):

    res, deleted = molRemover.StripMolWithDeleted(mol)

    if len(deleted) > 0: logger.info(f"Stripped salt {Chem.MolToSmiles(mol)}:  Salts:{' '.join([Chem.MolToSmiles(m) for m in deleted])}")

    return res

def removesalts(df: pd.DataFrame) -> pd.DataFrame:
    
    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'salt_stripping.txt'
    logger.info(f'Parsing salts SMARTS file: {smartsFile.resolve()}')

    remover = SaltRemover.SaltRemover(defnFilename=smartsFile)

    df['mol'] = df['mol'].apply(lambda x: stripSMILESsalt(x, remover))
    df['smiles'] = df['mol'].apply(lambda x: Chem.MolToSmiles(x))
    df=df[df['smiles']!=''] #Remove purely salt molecules
    return df


def detectPAINS(mol, catalog):
    entry = catalog.GetFirstMatch(mol)  # Get the first matching PAINS
    if entry is not None:
        return ('PAINS violation: ' + str(entry.GetDescription().capitalize()))  # Indicate a PAINS violation
    else:
        return 'OK'

def pains(df: pd.DataFrame, rejectedFile) -> pd.DataFrame:
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
    
    logger.info(f'Loading SMARTS from: {smartsFile}')
    if unwanted_option is not None: 
        # Using default substructure file 
        smarts_df = pd.read_csv(smartsFile, sep='\s+', names=['smarts','label', 'reason', 'mode', 'ref'])
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
    logger.info(f'Parsed {len(smarts_df)} substructures')
    return smarts_df

def filterbysmarts(mol, smarts_df: pd.DataFrame) -> str:
    for _, substructure in smarts_df.iterrows():
        if mol.HasSubstructMatch(substructure.mol):
            return substructure.label
    return 'OK'

def unwanted(df: pd.DataFrame, rejectedFile, unwanted_option) -> pd.DataFrame:

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

def custom(df: pd.DataFrame, rejectedFile, smartsFile) -> pd.DataFrame:

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


def stereoisomers(df: pd.DataFrame, max_isomers = 0) -> pd.DataFrame:
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
        isomers = generate_stereoisomers(row['mol'], max_isomers)
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
    reactions = []
    with open(file_path, 'r') as file:
        for line in file:
            smarts = line.strip().split()
            if smarts:
                reactions.append(AllChem.ReactionFromSmarts(smarts[0]))
    return reactions

# Function to apply reactions to a molecule
def recursive_reaction(mol, reactions):
    product = mol
    for rxn in reactions:
        outcomes = rxn.RunReactants((mol,))
        if len(outcomes) >= 1:
            product = outcomes[0][0]
            Chem.SanitizeMol(product)
            #print(Chem.MolToSmiles(product))
            return recursive_reaction(product, reactions)
    return product

def protonation(df: pd.DataFrame) -> pd.DataFrame:

    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'ionizations.txt'
    logger.info(f'Parsing ionization SMARTS file: {smartsFile.resolve()}')

    # Load reactions from file
    reactions = load_reactions(smartsFile)

    # Apply reactions to each SMILES in the DataFrame
    df['mol'] = df['mol'].apply(lambda x: recursive_reaction(x, reactions))

    return df



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
