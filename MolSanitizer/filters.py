from pathlib import Path

import pandas as pd

from rdkit import Chem

from rdkit.Chem import AllChem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

import logging
logger = logging.getLogger('molsani')

def stripSMILESsalt(mol, molRemover):

    res, deleted = molRemover.StripMolWithDeleted(mol)

    if len(deleted) > 0: logger.info(f"Stripped salt {Chem.MolToSmiles(mol)}:  Salts:{' '.join([Chem.MolToSmiles(m) for m in deleted])}")

    return res

def remove_invalid_SMILES(df):

    # Log rows where 'mol' is None before dropping
    invalid_rows = df[df['mol'].isna()]
    for index, row in invalid_rows.iterrows():
        logger.warning(f"INVALID SMILES: Removing row at index {index}: {row.to_dict()}")

    # Remove rows where 'mol' is None
    df_cleaned = df.dropna(subset=['mol'])

    return df_cleaned


def removeSalts(df: pd.DataFrame) -> pd.DataFrame:

    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'salt_stripping.txt'
    logger.info(f'Parsing salts SMARTS file: {smartsFile.resolve()}')

    remover = SaltRemover.SaltRemover(defnFilename=smartsFile)

    df['mol'] = df['mol'].apply(lambda x: stripSMILESsalt(x, remover))

    return df

def protonation(df: pd.DataFrame) -> pd.DataFrame:

    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'ionizations.txt'
    logger.info(f'Parsing ionization SMARTS file: {smartsFile.resolve()}')

    # Load reactions from file
    reactions = load_reactions(smartsFile)

    # Apply reactions to each SMILES in the DataFrame
    df['mol'] = df['mol'].apply(lambda x: apply_reactions(x, reactions))

    return df


# Function to read SMARTS reactions from a file
def load_reactions(file_path):
    reactions = []
    with open(file_path, 'r') as file:
        for line in file:
            smarts = line.strip()
            if smarts:
                reactions.append(AllChem.ReactionFromSmarts(smarts))
    return reactions

# Function to apply reactions to a molecule
def apply_reactions(mol, reactions):

    products = set()
    for rxn in reactions:
        outcomes = rxn.RunReactants((mol,))
        for outcome in outcomes:
            for prod in outcome:
                products.add(Chem.MolToSmiles(prod))
    return list(products)


def loadSMARTSdata(smartsFile: str):
    
    logger.info(f'Loading SMARTS from: {smartsFile}')
    
    smarts_df = pd.read_csv(smartsFile, delim_whitespace=True, header=None)
    
    smarts_lst = smarts_df[0].to_list()
    
    smarts_names = smarts_df[1].to_list()
    
    if smarts_df.shape[0] > 1: smarts_names = smarts_df[1].to_list()

    matchers = [Chem.MolFromSmarts(x, mergeHs=True) for x in smarts_lst]
    
    painsDefs = [[smart, id] for smart,id in zip(smarts_lst, smarts_names)]
    
    return matchers, painsDefs

def reactivityFilter(df: pd.DataFrame) -> pd.DataFrame:

    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / 'filter_out.csv'
    logger.info(f'Parsing reactive filters SMARTS file: {smartsFile.resolve()}')

    # Load smarts to clean  from file
    matchers, painsDefs = loadSMARTSdata(smartsFile.resolve())

    # Apply reactions to each SMILES in the DataFrame
    # df['mol'] = df['mol'].apply(lambda x: apply_reactions(x, reactions))

    return df


def generateFilteredFile(matchers: list, painsDefs: list, filtered_dict: dict):

    painsLOGdata = []
    filtered_pains = {}
    
    my_iterator = filtered_dict.items()
    
    for id, data in my_iterator:
        
        try:
    
            moleculeRDkit = Chem.MolFromSmiles(data)
            
        except:

            logger.info(f'Problem generating molecule RDkit: {id} ')
            continue

        if moleculeRDkit == None: 
            
                logger.warning(f'Molecule ID {id} has some structural problem. Molecule removed')
                continue

        matchFound = False
        painDef = ''
        
        for matcher in matchers:
            
            if moleculeRDkit.HasSubstructMatch(matcher):
                
                matchFound = True
                painDef = painsDefs[matchers.index(matcher)][1].replace('regId=', '')
                break                

        if matchFound: painsLOGdata.append(f"{' '.join([str(iddata) for iddata in id])} {painDef}")
        else: filtered_pains[id] = data
                        
            
    return filtered_pains, painsLOGdata


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
