from pathlib import Path

import pandas as pd

from rdkit import Chem, AllChem
from rdkit.Chem import SaltRemover

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


def removeSalts(df: pd.Dataframe) -> pd.DataFrame:

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

def cleanFilter(df: pd.DataFrame) -> pd.DataFrame:

    pass