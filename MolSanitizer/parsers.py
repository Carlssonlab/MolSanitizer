import argparse

import pandas as pd
from pathlib import Path

import logging
logger = logging.getLogger('molsani')



def parseArguments():
    info = """MolSanitizer - A package to prepare SMILES databases

    Ex. input file (space or tab-separated file):
        COCCC(=O)Nc1ncc(s1)Br  CP000000418470
        C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
        CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597

    Ex. run
    msani -i example.smi --removesalts --pains --unwanted all --stereoisomers --protonation
    msani -i example.smi --removesalts
    msani -i example.smi --pains --unwanted all --stereoisomers --protonation
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description= info, formatter_class=argparse.RawTextHelpFormatter)
    
    # Add the required input files argument
    parser.add_argument('-i', '--input_files', type=str, nargs='+', help='Input files containing chemical structures')

    # Add Boolean options
    parser.add_argument('--removesalts', action='store_true', help='Remove salts from the structures')
    parser.add_argument('--tautomers', action='store_true', help='Tautomers enumeration')
    parser.add_argument('--pains', action='store_true', help='Remove PAINS violations from the structures')
    parser.add_argument('--unwanted', choices=['all', 'regular', 'special', 'optional'], default=None, nargs='*', help='Filter out unwanted substructures using the default list')
    parser.add_argument('--custom', default=None, type=str, help='Filter out unwanted substructures using the customized list. To generate an example list, use --create_custom')
    parser.add_argument('--create_custom', action='store_true', help='Generate a template for customized substructure filtering')
    parser.add_argument('--stereoisomers', action='store_true', help='Stereoisomers enumeration (only consider unspecified chiral centers)')
    parser.add_argument('--protonation', action='store_true', help='Apply protonation to the structures')
    parser.add_argument('--neutralize', action='store_true', help='Neutralize the structures')
    #parser.add_argument('--flavioFilters', action='store_true', help='Filter using Flavio script (For databases based on Greg Landrum)')
    parser.add_argument('--debug', action='store_true', help='Debugging mode')

    # Add integer option
    parser.add_argument('--max_isomers', type=int, default=0, help='Maximum number of tautomers to consider (default: 0 = no limit)')
    
    # Parse the arguments
    args = parser.parse_args()
    
    return args

def template(templateFile: str, label: str) -> pd.DataFrame:

    # Get the absolute path to the template SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / templateFile
    logger.info(f'Parsing {label} SMARTS file: {smartsFile.resolve()}')

    df = pd.read_csv(smartsFile, sep=r'\s+', names=['smarts', 'ids'], header=None)

    return df


def parseDatabases(args):

    protonation_df = None
    cleanFilter_df = None

    if args.protonation: 
        protonation_df = template('ionizations.txt', 'ionization')
    if args.cleanFilter: 
        cleanFilter_df = template('filter-out.txt', 'CleanFilter')


    return cleanFilter_df, protonation_df