import argparse

import pandas as pd
from pathlib import Path

import logging
logger = logging.getLogger('molsani')

def parseArguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Process some chemical structures.")
    
    # Add the required input files argument
    parser.add_argument('input_files', metavar='input_files', type=str, nargs='+', help='Input files containing chemical structures')

    # Add Boolean options
    parser.add_argument('--removeSalts', action='store_true', help='Remove salts from the structures')
    parser.add_argument('--PAINSFilter', action='store_true', help='Remove PAINS violations from the structures')
    # TODO: for this argument, let the user choose which reactive functional group to retain (Target-specific)
    parser.add_argument('--reactivityFilter', action='store_true', help='Filter out reactive functional groups')
    parser.add_argument('--tautomers', action='store_true', help='Consider tautomers for the structures')
    parser.add_argument('--protonation', action='store_true', help='Apply protonation to the structures')
    parser.add_argument('--neutralize', action='store_true', help='Neutralize the structures')
    parser.add_argument('--flavioFilters', action='store_true', help='Filter using Flavio script (For databases based on Greg Landrum)')
    
    # Add integer option
    parser.add_argument('--maxTautomers', type=int, default=0, help='Maximum number of tautomers to consider (default: 0)')
    
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