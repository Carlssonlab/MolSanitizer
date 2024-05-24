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
    parser.add_argument('--protonation', action='store_true', help='Apply protonation to the structures')
    parser.add_argument('--removeSalts', action='store_true', help='Remove salts from the structures')
    parser.add_argument('--neutralize', action='store_true', help='Neutralize the structures')
    parser.add_argument('--tautomers', action='store_true', help='Consider tautomers for the structures')
    parser.add_argument('--cleanFilter', action='store_true', help='Filter BAD molecules (PAINS, Reactive,...)')
    
    # Add integer option
    parser.add_argument('--maxTautomers', type=int, default=0, help='Maximum number of tautomers to consider (default: 0)')
    
    # Parse the arguments
    args = parser.parse_args()
    
    return args

def template(templateFile: str, label: str) -> pd.DataFrame:

    # Get the absolute path to the protonation SMARTS file using pathlib
    smartsFile = Path(__file__).parent / 'Data' / templateFile
    logger.info(f'Parsing {label} SMARTS file: {smartsFile.resolve()}')

    df = pd.read_csv(smartsFile, sep=r'\s+', names=['smarts', 'ids'], header=None)

    return df


def parseDatabases(args):

    protonation_df = None
    saltStripping_df = None
    cleanFilter_df = None

    if args.protonation: 
        protonation_df = template('ionizations.txt', 'ionization')
    if args.removeSalts: 
        saltStripping_df = template('salt_stripping.txt', 'removeSalts')
    if args.cleanFilter: 
        cleanFilter_df = template('filter-out.txt', 'CleanFilter')


    return saltStripping_df, cleanFilter_df, protonation_df