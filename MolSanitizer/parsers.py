import argparse

import pandas as pd
from pathlib import Path

import logging
logger = logging.getLogger('molsani')



def parseArguments(args = None):
    info = """MolSanitizer - A package to prepare SMILES databases

    Ex. input file (space or tab-separated file):
        COCCC(=O)Nc1ncc(s1)Br  CP000000418470
        C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
        CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597

    Ex. run
    msani -i example.smi --removesalts --pains --unwanted all --stereoisomers --protonation
    msani -i example.smi --removesalts
    msani -i example.smi --pains --unwanted all --stereoisomers --protonation
    msani -i example.smi --pains --unwanted all --stereoisomers --protonation --db2
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description= info, formatter_class=argparse.RawTextHelpFormatter)
    
    # Add the required input files argument
    parser.add_argument('-i', '--input_files', type=str, nargs='+', help='Input files containing chemical structures')
    parser.add_argument('-s', '--smiles', default=None, type=str, nargs='+', help='Input SMILES strings')
    # Add Boolean options
    parser.add_argument('-e', '--enamine', action='store_true', help='Enamine input format (default: False)')
    parser.add_argument('--lazy', action='store_true', help='Implement all the processing and preparation steps (default: False)')
    parser.add_argument('--removesalts', action='store_true', help='Remove salts from the structures')
    parser.add_argument('--tautomers', action='store_true', help='Tautomers enumeration')
    parser.add_argument('--pains', action='store_true', help='Remove PAINS violations from the structures')
    parser.add_argument('--unwanted', choices=['all', 'regular', 'special', 'optional'], default=None, nargs='*', help='Filter out unwanted substructures using the default list')
    parser.add_argument('--create_custom', action='store_true', help='Generate a template for customized substructure filtering')
    parser.add_argument('--stereoisomers', action='store_true', help='Stereoisomers enumeration (only consider unspecified chiral centers)')
    parser.add_argument('--protonation', action='store_true', help='Apply protonation to the structures')
    parser.add_argument('--db2', action='store_true', help='Generate conformers and stored in the DB2 format for DOCK 3.8')
    parser.add_argument('--nocleanup', action='store_false', dest='cleanup', default = True, help='Do not clean up the temporary files (default: False)')

    #parser.add_argument('--neutralize', action='store_true', help='Neutralize the structures')
    #parser.add_argument('--flavioFilters', action='store_true', help='Filter using Flavio script (For databases based on Greg Landrum)')
    parser.add_argument('--debug', action='store_true', help='Debugging mode')
    parser.add_argument('--test', action='store_true', help='Test mode (silent mode)')

    # Add string option
    parser.add_argument('--custom', default=None, type=str, help='Filter out unwanted substructures using the customized list. To generate an example list, use --create_custom')
    parser.add_argument('-pre', '--prefix', default=None, type=str, help='Prefix for the output files. If not provided, the input file name will be used.')

    # Add integer option
    parser.add_argument('--max_isomers', type=int, default=0, help='Maximum number of stereoisomers to consider (default: 0 = no limit)')
    parser.add_argument('-nconf', '--numconfs', type=int, default=2000, help='Maximum number of conformers to generate (default: 2000)')
    parser.add_argument('-rs', '--randomSeed', type=int, default=42, help='Random seed for reproducibility (default: 42)')

    # Parse the arguments
    args = parser.parse_args()
    if args.input_files is not None and args.smiles is None:
        for inFile in args.input_files:
            if not Path(inFile).is_file():
                parser.error(f'The input file: {inFile} does not exist.')

    return args

def parseArguments_batch(args = None):
    info = """MolSanitizer - A package to prepare SMILES databases
    This is a batch version of the MolSanitizer package. It reads a list of input files and processes them in parallel.

    Ex. input file (space or tab-separated file):
        COCCC(=O)Nc1ncc(s1)Br  CP000000418470
        C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
        CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597

    Ex. run
    msani_batch -i example.smi --lazy --enamine --db2
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description= info, formatter_class=argparse.RawTextHelpFormatter)
    
    # Add the required input files argument
    parser.add_argument('-i', '--input_files', type=str, nargs='+', help='Input files containing chemical structures')

    # Add Boolean options
    parser.add_argument('-e', '--enamine', action='store_true', help='Enamine input format (default: False)')
    parser.add_argument('--lazy', action='store_true', help='Implement all the processing and preparation steps (default: False)')
    parser.add_argument('--removesalts', action='store_true', help='Remove salts from the structures')
    parser.add_argument('--tautomers', action='store_true', help='Tautomers enumeration')
    parser.add_argument('--pains', action='store_true', help='Remove PAINS violations from the structures')
    parser.add_argument('--unwanted', choices=['all', 'regular', 'special', 'optional'], default=None, nargs='*', help='Filter out unwanted substructures using the default list')
    parser.add_argument('--stereoisomers', action='store_true', help='Stereoisomers enumeration (only consider unspecified chiral centers)')
    parser.add_argument('--protonation', action='store_true', help='Apply protonation to the structures')
    parser.add_argument('--db2', action='store_true', help='Generate conformers and stored in the DB2 format for DOCK 3.8')
    parser.add_argument('--nocleanup', action='store_false', dest='cleanup', default = True, help='Do not clean up the temporary files (default: False)')

    # Add string option
    parser.add_argument('--custom', default=None, type=str, help='Filter out unwanted substructures using the customized list')
    parser.add_argument('-n', '--projectName', default='naiss2023-3-39', dest='proj_name', type=str, help='Project name for the SLURM script (default: naiss2023-3-39)')

    # Add integer option
    parser.add_argument('-l', '--lines_per_job', dest='lines', type=int, default=1000, help='Number of lines to process per job (default: 1000)')
    parser.add_argument('--max_isomers', type=int, default=0, help='Maximum number of stereoisomers to consider (default: 0 = no limit)')
    parser.add_argument('-nconf', '--numconfs', type=int, default=2000, help='Maximum number of conformers to generate (default: 2000)')
    parser.add_argument('-rs', '--randomSeed', type=int, default=42, help='Random seed for reproducibility (default: 42)')
    parser.add_argument('-t', '--time', type=int, default=24, help='Time limit for the SLURM job in hours (default: 24H)')
    parser.add_argument('--max_jobs', type=int, default=100, help='Maximum number of jobs to run simultaneously (default: 100)')

    # Parse the arguments
    args = parser.parse_args()
    for inFile in args.input_files:
        if not Path(inFile).is_file():
            parser.error(f'The input file: {inFile} does not exist.')

    return args
