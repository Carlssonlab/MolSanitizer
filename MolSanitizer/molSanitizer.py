import logging
logger = logging.getLogger('molsani')

import pandas as pd

import pathlib
import os
import time
import sys

from . import parsers
from . import loggers
from . import filters

from rdkit import Chem
from rdkit import rdBase

def process_chunk(chunk):
    # Split the 'smiles' column by space
    chunk['smiles'] = chunk['smiles'].apply(lambda x: x.split()[0])
    return chunk

def cleanData(args):

    start_time = time.time()

    if args.prefix is not None:
            outputFile = f"{args.prefix}_clean.txt"
            rejectedFile = f"{args.prefix}_rejected.txt"
            if os.path.exists(outputFile): os.remove(outputFile)
    logger.info(f'Rdkit version: {rdBase.rdkitVersion}')

    for inputFile in args.input_files:
        if args.enamine: 
            df_input = pd.read_csv(inputFile, sep='\t', names=['smiles', 'ids'], usecols=[0,1], header=None, chunksize=500_000)
            logger.info(f'Using Enamine format for parsing')
        else: 
            df_input = pd.read_csv(inputFile, sep=r'\s+', names=['smiles', 'ids'], usecols=[0,1], header=None, chunksize=1_000)
        
        inputFilePath = pathlib.Path(inputFile)
        logger.info(f'Processing: {inputFile}')

        if args.prefix is None:
            outputFile = inputFilePath.with_name(f"{inputFilePath.stem}_clean{inputFilePath.suffix}")
            rejectedFile = inputFilePath.with_name(f"{inputFilePath.stem}_rejected{inputFilePath.suffix}")
            if os.path.exists(outputFile): os.remove(outputFile)

        for step, chunk in enumerate(df_input, start=1):
            if args.enamine: chunk = process_chunk(chunk)
            chunk['mol'] = chunk['smiles'].apply(lambda x: Chem.MolFromSmiles(x))

            chunk = filters.remove_invalid_SMILES(chunk)
            # Remove salts
            if args.removesalts: chunk = filters.removesalts(chunk, args.debug)

            # Tautomers enumeration
            if args.tautomers: chunk = filters.tautomers(chunk, args.debug)

            # PAINS functional groups filtering
            if args.pains: chunk = filters.pains(chunk, rejectedFile, args.debug) 

            # Unwanted substructures filtering
            if args.unwanted is not None: chunk = filters.unwanted(chunk, rejectedFile, args.unwanted, args.debug) 
            if args.custom is not None: chunk = filters.custom(chunk, rejectedFile, args.custom, args.debug) 

            # Stereoisomers enumeration
            if args.stereoisomers: chunk = filters.stereoisomers(chunk, args.max_isomers, args.debug)

            # Protonation
            if args.protonation: chunk = filters.protonation(chunk, args.debug)

            #if args.standarizeFilters: chunk = filters.standarizeFilters(chunk)

            if not (args.test): print(f"Step {step} took {time.time() - start_time:.2f} seconds to complete.")
            chunk['smiles'] = chunk['mol'].apply(lambda x: Chem.MolToSmiles(x))
            chunk.to_csv(outputFile, index=False, mode='a', columns=['smiles','ids'], header=False, sep=' ')

def Sanitycheck(args: dict):
    """Sanity check for the unwanted flag

    Args:
        args (dict): Arguments from the command line

    Returns:
        dict: The updated arguments
    """
    if args.lazy:
        args.removesalts = True
        args.tautomers = True
        args.pains = True
        args.unwanted = ['all']
        args.stereoisomers = True
        args.protonation = True

    if args.unwanted is not None:
        if not args.unwanted: args.unwanted=['regular']
        if 'all' in args.unwanted: args.unwanted=['regular','special','optional']
        args.unwanted=[word.title() for word in args.unwanted]
    return args

def generateCustomTemplate(args):
    """Generate the custom template for substructure filtering 
    by copying the default template to the current directory.
    """
    file = os.path.join(os.path.dirname(__file__), 'Data', 'filter_out.csv')
    if args.prefix is not None: os.system(f"cp {file} {args.prefix}.tsv") 
    else: os.system(f"cp {file} template.tsv")
    if not (args.test): 
        print(f"Generated template substructure list as template.tsv")
        print(f"The first two columns (SMARTS and LABEL) are required for substructure filtering.")
        print(f"Other arguments are skipped, the program exitted normally.")

def main():

    args = parsers.parseArguments(sys.argv[1:])
    args = Sanitycheck(args)

   
    if args.create_custom: 
        generateCustomTemplate(args)
    else:
        input_path = pathlib.Path(args.input_files[0])
        if args.prefix is not None: log_file = f'{args.prefix}.log' 
        else: log_file = input_path.with_suffix('.log')
        loggers.setup_logger(log_file)
        original_command = ' '.join(sys.argv)
        logger.info(f"STARTING MOLSANITIZER")
        logger.info(f"Input: {original_command}")    
        loggers.arguments(args)
        cleanData(args)



if __name__=="__main__":
    main()