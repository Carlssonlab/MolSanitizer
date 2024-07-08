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

def cleanData(args):

    start_time = time.time()

    for inputFile in args.input_files:

        df_input = pd.read_csv(inputFile, sep=r'\s+', names=['smiles', 'ids'], header=None, chunksize=500_000)

        inputFilePath = pathlib.Path(inputFile)
        outputFile = inputFilePath.with_name(f"{inputFilePath.stem}_clean{inputFilePath.suffix}")
        rejectedFile = inputFilePath.with_name(f"{inputFilePath.stem}_rejected{inputFilePath.suffix}")

        if os.path.exists(outputFile): 
            os.remove(outputFile)

        logger.info(f'Rdkit version: {rdBase.rdkitVersion}')

        for step, chunk in enumerate(df_input, start=1):

            chunk['mol'] = chunk['smiles'].apply(lambda x: Chem.MolFromSmiles(x))

            chunk = filters.remove_invalid_SMILES(chunk)
            # Remove salts
            if args.removesalts: chunk = filters.removesalts(chunk)

            # Tautomers enumeration
            if args.tautomers: chunk = filters.tautomers(chunk)

            # PAINS functional groups filtering
            if args.pains: chunk = filters.pains(chunk, rejectedFile) 

            # Unwanted substructures filtering
            if args.unwanted is not None: chunk = filters.unwanted(chunk, rejectedFile, args.unwanted) 
            if args.custom is not None: chunk = filters.custom(chunk, rejectedFile, args.custom) 


            # Stereoisomers enumeration
            if args.stereoisomers: chunk = filters.stereoisomers(chunk, args.max_isomers)

            # Protonation
            if args.protonation: chunk = filters.protonation(chunk)

            #if args.standarizeFilters: chunk = filters.standarizeFilters(chunk)

            print(f"Step {step} took {time.time() - start_time:.2f} seconds to complete.")
            chunk['smiles'] = chunk['mol'].apply(lambda x: Chem.MolToSmiles(x))
            chunk.to_csv(outputFile, index=False, mode='a', columns=['smiles','ids'], header=False, sep=' ')

def Sanitycheck(args: dict):
    if args.unwanted is not None:
        if not args.unwanted: args.unwanted=['regular']
        if 'all' in args.unwanted: args.unwanted=['regular','special','optional']
        args.unwanted=[word.title() for word in args.unwanted]
    return args

def generateCustomTemplate():
    file = os.path.join(os.path.dirname(__file__), 'Data', 'filter_out.csv')
    os.system(f"cp {file} template.tsv") 
    print(f"Generated template substructure list as template.tsv")
    print(f"The first two columns (SMARTS and LABEL) are required for substructure filtering.")
    print(f"Other arguments are skipped, the program exitted normally.")
    exit()

def main():

    args = parsers.parseArguments()
    args = Sanitycheck(args)

   
    if args.create_custom: 
        generateCustomTemplate()
    else:
        input_path = pathlib.Path(args.input_files[0])
        log_file = input_path.with_suffix('.log')
        loggers.setup_logger(log_file)
        original_command = ' '.join(sys.argv)
        logger.info(f"Starting MolSanitizer: {original_command}")    
        loggers.arguments(args)
        cleanData(args)



if __name__=="__main__":
    main()