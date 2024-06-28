import logging
logger = logging.getLogger('molsani')

import pandas as pd

import pathlib
import os
import time

from . import parsers
from . import loggers
from . import filters

from rdkit import Chem
from rdkit import rdBase

def cleanData(args):

    start_time = time.time()

    for inputFile in args.input_files:

        df_input = pd.read_csv(inputFile, sep=r'\s+', names=['smiles', 'ids'], header=None, chunksize=1_000_000)

        inputFilePath = pathlib.Path(inputFile)
        outputFile = inputFilePath.with_name(f"{inputFilePath.stem}_clean{inputFilePath.suffix}")

        if os.path.exists(outputFile): 
            os.remove(outputFile)

        logger.info(rdBase.rdkitVersion)

        for step, chunk in enumerate(df_input, start=1):

            chunk['mol'] = chunk['smiles'].apply(lambda x: Chem.MolFromSmiles(x))

            chunk = filters.remove_invalid_SMILES(chunk)

            if args.removeSalts: chunk = filters.removeSalts(chunk)

            if args.protonation: chunk = filters.protonation(chunk)

            if args.reactivityFilter: chunk = filters.reactivityFilter(chunk)

            if args.standarizeFilters: chunk = filters.standarizeFilters(chunk)

            print(f"Step {step} took {time.time() - start_time:.2f} seconds to complete.")

            chunk.to_csv(outputFile, index=False, mode='a', header=False, sep=' ')


def main():

    args = parsers.parseArguments()

    input_path = pathlib.Path(args.input_files[0])
    log_file = input_path.with_suffix('.log')
    loggers.setup_logger(log_file)
    
    loggers.arguments(args)


    cleanData(args)



if __name__=="__main__":
    main()