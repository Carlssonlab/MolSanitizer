import logging
logger = logging.getLogger('molsani')

import pandas as pd

import pathlib
import os

from . import parsers
from . import loggers


def cleanData(args):

    saltStripping_df, cleanFilter_df, protonation_df = parsers.parseDatabases(args)

    for inputFile in args.input_files:

        df_input = pd.read_csv(inputFile, sep=r'\s+', names=['smarts', 'ids'], header=None, chunksize=1_000_000)

        inputFilePath = pathlib.Path(inputFile)
        outputFile = inputFilePath.with_name(f"{inputFilePath.stem}_clean{inputFilePath.suffix}")

        if os.path.exists(outputFile): os.remove(outputFile)

        for chunk in df_input:

            #DO something with the chunk

            print(chunk)

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