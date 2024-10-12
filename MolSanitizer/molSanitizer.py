"""
MolSanitizer.
"""

__author__ = "Thua-Phong Lam, Israel Cabeza de Vaca Lopez, Szymon Pach"
__place__ = "Jens Carlsson lab, Uppsala University, Sweden"
__license__ = "MIT"
__version__ = "0.1.3"

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
from . import smi2db2

from rdkit import Chem
from rdkit import rdBase

def process_enamine_name(chunk):
    # Split the 'smiles' column by space
    chunk['smiles'] = chunk['smiles'].apply(lambda x: x.split()[0])
    return chunk

def apply_filters(chunk, args, rejected_file):
    chunk = filters.remove_invalid_SMILES(chunk)

    if args.removesalts:
        chunk = filters.removesalts(chunk, args.debug)
    if args.tautomers:
        chunk = filters.tautomers(chunk, args.debug)
    if args.pains:
        chunk = filters.pains(chunk, rejected_file, args.debug)
    if args.unwanted:
        chunk = filters.unwanted(chunk, rejected_file, args.unwanted, args.debug)
    if args.custom:
        chunk = filters.custom(chunk, rejected_file, args.custom, args.debug)
    if args.protonation:
        chunk = filters.protonation(chunk, args.debug)
    if args.stereoisomers:
        chunk = filters.stereoisomers(chunk, args.max_isomers, args.numcores, args.debug)

    return chunk

def log_step_time(elapsed_time, step):
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours != 0:
        print(f"Step {step} took {int(hours):02}:{int(minutes):02}:{int(seconds):02} hours to complete.")
    elif minutes != 0:
        print(f"Step {step} took {int(minutes):02}:{int(seconds):02} minutes to complete.")
    else:
        print(f"Step {step} took {elapsed_time:.2f} seconds to complete.")

def log_execution_time(start_time, is_test):
    if not is_test:
        elapsed_time = time.time() - start_time
        minutes, seconds = divmod(elapsed_time, 60)
        if minutes != 0:
            print(f"MolSanitizer took {int(minutes):02}:{int(seconds):02} minutes to complete.")
        else:
            print(f"MolSanitizer took {elapsed_time:.2f} seconds to complete.")

def get_output_files(args, input_file_path):
    if args.prefix:
        output_file = f"{args.prefix}_clean.txt"
        rejected_file = f"{args.prefix}_rejected.txt"
    else:
        output_file = input_file_path.with_name(f"{input_file_path.stem}_clean{input_file_path.suffix}")
        rejected_file = input_file_path.with_name(f"{input_file_path.stem}_rejected{input_file_path.suffix}")

    if os.path.exists(output_file):
        os.remove(output_file)

    return output_file, rejected_file

def read_input_file(input_file, is_enamine):
    if is_enamine:
        logger.info('Using Enamine format for parsing')
        return pd.read_csv(input_file, sep='\t', names=['smiles', 'ids'], usecols=[0, 1], header=None, chunksize=500_000)
    else:
        return pd.read_csv(input_file, sep=r'\s+', names=['smiles', 'ids'], usecols=[0, 1], header=None, chunksize=500_000)

def process_files(args, start_time: int):
    for input_file in args.input_files:
        input_file_path = pathlib.Path(input_file)
        logger.info(f'Processing: {input_file}')

        output_file, rejected_file = get_output_files(args, input_file_path)

        df_input = read_input_file(input_file, args.enamine)

        for step, chunk in enumerate(df_input, start=1):
            if args.enamine:
                chunk = process_enamine_name(chunk)
            
            chunk['mol'] = chunk['smiles'].apply(Chem.MolFromSmiles)
            chunk['ids'] = chunk['ids'].astype(str)

            chunk = apply_filters(chunk, args, rejected_file)
            
            if not chunk.empty:
                chunk.to_csv(output_file, index=False, mode='a', columns=['smiles', 'ids'], header=False, sep=' ')


            if args.db2:
                smi2db2.gen_conf_chunk_ver2(chunk, args)
                       
            if not args.test:
                if step == 1: time_step1 = time.time()-start_time
                if step == 2:
                    log_step_time(time_step1, 1)
                    log_step_time(time.time()-start_time, 2)
                elif step > 2:
                    log_step_time(time.time()-start_time, step)
                start_time = time.time()

def process_smiles(args):
    rejected_file = "msani_rejected.txt"
    chunk = pd.DataFrame({'smiles': args.smiles, 'ids': range(len(args.smiles))})
    chunk['ids'] = chunk['ids'].astype(str)
    chunk['mol'] = chunk['smiles'].apply(Chem.MolFromSmiles)
    
    chunk = apply_filters(chunk, args, rejected_file)

    if args.db2: 
        if os.path.exists('db2/0.db2'): os.remove('db2/0.db2') # 0 is the default name
        smi2db2.gen_conf_chunk_ver2(chunk, args)
    
    print('Processed SMILES:')
    for i, row in chunk.iterrows():
        print(row['smiles'])

def clean_data(args):
    start_time = time.time()
    rdkit_version = rdBase.rdkitVersion
    logger.info(f'RDKit version: {rdkit_version}')
    if rdkit_version != '2024.09.1':
        print('\n########################################################')
        print('RDKit version 2024.09.1 is recommended for MolSanitizer.')
        print("Use 'conda install rdkit==2024.9.1' to avoid potential issues.")
        print('########################################################\n')
        
    if args.smiles:
        process_smiles(args)
    else:
        process_files(args, start_time)

    log_execution_time(start_time, args.test)


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
    if args.db2:
        # Always enumerate stereoisomers for before generating DB2 files
        # Maximum number of stereoisomers is set to in parser
        args.stereoisomers = True
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
        if args.input_files is not None and args.smiles is None:
            input_path = pathlib.Path(args.input_files[0])
            if args.prefix is not None: log_file = f'{args.prefix}.log' 
            else: log_file = input_path.with_suffix('.log')
        else: log_file = 'molsani.log'
        loggers.setup_logger(log_file)
        original_command = ' '.join(sys.argv)
        logger.info(f"#######  STARTING MOLSANITIZER  #######")
        logger.info(f"Input: {original_command}")    
        loggers.arguments(args)
        clean_data(args)
        logger.info(f"*******  MOLSANITIZER FINISHED *******")



if __name__=="__main__":
    main()