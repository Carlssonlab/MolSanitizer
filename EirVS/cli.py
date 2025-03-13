"""
EirVS in a standalone mode.
"""

__author__ = "Thua-Phong Lam, Israel Cabeza de Vaca Lopez, Szymon Pach"
__place__ = "Jens Carlsson lab, Uppsala University, Sweden"
__license__ = "GPLv2"
__version__ = "0.2.3"

import logging
logger = logging.getLogger('eirvs')

import pandas as pd

import pathlib
import os
import time
import sys

from . import parsers
from . import loggers
from . import smiles_sanitizer
from . import smi2db2

from rdkit import Chem
from rdkit import rdBase

def process_enamine_name(chunk):
    # Split the 'smiles' column by space
    chunk['smiles'] = chunk['smiles'].apply(lambda x: x.split()[0])
    return chunk

def apply_filters(chunk, args, rejected_file):
    sanitizer = smiles_sanitizer.SmilesSanitizer(
        removesalts=args.removesalts, custom= args.custom, unwanted=args.unwanted,
        pains=args.pains, ha=args.ha, logp=args.logp, hba=args.hba, hbd=args.hbd, 
        mw=args.mw, tautomers=args.tautomers, taurdkit=args.taurdkit, neutralize=args.neutralize,
        stereoisomers=args.stereoisomers, max_stereoisomers=args.max_stereoisomers,
        protonation=args.protonation, pH=args.pH, pH_range=args.pH_range, 
        numcores=args.numcores, conformal=args.conformal, db2 = args.db2, debug=args.debug)
    chunk = sanitizer.Sanitize(chunk, rejected_file)
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
            print(f"EirVS took {int(minutes):02}:{int(seconds):02} minutes to complete.")
        else:
            print(f"EirVS took {elapsed_time:.2f} seconds to complete.")

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

def read_input_file(input_file, is_enamine, is_synthon):
    if is_synthon:
        logger.info('Using Synthon format for parsing')
        return pd.read_csv(
            input_file,
            sep=r'\s+',
            names=['smiles', 'ids', 'highlights'],
            usecols=[0, 1, 2],
            header=None,
            chunksize=250_000,
            dtype={'smiles': str, 'ids': str, 'highlights': str}  # Enforce string types
        )
    if is_enamine:
        logger.info('Using Enamine format for parsing')
        return pd.read_csv(
            input_file,
            sep='\t',
            names=['smiles', 'ids'],
            usecols=[0, 1],
            header=None,
            chunksize=250_000,
            dtype={'smiles': str, 'ids': str}  # Enforce string types
        )
    else:
        return pd.read_csv(
            input_file,
            sep=r'\s+',
            names=['smiles', 'ids'],
            usecols=[0, 1],
            header=None,
            chunksize=250_000,
            dtype={'smiles': str, 'ids': str}  # Enforce string types
        )
    
def process_files(args, start_time: int):
    if args.conformal:
        logger.warning('Conformal predictor format preparation selected. Will skip all other flags and only standardize the molecules using RDKit default functions.')

    for input_file in args.input_files:
        input_file_path = pathlib.Path(input_file)
        logger.info(f'Processing: {input_file}')

        output_file, rejected_file = get_output_files(args, input_file_path)
        if os.path.exists(rejected_file): os.remove(rejected_file)
        
        df_input = read_input_file(input_file, args.enamine, args.synthon)

        for step, chunk in enumerate(df_input, start=1):
            if args.enamine: chunk = process_enamine_name(chunk)
            chunk = apply_filters(chunk, args, rejected_file)
            if not chunk.empty:
                if args.synthon and not(args.conformal):
                    chunk.to_csv(output_file, index=False, mode='a', columns=['smiles', 'ids', 'highlights'], header=False, sep=' ')
                else:
                    chunk.to_csv(output_file, index=False, mode='a', columns=['smiles', 'ids'], header=False, sep=' ')


            if args.db2:
                if args.corina:
                    smi2db2.gen_conf_chunk_corina(chunk, args, input_file_path.stem)
                else:
                    smi2db2.gen_conf_chunk(chunk, args, input_file_path.stem)
            
            if args.pdbqt:
                from . import smi2pdbqt
                smi2pdbqt.gen_conf_chunk(chunk, args)
            if not args.test:
                if step == 1: time_step1 = time.time()-start_time
                if step == 2:
                    log_step_time(time_step1, 1)
                    log_step_time(time.time()-start_time, 2)
                elif step > 2:
                    log_step_time(time.time()-start_time, step)
                start_time = time.time()

def process_smiles(args):
    rejected_file = "eirvs_rejected.txt"
    chunk = pd.DataFrame({'smiles': args.smiles, 'ids': range(len(args.smiles))})
    chunk['ids'] = chunk['ids'].astype(str)
    chunk['mol'] = chunk['smiles'].apply(Chem.MolFromSmiles)
    
    chunk = apply_filters(chunk, args, rejected_file)

    if args.db2: 
        if os.path.exists('db2/0.db2'): os.remove('db2/0.db2') # 0 is the default name
        if args.corina:
            smi2db2.gen_conf_chunk_corina(chunk, args)
        else:
            smi2db2.gen_conf_chunk(chunk, args)

    if args.pdbqt:
        from . import smi2pdbqt
        smi2pdbqt.gen_conf_chunk(chunk, args)
    
    print('Processed SMILES:')
    for i, row in chunk.iterrows():
        print(row['smiles'])

def clean_data(args):
    start_time = time.time()
    rdkit_version = rdBase.rdkitVersion
    logger.info(f'RDKit version: {rdkit_version}')
    if rdkit_version != '2024.09.1':
        print('\n########################################################')
        print('RDKit version 2024.09.1 is recommended for EirVS.')
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
    if args.db2 or args.pdbqt:
        # Always enumerate stereoisomers for before generating DB2 and PDBQT files
        # Maximum number of stereoisomers is set to in parser
        args.stereoisomers = True
    
    if (args.pH != 7 or args.pH_range != 0) and not args.protonation:
        print("It seems like you forget the --protonation flag. We turned it on for you.")
        args.protonation = True
    return args

def generateCustomTemplate(args):
    """Generate the custom template for substructure filtering 
    by copying the default template to the current directory.
    """
    file = os.path.join(os.path.dirname(__file__), 'Data', 'filter_out.csv')
    if args.prefix is not None: os.system(f"cp {file} {args.prefix}.txt") 
    else: os.system(f"cp {file} template.txt")
    if not (args.test): 
        print(f"Generated template substructure list as template.txt")
        print(f"The first two columns (SMARTS and LABEL) are required for substructure filtering.")
        print(f"Other arguments are skipped, the program exitted normally.")

def main():

    args = parsers.parseArguments(sys.argv[1:])
    args = Sanitycheck(args)
    if args.version:
        print(f"EirVS version: {__version__}")
        print(f"RDKit version: {rdBase.rdkitVersion}")
        return
   
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
        logger.info(f"#######  STARTING EirVS {__version__} #######")
        logger.info(f"{original_command}")    
        loggers.arguments(args)
        clean_data(args)
        logger.info(f"***********  EirVS FINISHED *****************")



if __name__=="__main__":
    main()