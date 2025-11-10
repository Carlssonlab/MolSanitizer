"""
MolSanitizer in a standalone mode.
"""

__author__ = "Thua-Phong Lam, Szymon Pach, Israel Cabeza de Vaca"
__place__ = "Jens Carlsson lab, Uppsala University, Sweden"
__license__ = "GPLv2"
__version__ = "0.5.0"

import pathlib
import os
import time
import sys
import logging
import shlex

from pandas import DataFrame, read_csv  # only what you use
from rdkit import Chem, rdBase

from msani.io import parsers, loggers
from msani.api import Msani  

logger = logging.getLogger('msani')

version_text = f"""Python version: {sys.version.split('|')[0]}
MolSanitizer version: {__version__}
RDKit version: {rdBase.rdkitVersion}"""

logo=r""" __  __         _  _____                _  _    _                 
|  \/  |       | |/  ___|              (_)| |  (_)                
| .  . |  ___  | |\ `--.   __ _  _ __   _ | |_  _  ____ ___  _ __ 
| |\/| | / _ \ | | `--. \ / _` || '_ \ | || __|| ||_  // _ \| '__|
| |  | || (_) || |/\__/ /| (_| || | | || || |_ | | / /|  __/| |   
\_|  |_/ \___/ |_|\____/  \__,_||_| |_||_| \__||_|/___|\___||_|   
"""


def process_enamine_name(chunk):
    # Split the 'smiles' column by space
    chunk['smiles'] = chunk['smiles'].apply(lambda x: x.split()[0])
    return chunk

def apply_processes(chunk, args, rejected_file):
    processor = Msani(
        removesalts=args.removesalts, custom= args.custom, unwanted=args.unwanted,
        pains=args.pains, ha=args.ha, logp=args.logp, hba=args.hba, hbd=args.hbd, 
        mw=args.mw, chiral = args.chiral, tautomers=args.tautomers, taurdkit=args.taurdkit, 
        neutralize=args.neutralize, stereoisomers=args.stereoisomers, 
        max_stereoisomers=args.max_isomers, protonation=args.protonation, pH=args.pH, 
        pH_range=args.pH_range, numcores=args.numcores, randomSeed=args.randomSeed, 
        standardize=args.standardize, protonation_library=args.protlib, tautomer_library=args.taulib,  
        debug=args.debug)
    chunk = processor.run(chunk, rejected_file)
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
        output_file = pathlib.Path(f"{input_file_path.stem}_clean{input_file_path.suffix}")
        rejected_file = pathlib.Path(f"{input_file_path.stem}_rejected{input_file_path.suffix}")

    if os.path.exists(output_file):
        os.remove(output_file)

    return output_file, rejected_file

def read_input_file(input_file, is_enamine, is_synthon):
    if is_synthon:
        logger.info('Using Synthon format for parsing')
        return read_csv(
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
        return read_csv(
            input_file,
            sep='\t',
            names=['smiles', 'ids'],
            usecols=[0, 1],
            header=None,
            chunksize=250_000,
            dtype={'smiles': str, 'ids': str}  # Enforce string types
        )
    else:
        return read_csv(
            input_file,
            sep=r'\s+',
            names=['smiles', 'ids'],
            usecols=[0, 1],
            header=None,
            chunksize=250_000,
            dtype={'smiles': str, 'ids': str}  # Enforce string types
        )
    
def process_files(args, start_time: int):
    if args.standardize:
        logger.warning('standardize predictor format preparation selected. Will skip all other flags and only standardize the molecules using RDKit default functions.')

    for input_file in args.input_files:
        input_file_path = pathlib.Path(input_file)
        logger.info(f'Processing: {input_file}')

        output_file, rejected_file = get_output_files(args, input_file_path)
        if os.path.exists(rejected_file): os.remove(rejected_file)
        
        df_input = read_input_file(input_file, args.extended, args.synthon)

        for step, chunk in enumerate(df_input, start=1):
            # if args.enamine: chunk = process_enamine_name(chunk)
            chunk = apply_processes(chunk, args, rejected_file)
            if not chunk.empty:
                if args.synthon and not(args.standardize):
                    chunk.to_csv(output_file,
                                 index=False,
                                 mode='a',
                                 columns=['smiles', 'ids', 'highlights'],
                                 header=False,
                                 sep=' ')
                else:
                    chunk.to_csv(output_file,
                                 index=False,
                                 mode='a',
                                 columns=['smiles', 'ids'],
                                 header=False,
                                 sep=' ')

            if args.gen3d:
                from msani.conformers import conformers
                conformers.gen_conf_chunk(chunk, args, input_file_path.stem)
            
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
    smiles_list, mols, names = [], [], []
    for idx, smiles in enumerate(args.smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return
        if mol.HasProp('_Name'): name = str(mol.GetProp('_Name')) 
        else: name = str(idx + 1)  # Use index as name if no _Name property
        smiles_list.append(smiles)
        mols.append(mol)
        names.append(name)
    chunk = DataFrame({'smiles': smiles_list, 'ids': names, 'mol': mols})
    
    chunk = apply_processes(chunk, args, rejected_file)

    if args.gen3d:
        from msani.conformers import conformers
        if os.path.exists('db2/0.db2'): os.remove('db2/0.db2') # 0 is the default name
        conformers.gen_conf_chunk(chunk, args)
    
    print('Processed SMILES:')
    for _, row in chunk.iterrows():
        print(f"{row['smiles']} {row['ids']}")

def clean_data(args):
    start_time = time.time()
    rdkit_version = rdBase.rdkitVersion
    logger.info(f'RDKit version: {rdkit_version}')        
    if args.smiles:
        process_smiles(args)
    elif args.input_files:
        process_files(args, start_time)

    log_execution_time(start_time, args.test)



def generateCustomTemplate(args, filename):
    """Generate the custom template for substructure filtering 
    by copying the default template to the current directory.
    """
    import pkgutil
    data = pkgutil.get_data('msani', f'Data/{filename}')
    if data is None:
        raise FileNotFoundError(f"Could not find Data/{filename} inside the binary.")
    if filename == 'argument_config.yaml': output_filename = 'config.yaml'
    else: output_filename = f"{args.prefix}.txt" if args.prefix else filename
    with open(output_filename, 'wb') as f:
        f.write(data)

    if not args.test:
        print(f"Generated template file: {output_filename}")

def main():
    if os.getenv('SLURM_JOB_ID') is None: print(logo)
    if len(sys.argv) == 1:
        print(version_text)
        print("No arguments provided. Use -h or --help for usage instructions.")
        sys.exit(0)
    args = parsers.parseArguments(sys.argv[1:])
    args = parsers.Sanitycheck(args)
    if args.version:
        print(version_text)
        return
   
    if args.create_custom or args.create_protlib or args.create_taulib or args.create_torsion or args.create_config: 
        if args.create_custom:
            generateCustomTemplate(args, filename = 'filter_out.txt')
        if args.create_protlib:
            generateCustomTemplate(args, filename = 'ionizations_v3.txt')
        if args.create_taulib:
            generateCustomTemplate(args, filename = 'tautomers_v3.txt')
        if args.create_torsion:
            generateCustomTemplate(args, filename = 'custom_torsion_templates.txt')
        if args.create_config:
            generateCustomTemplate(args, filename = 'argument_config.yaml')
        print("MolSanitizer templates have been generated. The program exits normally.")
        
    else:
        if args.input_files is not None and args.smiles is None:
            input_path = pathlib.Path(args.input_files[0])
            if args.prefix is not None: log_file = f'{args.prefix}.log' 
            else: log_file = pathlib.Path(f"{input_path.stem}.log")
        else: log_file = 'msani.log'
        loggers.setup_logger(log_file)
        original_command = ' '.join(shlex.quote(arg) for arg in sys.argv)
        logger.info(f"#######  STARTING MOLSANITIZER {__version__}  #######")
        logger.info(f"{original_command}")    
        loggers.arguments(args)
        clean_data(args)
        logger.info(f"***********  MOLSANITIZER FINISHED  ***************")

if __name__=="__main__":
    main()