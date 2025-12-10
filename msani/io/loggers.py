import logging
from pathlib import Path
import os
logger = logging.getLogger('msani')

def setup_logger(log_file):
    # Create a logger object
    logger = logging.getLogger('msani')
    logger.setLevel(logging.DEBUG)  # Set the default logging level

    # Create a file handler for logging to a file
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)  # Set the logging level for the file handler

    # Create a stream handler for logging errors to stderr
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.ERROR)  # Set the logging level for the stream handler

    # Create a formatter with the desired format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    # Add the formatter to both handlers
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    # Add both handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    # return logger

def arguments(args):
    # Print the parsed arguments 
    job_id = os.getenv('SLURM_JOB_ID')
    if job_id:
        logger.info(f"Job ID: {job_id}")
    if args.input_files:
        logger.info(f"Input files: {[str(input_file) for input_file in args.input_files]}")

    if args.extended: 
        logger.info(f"Parsing the extended SMILES format: {args.extended}")

    if args.removesalts:
        logger.info(f"Remove salts and retain largest fragments: {args.removesalts}")
        if args.debug:
            smartsFile = Path(__file__).parent.parent / 'Data' / 'salt_stripping.txt'
            logger.info(f'Loading salt stripping rules from: {smartsFile.resolve()}')

    if args.protonation or args.tautomers or args.removesalts:
        logger.info(f"Neutralize after removesalts and before tautomerization/protonation: {args.neutralize}")

    if args.protonation:
        logger.info(f"Protonation: {args.protonation}")
        logger.info(f"pH: {args.pH}")
        logger.info(f"pH range: {args.pH_range}")
        if args.protlib:
            logger.info(f"Using protonation library: {args.protlib}")


    if args.tautomers:
        logger.info(f"Tautomers enumeration: {args.tautomers}")
        if args.taulib:
            logger.info(f"Using tautomer library: {args.taulib}")

    if args.stereoisomers:
        logger.info(f"Stereoisomers enumeration: {args.stereoisomers}")
        logger.info(f"Max stereoisomers: {args.max_isomers}")

    if args.pains:
        logger.info(f"PAINS filter: {args.pains}")

    
    if args.ha: logger.info(f"HA filter: {args.ha}")
    if args.logp: logger.info(f"LogP filter: {args.logp}")
    if args.hba: logger.info(f"HBA filter: {args.hba}")
    if args.hbd: logger.info(f"HBD filter: {args.hbd}")
    if args.mw: logger.info(f"MW filter: {args.mw}")
    if args.chiral: logger.info(f"Chiral filter: {args.chiral}")

    if args.gen3d:
        logger.info(f"Generate 3D conformers: {args.gen3d}")
        logger.info(f"Output format: {args.format}")
        logger.info(f"Sampling mode: {args.mode}")
        if args.mode == 'random':
            logger.info(f"Dihedral tolerance: {args.tolerance}")
        logger.info(f"Number of conformers: {args.numconfs}")
        logger.info(f"Cleanup: {args.cleanup}")
        logger.info(f"Random seed: {args.randomSeed}")
        logger.info(f"Energy window: {args.energywindow}")
        logger.info(f"Number of ring conformations: {args.nringconfs}")
        logger.info(f"Embedding method: {args.method}")
        if args.method == 'rdkit': logger.info(f"Timelimit for initial embedding using RDKit: {args.timeout}")
        if args.rigid: logger.info(f"Only align based on: {args.rigid}")

        

