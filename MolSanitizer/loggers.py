import logging
from pathlib import Path
from .filters import loadSMARTSdata, load_reactions
logger = logging.getLogger('molsani')

def setup_logger(log_file):
    # Create a logger object
    logger = logging.getLogger('molsani')
    logger.setLevel(logging.DEBUG)  # Set the default logging level

    # Create a file handler for logging to a file
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)  # Set the logging level for the file handler

    # Create a stream handler for logging errors to stderr
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.ERROR)  # Set the logging level for the stream handler

    # Create a formatter with the desired format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')

    # Add the formatter to both handlers
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    # Add both handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    # return logger

def arguments(args):
    # Print the parsed arguments 
    logger.info(f"Input files: {args.input_files}")
    logger.info(f"Enamine format: {args.enamine}")
    logger.info(f"Remove Salts: {args.removesalts}")
    logger.info(f"Neutralize before tautomerization and protonation: {args.neutralize}")
    logger.info(f"Tautomers enumeration: {args.tautomers}")
    logger.info(f"PAINS filter: {args.pains}")
    logger.info(f"Unwanted filter: {args.unwanted}")
    logger.info(f"Customized filter: {args.custom}")

    if args.protonation:
        logger.info(f"Protonation: {args.protonation}")
        logger.info(f"pH: {args.pH}")
        logger.info(f"pH range: {args.pH_range}")

    if args.stereoisomers:
        logger.info(f"Stereoisomers: {args.stereoisomers}")
        logger.info(f"Max stereoisomers: {args.max_stereoisomers}")

    if args.db2:
        logger.info(f"Generate DB2 files for DOCK 3.8: {args.db2}")
        logger.info(f"Number of conformers: {args.numconfs}")
        logger.info(f"Cleanup: {args.cleanup}")
        logger.info(f"Random seed: {args.randomSeed}")
        logger.info(f"Energy window: {args.energywindow}")
        logger.info(f"Timelimit for initial embedding using Rdkit: {args.timeout}")
        logger.info(f"Use CORINA for initial embedding: {args.corina}")
        if args.enrichment:
            logger.info(f"Enrichment mode: {args.enrichment}")

    if args.tautomers:
        smartsFile = Path(__file__).parent / 'Data' / 'tautomers.txt'
        temp_df = loadSMARTSdata(smartsFile.resolve())
        logger.info(f'Parsed {len(temp_df)} tautomerization rules from: {smartsFile}')

    if args.unwanted: 
        smartsFile = Path(__file__).parent / 'Data' / 'filter_out.csv'
        temp_df = loadSMARTSdata(smartsFile.resolve(), args.unwanted)
        logger.info(f'Parsed {len(temp_df)} substructures from: {smartsFile}')
    
    if args.protonation:
        smartsFile = Path(__file__).parent / 'Data' / 'ionizations.txt' 
        reactions = load_reactions(smartsFile)
        logger.info(f'Parsed {len(reactions)} ionization reactions from: {smartsFile}')
        
    if args.custom is not None: 
        temp_df = loadSMARTSdata(args.custom)
        logger.info(f'Parsed {len(temp_df)} substructures from: {args.custom}')

