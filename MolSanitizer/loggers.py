import logging
from pathlib import Path
from .smiles_sanitizer import loadSMARTSdata
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
    logger.info(f"Input files: {args.input_files}")

    if args.enamine: 
        logger.info(f"Enamine format: {args.enamine}")

    if args.removesalts:
        logger.info(f"Remove salts and retain largest fragments: {args.removesalts}")

    if args.protonation or args.tautomers:
        logger.info(f"Neutralize before tautomerization and protonation: {args.neutralize}")

    if args.protonation:
        logger.info(f"Protonation: {args.protonation}")
        logger.info(f"pH: {args.pH}")
        logger.info(f"pH range: {args.pH_range}")


    if args.tautomers:
        logger.info(f"Tautomers enumeration: {args.tautomers}")
        smartsFile = Path(__file__).parent / 'Data' / 'tautomers_v2.txt'
        temp_df = loadSMARTSdata(smartsFile.resolve())
        logger.info(f'Parsed {len(temp_df)} tautomerization rules from: {smartsFile}')

    if args.stereoisomers:
        logger.info(f"Stereoisomers enumeration: {args.stereoisomers}")
        logger.info(f"Max stereoisomers: {args.max_stereoisomers}")

    if args.pains:
        logger.info(f"PAINS filter: {args.pains}")

    if args.unwanted: 
        logger.info(f"Unwanted filter: {args.unwanted}")
        smartsFile = Path(__file__).parent / 'Data' / 'filter_out.csv'
        temp_df = loadSMARTSdata(smartsFile.resolve(), args.unwanted)
        logger.info(f'Parsed {len(temp_df)} substructures from: {smartsFile}')
    
    if args.ha: logger.info(f"HA filter: {args.ha}")
    
    if args.logp: logger.info(f"LogP filter: {args.logp}")

    if args.custom is not None: 
        logger.info(f"Customized filter: {args.custom}")
        temp_df = loadSMARTSdata(args.custom)
        logger.info(f'Parsed {len(temp_df)} substructures from: {args.custom}')

    if args.db2:
        logger.info(f"Generate DB2 files for DOCK3.8: {args.db2}")
        logger.info(f"Number of conformers: {args.numconfs}")
        logger.info(f"Cleanup: {args.cleanup}")
        logger.info(f"Random seed: {args.randomSeed}")
        logger.info(f"Energy window: {args.energywindow}")
        logger.info(f"Use CORINA for initial embedding: {args.corina}")
        if not(args.corina): logger.info(f"Timelimit for initial embedding using RDKit: {args.timeout}")
        if args.enrichment: logger.info(f"Enrichment mode: {args.enrichment}")

