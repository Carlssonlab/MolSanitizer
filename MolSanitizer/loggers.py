import logging
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
    logger.info(f"Protonation: {args.protonation}")
    logger.info(f"Remove Salts: {args.removeSalts}")
    logger.info(f"Neutralize: {args.neutralize}")
    logger.info(f"Tautomers: {args.tautomers}")
    logger.info(f"Max Tautomers: {args.maxTautomers}")