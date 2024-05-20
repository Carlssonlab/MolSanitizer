import logging
logger = logging.getLogger('molsani')

import pathlib

from . import parsers
from . import loggers


def main():

    args = parsers.parseArguments()

    input_path = pathlib.Path(args.input_files[0])
    log_file = input_path.with_suffix('.log')
    loggers.setup_logger(log_file)
    
    loggers.arguments(args)





if __name__=="__main__":
    main()