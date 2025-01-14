import argparse

import os
from pathlib import Path
try:
    import yaml
except ImportError:
    print("Please install pyaml using 'pip install pyaml'")
    exit(1)

import logging
logger = logging.getLogger('molsani')

with open(os.path.join(os.path.dirname(__file__), 'msani_configurations.yaml')) as confFile:
    configurations = yaml.full_load(confFile)
    slurm_account = configurations['SLURM_ACCOUNT']
    time_limit = configurations['TIME_LIMIT']
    lines_per_job = configurations['LINES_PER_JOB']
    max_jobs = configurations['MAX_JOBS']
    timeout = configurations['TIMEOUT']
    use_corina = configurations['USE_CORINA']
    corina_exe = configurations['CORINA']
    energy_window = configurations['ENERGY_WINDOW']
    numconfs = configurations['NUMCONFS']
    max_stereoisomers = configurations['MAX_STEREOISOMERS']
    pH = configurations['PH']
    pH_range = configurations['PH_RANGE']

class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        """
        Override to customize the argument display in the help message.
        Suppress the metavar formatting like `$short $metavar, $long=$metavar`.
        """
        if not action.option_strings:
            # Positional argument, return it as-is
            return super()._format_action_invocation(action)
        
        parts = []
        for option_string in action.option_strings:
            # Only display the option string (short/long flag) without metavar
            parts.append(option_string)
        return ', '.join(parts)
    
def parseArguments(args = None):
    info = """MolSanitizer - A package to prepare SMILES databases

    Ex. input file (space or tab-separated file):
        COCCC(=O)Nc1ncc(s1)Br  CP000000418470
        C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
        CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597

    Ex. run
    msani -i example.smi --removesalts --pains --unwanted all --stereoisomers --protonation
    msani -i example.smi --removesalts
    msani -i example.smi --pains --unwanted all --stereoisomers --protonation
    msani -i example.smi --pains --unwanted all --stereoisomers --protonation --db2
    """

    # Create the argument parser
    parser = argparse.ArgumentParser(description=info,
                                     formatter_class=CustomHelpFormatter,
                                     add_help=False)  # Suppress default -h/--help)
    
    # Create a mutually exclusive group
    #group = parser.add_mutually_exclusive_group(required=True)

    # Add the required input files argument

    # Group 1: Input and output options
    io_group = parser.add_argument_group("Input and output options")
    io_group.add_argument('-i', '--input_files', type=str,  default=None, nargs='+', help='Input files containing chemical structures')
    io_group.add_argument('-s', '--smiles', default=None, type=str, nargs='+', help='Input SMILES strings')
    io_group.add_argument('-e', '--enamine', action='store_true', help='Enamine input format (default: False)')
    io_group.add_argument('-pre', '--prefix', default=None, type=str, help='Prefix for the output files. If not provided, the input file name will be used.')
    io_group.add_argument('-enrich', '--enrichment', action='store_true', help='Enrichment mode (do not put in db2.tgz files)')
    io_group.add_argument('-synthon', '--synthon', action='store_true', help='Synthon mode (Additional metadata about the capping groups required)')

    # Group 2: Filtering options
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument('--removesalts', action='store_true', help='Remove salts from the structures. The small fragments within the same molecule are also removed.')
    filter_group.add_argument('--create_custom', action='store_true', help='Generate a template for customized substructure filtering')
    filter_group.add_argument('--custom', default=None, type=str, help='Filter out unwanted substructures using the customized list. To generate an example list, use --create_custom')
    filter_group.add_argument('--unwanted', choices=['all', 'regular', 'special', 'optional'], default=None, nargs='*', help='Filter out unwanted substructures using the default list (choose from all, regular, special, optional)')
    filter_group.add_argument('--pains', action='store_true', help='Remove PAINS violations from the structures')
    filter_group.add_argument('--ha', '-ha', default=None, type=str, help='Retain only compounds with a specified number of heavy atoms. Supported formats: range (e.g., 17-25), less than or equal to (<=25), greater than (>17), etc.')
    filter_group.add_argument('--logp', '-lp', default = None, type = str, help='Retain only the compounds with a specified value of cLogP (UCSF format: cLogP 3.5->350). Supported formats: range (e.g., -350-500), less than or equal to (<=350 or single value 350), greater than (>=200), etc.')

    # Group 3: SMILES processing options
    smiles_group = parser.add_argument_group("SMILES processing options")
    smiles_group.add_argument('--tautomers', action='store_true', help='Tautomers enumeration')
    smiles_group.add_argument('--noneutralize','-noneu',  action='store_false', dest='neutralize', default = True, help='Do not neutralize the molecule before tautomerization')
    smiles_group.add_argument('--notaurdkit', action='store_false', dest='taurdkit', default = True, help='Do not use RDKit to canonicalize the input SMILES')
    smiles_group.add_argument('--stereoisomers', action='store_true', help='Stereoisomers enumeration (only consider unspecified chiral centers)')
    smiles_group.add_argument('--max_stereoisomers','-max_stereo', type=int, default=max_stereoisomers, help=f'Maximum number of stereoisomers to consider (default: {max_stereoisomers} = 3 stereocenters)')
    smiles_group.add_argument('--protonation', action='store_true', help='Apply protonation to the structures')
    smiles_group.add_argument('--pH', '-p', type=int, default=pH, help='pH for the protonation (default: 7)')
    smiles_group.add_argument('--pH_range', '-r', type=int, default=pH_range, help='pH range for the protonation (default: 0)')
    smiles_group.add_argument('--conformal', '-cp', action='store_true', help='Standardize structures for conformal predictors using RDKit')

    # Group 4: DB2 related options
    db2_group = parser.add_argument_group("DB2 related options")
    db2_group.add_argument('--db2', '-db2',  action='store_true', help='Generate conformers and stored in the DB2 format for DOCK 3.8')
    db2_group.add_argument('--corina', '-c', action='store_true', default = use_corina, help=f'Use Corina for 3D structure generation (default: {use_corina})')
    db2_group.add_argument('--långben', '-igtor', action='store_true', dest='ignoretorlib', default = False, help='Ignore the Torsion Library - generate every possible conformer')
    db2_group.add_argument('--numconfs', '-nconfs', type=int, default=2000, help='Maximum number of conformers to generate (default: 2000)')
    db2_group.add_argument('--randomSeed', '-rs', type=int, default=42, help='Seed for reproducibility (default: 42)')
    db2_group.add_argument('--numcores', '-j', type=int, default=4, help='Number of cores to use for parallel processing (default: 4)')
    db2_group.add_argument('--timeout', '-t', type=int, default=2, help='Timeout for the initial embedding for each SMILES entry before using OpenBabel in minutes (default: 2)')
    db2_group.add_argument('--nocleanup', action='store_false', dest='cleanup', default = True, help='Do not clean up the temporary files')
    db2_group.add_argument('--energywindow', '-w', type=float, default=energy_window, help=f'Energy window for sampling the conformations (default: {energy_window} kcal/mol)')

    # Group 5: Miscellaneous
    misc_group = parser.add_argument_group("Miscellaneous")
    misc_group.add_argument("--debug", "-d", action="store_true", help="Enable debugging mode")
    misc_group.add_argument('--lazy', action='store_true', help='Implement all the processing and preparation steps')
    misc_group.add_argument("--help", "-h", action="help", help="Show this help message and exit")
    misc_group.add_argument('--timing', action='store_true', help='Time the process')
    misc_group.add_argument('--test', action='store_true', help=argparse.SUPPRESS)
    misc_group.add_argument('--version', '-v', action='store_true', help = 'Show the current version of MolSanitizer')

    
    # Parse the arguments
    args = parser.parse_args()
    if args.input_files and args.smiles:
        parser.error('Please provide either input files or SMILES strings, not both.')
    if args.input_files is not None:
        for inFile in args.input_files:
            if not Path(inFile).is_file():
                parser.error(f'The input file: {inFile} does not exist.')
    if args.corina:
        if not Path(corina_exe).is_file:
            parser.error('Corina path is not correct or corina not found. Please check the configuration file.')
    return args


def parseArguments_batch(args = None):


    info = f"""MolSanitizer - A package to prepare SMILES databases
    This is a batch version of the MolSanitizer package. 
    It reads a list of input files, splits the files into chunks of "--lines_per_job" and processes them parallelly on the HPC.
    For more information, use msani_batch -h

    Default settings (modifiable in msani_configurations.yaml):
    Project name: {slurm_account}
    Time limit: {time_limit} hour(s)
    Number of compounds per job: {lines_per_job} 
    Maximum jobs at the same time: {max_jobs} job(s)

    Ex. run
    msani_batch -i example.smi -l 50 --db2
    msani_batch -i example.smi -l 50 --stereosiomers --protonation --db2 --nocleanup
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description= info,
                                     formatter_class=CustomHelpFormatter,
                                     add_help=False)  # Suppress default -h/--help)
    
    # Group 1: Input and output options
    io_group = parser.add_argument_group("Input and output options")
    io_group.add_argument('-i', '--input_files', type=str,  default=None, nargs='+', help='Input files containing chemical structures')
    io_group.add_argument('-e', '--enamine', action='store_true', help='Enamine input format (default: False)')
    io_group.add_argument('-pre', '--prefix', default=None, type=str, help='Prefix for the output files. If not provided, the input file name will be used.')
    io_group.add_argument('-enrich', '--enrichment', action='store_true', help='Enrichment mode (do not put in db2.tgz files)')
    io_group.add_argument('-synthon', '--synthon', action='store_true', help='Synthon mode (Additional metadata about the capping groups required)')

    # Group 2: Filtering options
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument('--removesalts', action='store_true', help='Remove salts from the structures. The small fragments within the same molecule are also removed.')
    filter_group.add_argument('--custom', default=None, type=str, help='Filter out unwanted substructures using the customized list. To generate an example list, use --create_custom')
    filter_group.add_argument('--unwanted', choices=['all', 'regular', 'special', 'optional'], default=None, nargs='*', help='Filter out unwanted substructures using the default list (choose from all, regular, special, optional)')
    filter_group.add_argument('--pains', action='store_true', help='Remove PAINS violations from the structures')
    filter_group.add_argument('--ha', '-ha', default=None, type=str, help='Retain only compounds with a specified number of heavy atoms. Supported formats: range (e.g., 17-25), less than or equal to (<=25), greater than (>17), etc.')
    filter_group.add_argument('--logp', '-lp', default = None, type = str, help='Retain only the compounds with a specified value of cLogP (UCSF format: cLogP 3.5->350). Supported formats: range (e.g., -350-500), less than or equal to (<=350 or single value 350), greater than (>=200), etc.')

    # Group 3: SMILES processing options
    smiles_group = parser.add_argument_group("SMILES processing options")
    smiles_group.add_argument('--tautomers', action='store_true', help='Tautomers enumeration')
    smiles_group.add_argument('--noneutralize','-noneu',  action='store_false', dest='neutralize', default = True, help='Do not neutralize the molecule before tautomerization and protonation')
    smiles_group.add_argument('--notaurdkit', action='store_false', dest='taurdkit', default = True, help='Do not use RDKit to canonicalize the input SMILES')
    smiles_group.add_argument('--stereoisomers', action='store_true', help='Stereoisomers enumeration (only consider unspecified chiral centers)')
    smiles_group.add_argument('--max_stereoisomers','-max_stereo', type=int, default=max_stereoisomers, help=f'Maximum number of stereoisomers to consider (default: {max_stereoisomers} = 3 stereocenters)')
    smiles_group.add_argument('--protonation', action='store_true', help='Apply protonation to the structures')
    smiles_group.add_argument('--pH', '-p', type=int, default=pH, help='pH for the protonation (default: 7)')
    smiles_group.add_argument('--pH_range', '-r', type=int, default=pH_range, help='pH range for the protonation (default: 0)')
    smiles_group.add_argument('--conformal', '-cp', action='store_true', help='Standardize structures for conformal predictors using RDKit')

    # Group 4: DB2 related options
    db2_group = parser.add_argument_group("DB2 related options")
    db2_group.add_argument('--db2', '-db2',  action='store_true', help='Generate conformers and stored in the DB2 format for DOCK 3.8')
    db2_group.add_argument('--corina', '-c', action='store_true', default = use_corina, help=f'Use Corina for 3D structure generation (default: {use_corina})')
    db2_group.add_argument('--långben', '-igtor', action='store_true', dest='ignoretorlib', default = False, help='Ignore the Torsion Library - generate every possible conformer')
    db2_group.add_argument('--numconfs', '-nconfs', type=int, default=2000, help='Maximum number of conformers to generate (default: 2000)')
    db2_group.add_argument('--randomSeed', '-rs', type=int, default=42, help='Seed for reproducibility (default: 42)')
    db2_group.add_argument('--numcores', '-j', type=int, default=4, help='Number of cores to use for parallel processing (default: 4)')
    db2_group.add_argument('--timeout', '-t', type=int, default=2, help='Timeout for the initial embedding for each SMILES entry before using OpenBabel in minutes (default: 2)')
    db2_group.add_argument('--nocleanup', action='store_false', dest='cleanup', default = True, help='Do not clean up the temporary files')
    db2_group.add_argument('--energywindow', '-w', type=float, default=energy_window, help=f'Energy window for sampling the conformations (default: {energy_window} kcal/mol)')

    # Group 5: Miscellaneous
    misc_group = parser.add_argument_group("Miscellaneous")
    misc_group.add_argument("--debug", "-d", action="store_true", help="Enable debugging mode")
    misc_group.add_argument('--lazy', action='store_true', help='Implement all the processing and preparation steps (default: False)')
    misc_group.add_argument("--help", "-h", action="help", help="Show this help message and exit")
    misc_group.add_argument('--timing', action='store_true', help='Time the process')

    # Group 6: Batch mode options
    batch_group = parser.add_argument_group("Batch mode options")
    batch_group.add_argument('-n', '--projectName', default=slurm_account, dest='proj_name', type=str, help=f'Project name for the SLURM script (default: {slurm_account})')
    batch_group.add_argument('-l', '--lines_per_job', dest='lines', type=int, default=lines_per_job, help=f'Number of lines to process per job (default: {lines_per_job})')
    batch_group.add_argument('-tl', '--timelimit', type=int, default=time_limit, help=f'Time limit for the SLURM job in hours (default: {time_limit})')
    batch_group.add_argument('-mj', '--max_jobs', type=int, default=max_jobs, help=f'Maximum number of jobs to run simultaneously (default: {max_jobs})')

    # Parse the arguments
    args = parser.parse_args()
    for inFile in args.input_files:
        if not Path(inFile).is_file():
            parser.error(f'The input file: {inFile} does not exist.')
    if args.corina:
        if not Path(corina_exe).is_file:
            parser.error('Corina path is not correct or corina not found. Please check the msani_configurations.yaml file.')
    return args
