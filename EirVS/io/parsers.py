import argparse

from pathlib import Path
import yaml


import logging
logger = logging.getLogger('eirvs')


with open(Path(__file__).parent.parent / 'eirvs_configurations.yaml') as confFile:
    configurations = yaml.full_load(confFile)
    slurm_account = configurations['SLURM_ACCOUNT']
    time_limit = configurations['TIME_LIMIT']
    lines_per_job = configurations['LINES_PER_JOB']
    max_jobs = configurations['MAX_JOBS']
    timeout = configurations['TIMEOUT']
    embed_method = configurations['EMBED_METHOD']
    corina_exe = configurations['CORINA']
    energy_window = configurations['ENERGY_WINDOW']
    numconfs = configurations['NUMCONFS']
    max_stereoisomers = configurations['MAX_STEREOISOMERS']
    pH = configurations['PH']
    pH_range = configurations['PH_RANGE']

info_batch = f"""EirVS - A package to prepare SMILES databases
        This is a batch version of the EirVS package. 
        It reads a list of input files, splits the files into chunks of "--lines_per_job" and processes them parallelly on the HPC.
        For more information, use eirvs_batch -h

        Default settings (modifiable in eirvs_configurations.yaml):
        Project name: {slurm_account}
        Time limit: {time_limit} hour(s)
        Number of compounds per job: {lines_per_job} 
        Maximum jobs at the same time: {max_jobs} job(s)

        Ex. run
        eirvs_batch -i example.smi -l 50 --db2
        eirvs_batch -i example.smi -l 50 --stereoisomers --protonation --db2 --nocleanup
        """

info_standalone = """EirVS - A package to prepare SMILES databases

        Ex. input file (space or tab-separated file):
            COCCC(=O)Nc1ncc(s1)Br  CP000000418470
            C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
            CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597

        Ex. run
        eirvs -i example.smi --removesalts --pains --unwanted all --stereoisomers --protonation
        eirvs -i example.smi --logp "<=500" --hba "<=10" --hbd "<=5" --mw "<=500" -3d -f pdbqt
        eirvs -i example.smi --pains --unwanted regular optional --stereoisomers --protonation
        eirvs -i example.smi --pains --unwanted all --protonation -p 7 -r 1 --tautomers --stereoisomers -3d -f db2.tgz
        """

class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def _format_action_invocation(self, action):
        """
        Override to customize the argument display in the help message.
        Suppress the metavar formatting like `$short $metavar, $long=$metavar`.
        """
        if not action.option_strings:
            return super()._format_action_invocation(action)

        parts = []
        for option_string in action.option_strings:
            parts.append(option_string)
        return ', '.join(parts)

def parseArguments(args = None, batch_mode = False):
    if batch_mode: info = info_batch
    else: info = info_standalone
    if (set(['--db2', '-db2','--pdbqt','--långben']).intersection(set(args))):
        print('\nThe flags --db2,--pdbqt,--långben are deprecated.\nCheck the help message for more information.\n')
        exit(1)
        
    # Create the argument parser
    parser = argparse.ArgumentParser(description=info,
                                     formatter_class=CustomHelpFormatter,
                                     add_help=False)  # Suppress default -h/--help)
    

    # Group 1: Input and output options
    io_group = parser.add_argument_group("Input and output options")
    io_group.add_argument('--input_files', '-i',  type=str,  default=None, nargs='+', help='Input files containing chemical structures')
    io_group.add_argument('--smiles', '-s', default=None, type=str, nargs='+', help='Input SMILES strings')
    io_group.add_argument('--enamine', '-e', action='store_true', help='Enamine input format (default: False)')
    io_group.add_argument('--prefix', '-pre', default=None, type=str, help='Prefix for the output files. (defalt: input file name).')
    io_group.add_argument('--synthon', '-stn',  action='store_true', help='Synthon mode (Additional metadata about the capping groups required)')

    # Group 2: Filtering options
        
    filter_group = parser.add_argument_group(
    "Filtering options", 
    description="""Supported formats for descriptor-based filters (ha, logp, hba, hbd, mw):
    Range: Specify a range using two values (e.g., "17-25").
    Greater / Less than or equal to: Use >= or <= (e.g., ">=17", "<=25").
    Greater than / Less than: Use > or < (e.g., ">17", "<25").
    Exact match: Match a specific value (e.g., 17).
    For logP, the exact match format applies as 'less than or equal to'."""
    )
    filter_group.add_argument('--removesalts', action='store_true', help='Remove salts from the structures. Small fragments within the same molecule are also removed.')
    filter_group.add_argument('--create_custom', action='store_true', help='Generate a template for customized substructure filtering.')
    filter_group.add_argument('--custom', default=None, type=str, help='Filter out unwanted substructures using a customized list.\nTo generate an example list, use --create_custom.')
    filter_group.add_argument('--unwanted', choices=['all', 'regular', 'special', 'optional'], default=None, nargs='*', help='Filter out unwanted substructures using the default list\n(Options: all, regular, special, optional).')
    filter_group.add_argument('--pains', action='store_true', help='Remove PAINS violations from the structures.')
    filter_group.add_argument('--ha', default=None, type=str, help='Filter by the number of heavy atoms.')
    filter_group.add_argument('--logp', default=None, type=str, help='Filter by the value of cLogP*100 (UCSF format: cLogP 3.5->350).')
    filter_group.add_argument('--hba', default=None, type=str, help='Filter by the number of hydrogen bond acceptors.')
    filter_group.add_argument('--hbd', default=None, type=str, help='Filter by the number of hydrogen bond donors.')
    filter_group.add_argument('--mw', default=None, type=str, help='Filter by  molecular weight.')
    filter_group.add_argument('--chiral', default = None, type = str, help='Filter by the number of UNSPECIFIED chiral centers.')

    # Group 3: SMILES processing options
    smiles_group = parser.add_argument_group("SMILES processing options")
    smiles_group.add_argument('--tautomers', '-tau', action='store_true', help='Tautomers enumeration')
    smiles_group.add_argument('--stereoisomers', '-ste', action='store_true', help='Stereoisomers enumeration (only consider unspecified chiral centers)')
    smiles_group.add_argument('--max_stereoisomers','-max_stereo', type=int, default=max_stereoisomers, help=f'Maximum number of stereoisomers to consider (default: {max_stereoisomers} = 3 stereocenters)')
    smiles_group.add_argument('--protonation', '-prot', action='store_true', help='Apply protonation to the structures')
    smiles_group.add_argument('--pH', '-p', type=int, default=pH, help='pH for the protonation (default: 7)')
    smiles_group.add_argument('--pH_range', '-r', type=int, default=pH_range, help='pH range for the protonation (default: 0)')
    smiles_group.add_argument('--noneutralize',  action='store_false', dest='neutralize', default = True, help='Do not neutralize the molecule before tautomerization')
    smiles_group.add_argument('--notaurdkit', action='store_false', dest='taurdkit', default = True, help='Do not use RDKit to canonicalize the tautomeric form of the input SMILES')
    smiles_group.add_argument('--standardize', '-std', action='store_true', dest='standardize', help='Standardize structures for machine learning using RDKit')

    # Group 4: 3D related options
    gen3d = parser.add_argument_group("Generate 3D conformers options")
    gen3d.add_argument('--gen3d', '-3d',  action='store_true', help='Generate 3D conformers')
    gen3d.add_argument('--format', '-f', choices=['db2', 'db2.tgz', 'pdbqt', 'sdf', 'mol2'], default=['db2.tgz'], nargs='*', help='Output file format. Multiple formats simultaneously supported.\n(Default: db2.tgz - Options: sdf, db2, db2.tgz, mol2, pdbqt.)')
    gen3d.add_argument('--method', '-m', choices=['rdkit', 'obabel', 'corina'], dest='method', default = 'rdkit', help=f'Embedding method (default: {embed_method} - options: rdkit, obabel, corina)')
    gen3d.add_argument('--numconfs', '-nconfs', type=int, default=2000, help='Maximum number of conformers to generate (default: 2000)')
    gen3d.add_argument('--randomSeed', '-rs', type=int, default=42, help='Seed for reproducibility (default: 42)')
    gen3d.add_argument('--timeout', '-to', type=float, default=2, help='Timeout for the initial embedding for each SMILES entry before using OpenBabel\nDefault: 2 minutes')
    gen3d.add_argument('--energywindow', '-w', type=float, default=energy_window, help=f'Energy window for sampling the conformations (default: {energy_window} kcal/mol)')
    gen3d.add_argument('--rigid', type = str, default = None, help='Only align the DB2 on this rigid scaffold in SMARTS format. All rings if not provided.')
    gen3d.add_argument('--nringconfs', '-nr', type=int, default=1,
                                            help='Maximum number of ring conformers to generate (default: 1)')
    gen3d.add_argument('--mode', '-mode', choices=['fixed', 'random', 'ignoretorlib'], default='fixed', help='Mode for generating conformers\nDefault: fixed - Options: fixed, random, ignoretorlib')
    gen3d.add_argument('--tolerance', '-tol', type=float, default=30, help='Minimum angle for differentiating two conformers (default: 30)')
    gen3d.add_argument('--nocleanup', action='store_false', dest='cleanup', default = True, help='Do not clean up the temporary files')

    # Group 5: Miscellaneous
    misc_group = parser.add_argument_group("Miscellaneous")
    misc_group.add_argument("--debug", "-d", action="store_true", help="Enable debugging mode")
    misc_group.add_argument('--lazy', action='store_true', help='Implement all the processing and preparation steps')
    misc_group.add_argument('--numcores', '-j', type=int, default=4, help='Number of cores to use for parallel processing (default: 4)')
    misc_group.add_argument("--help", "-h", action="help", help="Show this help message and exit")
    misc_group.add_argument('--timing', action='store_true', help=argparse.SUPPRESS)
    misc_group.add_argument('--test', action='store_true', help=argparse.SUPPRESS)
    misc_group.add_argument('--version', '-v', action='store_true', help = 'Show the current version of EirVS')

    if batch_mode:
        # Group 6: Batch mode options
        batch_group = parser.add_argument_group("Batch mode options")
        batch_group.add_argument('--projectName', '-A', default=slurm_account, dest='proj_name', type=str, help=f'Project name for the SLURM script (default: {slurm_account})')
        batch_group.add_argument('--lines_per_job', '-l', dest='lines', type=int, default=lines_per_job, help=f'Number of lines to process per job (default: {lines_per_job})')
        batch_group.add_argument('--timelimit', '-tl', type=int, default=time_limit, help=f'Time limit for the SLURM job in hours (default: {time_limit})')
        batch_group.add_argument('--max_jobs','-mj', type=int, default=max_jobs, help=f'Maximum number of jobs to run simultaneously (default: {max_jobs})')

    
    # Parse the arguments
    args = parser.parse_args(args if args is not None else [])
    if args.input_files and args.smiles:
        parser.error('Please provide either input files or SMILES strings, not both.')
    if args.input_files is not None:
        for inFile in args.input_files:
            if not Path(inFile).is_file():
                parser.error(f'The input file: {inFile} does not exist.')
        args.input_files = [Path(inFile).resolve() for inFile in args.input_files]
    if args.method == 'corina':
        if not Path(corina_exe).is_file:
            parser.error('Corina path is not correct or corina not found. Please check the configuration file.')
    if args.custom:
        if not Path(args.custom).is_file():
            parser.error(f'The custom file: {args.custom} does not exist.')
        else:
            args.custom = Path(args.custom).resolve()
    if batch_mode: return args, parser
    else: return args

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
    if args.gen3d:
        # Always enumerate stereoisomers for before generating DB2 and PDBQT files
        # Maximum number of stereoisomers is set to in parser
        args.stereoisomers = True
    
    if (args.pH != 7 or args.pH_range != 0) and not args.protonation:
        print("It seems like you forget the --protonation flag. We turned it on for you.")
        args.protonation = True
    return args

