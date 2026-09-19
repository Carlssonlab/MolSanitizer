import argparse
import sys

from pathlib import Path
from yaml import safe_load



from msani.config import ConfigError, load_defaults

info_standalone = """MolSanitizer - A package to prepare SMILES databases
"""

epilog ="""  Example input file (space or tab-separated file):
        COCCC(=O)Nc1ncc(s1)Br  CP000000418470
        C1CC(C(=O)NC1)SCCC=CBr  CP000000432409
        CC(C)(C)CNC(=O)c1ccsc1Br  CP000001634597

  Extended SMILES is supported with the -e flag (SMILES and IDs need to be tab-separated):        
        CC[C@H]1[C@H](C(=O)N[C@H](C)CCCC(=O)NOCC(F)(F)F)CCN1C |&1:2,3|	Cmp0001
        CCC(CC(=O)N(CC)CCC(=O)N1CCO[C@H]2COC[C@H]21)C(F)F |&1:17,21|	Cmp0002
        CC(C)CC(CNC(=O)C1CSC1)C(=O)N[C@H]1C[C@@H](O)[C@H](F)C1 |&1:16,18,20|	Cmp0003
            
  Personal defaults:
    Bundled defaults are used when no personal configuration exists.
    Run msani config init to customize pH, energywindow, numConfs,
    CORINA path, project settings, and other defaults.
    Run msani config show to inspect defaults and their sources.
    Run msani config -h to see all configuration commands.

  Example usage:
    msani -i example.smi --removesalts --pains --unwanted all --stereoisomers --protonation
    msani -i example.smi --logp "<=500" --hba "<=10" --hbd "<=5" --mw "<=500" -3d -f pdbqt
    msani -i example.smi --pains --unwanted regular optional --stereoisomers --protonation
    msani -i example.smi --pains --unwanted all -prot -p 7 -tau -st -3d -f db2.tgz
    
"""
class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def _format_action_invocation(self, action):
        """
        Override to customize the argument display in the help message.
        Suppress the metavar formatting like `$short $metavar, $long=$metavar`.
        Also hides any automatically generated '--no-*' flags from BooleanOptionalAction.
        """
        if not action.option_strings:
            return super()._format_action_invocation(action)

        # Hide any '--no-*' option variants
        filtered_opts = [
            opt for opt in action.option_strings if not opt.startswith('--no-')
        ]

        return ', '.join(filtered_opts)

def none_or_str(value):
    if value.lower() in ("none", "null", ""):
        return None
    return value

def parseArguments(args = None, batch_mode = False):
    try:
        configurations, _ = load_defaults()
    except ConfigError as exc:
        raise SystemExit(f'Configuration error: {exc}') from exc
    slurm_account = configurations.get('SLURM_ACCOUNT', None)
    slurm_partition = configurations.get('SLURM_PARTITION', None)
    time_limit = configurations.get('TIME_LIMIT', 96)
    lines_per_job = configurations.get('LINES_PER_JOB', 200)
    max_jobs = configurations.get('MAX_JOBS', 1000)
    whole_node = configurations.get('WHOLE_NODE', False)
    whole_node_cores = configurations.get('WHOLE_NODE_CORES', 72)
    timeout = configurations.get('TIMEOUT', 2)
    embed_method = configurations.get('EMBED_METHOD', 'rdkit')
    corina_exe = configurations.get('CORINA') or 'corina'
    energy_window = configurations.get('ENERGY_WINDOW', 25)
    numconfs = configurations.get('NUMCONFS', 2000)
    max_isomers = configurations.get('MAX_STEREOISOMERS', 8)
    pH = configurations.get('PH', 7)
    pH_range = configurations.get('PH_RANGE', 0)

    info_batch = f"""MolSanitizer - A package to prepare SMILES databases
            This is a batch version of the MolSanitizer package.
            It reads a list of input files, splits the files into chunks of "--lines_per_job" and processes them parallelly on the HPC.
            For more information, use msani_batch -h

            Default settings (modifiable in msani_configurations.yaml):
            Project name: {slurm_account}
            Time limit: {time_limit} hour(s)
            Number of compounds per job: {lines_per_job}
            Maximum jobs at the same time: {max_jobs} job(s)

            """

    # Support the parsing of a YAML configuration file
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('--config', '-c', type=str, help='YAML configuration file')
    config_args, _ = pre_parser.parse_known_args(args)

    # Load YAML defaults if provided
    defaults = {}
    if config_args.config:
        try:
            with open(config_args.config, 'r') as f:
                defaults = safe_load(f) or {}
        except Exception as e:
            print(f"Error reading YAML config: {e}")
            sys.exit(1)

    if batch_mode: info = info_batch
    else: info = info_standalone
    if (set(['--db2', '-db2','--pdbqt','--långben']).intersection(set(args))):
        print('\nThe flags --db2,--pdbqt,--långben are deprecated.\nCheck the help message for more information.\n')
        exit(1)
    
    if '--help_advanced' in args or '-xh' in args: show_advanced_help = True
    else: show_advanced_help = False

    # Create the argument parser
    parser = argparse.ArgumentParser(description=info,
                                     formatter_class=CustomHelpFormatter,
                                     add_help=False,
                                     epilog=epilog)  # Suppress default -h/--help)

    
    # Group 1: Input and output options
    io_group = parser.add_argument_group("Input and output options")
    io_group.add_argument(
        '--input_files', '-i',  
        type=str, 
        default=defaults.get('input_files', None), 
        nargs='+', 
        help='Input files containing chemical structures (plain, gzip, bzip2, or xz)')
    io_group.add_argument(
        '--input_list', '-il', 
        type=str, 
        default=defaults.get('input_list', None), 
        help='Path to a text file containing one or more input file paths (one per line).')
    io_group.add_argument(
        '--smiles', '-s', 
        default=defaults.get('smiles', None), 
        type=str, 
        nargs='+', 
        help='Input SMILES strings')
    io_group.add_argument(
        '--extended', '-e', 
        action='store_true', 
        default=defaults.get('extended', False), 
        help='Extended SMILES reading (tab-separated files supported only).')
    io_group.add_argument(
        '--prefix', '-pre', 
        default=defaults.get('prefix', None), 
        type=str, help='Prefix for the output files. (default: input file name).')
    io_group.add_argument(
        '--synthon', '-stn',  
        action='store_true', 
        default=defaults.get('synthon', False), 
        help='Synthon mode (Additional metadata about the capping groups required)')
    

    # Group 2: Filtering options
        
    filter_group = parser.add_argument_group(
    title = "Filtering options", 
    description = 
    """Supported formats for descriptor-based filters (ha, logp, hba, hbd, mw, tpsa, fsp3, chiral):
    Range: Specify a range using two values (e.g., "17-25").
    Greater / Less than or equal to: Use >= or <= (e.g., ">=17", "<=25").
    Greater than / Less than: Use > or < (e.g., ">17", "<25").
    Exact match: Match a specific value (e.g., 17).
    For logP, the exact match format applies as 'less than or equal to'.
    
    Use --ha, --logp, --hba, --hbd, --mw, --tpsa, --fsp3, --chiral to apply these filters."""
    )
    filter_group.add_argument(
        '--removesalts', 
        action='store_true', 
        default=defaults.get('removesalts', False), 
        help='Remove salts from the structures.\nSmall fragments within the same molecule are also removed.')
    filter_group.add_argument(
        '--create_custom', 
        action='store_true', 
        help='Generate a template for customized substructure filtering.')
    filter_group.add_argument(
        '--custom', 
        default=defaults.get('custom', None), 
        type=str, 
        help='Filter out unwanted substructures using a customized list.\nTo generate an example list, use --create_custom.')
    filter_group.add_argument(
        '--unwanted', 
        choices=['all', 'regular', 'special', 'optional'], 
        default=defaults.get('unwanted', None), 
        nargs='*', 
        help='Filter out unwanted substructures using the default list\n(Options: all, regular, special, optional).')
    filter_group.add_argument(
        '--pains', 
        default=defaults.get('pains', False), 
        action='store_true', 
        help='Remove PAINS violations from the structures.')
    filter_group.add_argument(
        '--ha', 
        default=defaults.get('ha', None), 
        type=str, 
        help='Filter by the number of heavy atoms.' if show_advanced_help else argparse.SUPPRESS)
    filter_group.add_argument(
        '--logp', 
        default=defaults.get('logp', None),  
        type=str, 
        help='Filter by cLogP (for example, 3.5). Legacy UCSF-style values multiplied by 100 (for example, 350) are also accepted.' if show_advanced_help else argparse.SUPPRESS)
    filter_group.add_argument(
        '--hba', 
        default=defaults.get('hba', None),  
        type=str, 
        help='Filter by the number of hydrogen bond acceptors.' if show_advanced_help else argparse.SUPPRESS)
    filter_group.add_argument(
        '--hbd', 
        default=defaults.get('hbd', None),  
        type=str, 
        help='Filter by the number of hydrogen bond donors.' if show_advanced_help else argparse.SUPPRESS)
    filter_group.add_argument(
        '--mw', 
        default=defaults.get('mw', None), 
        type=str, 
        help='Filter by molecular weight.' if show_advanced_help else argparse.SUPPRESS)
    filter_group.add_argument(
        '--tpsa',
        default=defaults.get('tpsa', None),
        type=str,
        help='Filter by topological polar surface area (TPSA).' if show_advanced_help else argparse.SUPPRESS)
    filter_group.add_argument(
        '--fsp3',
        default=defaults.get('fsp3', None),
        type=str,
        help='Filter by the fraction of sp3 carbon atoms (range 0-1).' if show_advanced_help else argparse.SUPPRESS)
    filter_group.add_argument(
        '--chiral', 
        default=defaults.get('chiral', None), 
        type = str, 
        help='Filter by the number of UNSPECIFIED chiral centers.' if show_advanced_help else argparse.SUPPRESS)

    # Group 3: SMILES processing options
    smiles_group = parser.add_argument_group("SMILES processing options")
    smiles_group.add_argument(
        '--neutralize', '-neu',
        action=argparse.BooleanOptionalAction,
        default=defaults.get('neutralize', None),
        help='Neutralize molecules.\nWill be applied after removesalts and before tautomerization/protonation (use --no-neutralize to disable).'
    )
    smiles_group.add_argument(
        '--stereoisomers', '-st',
        action=argparse.BooleanOptionalAction,
        default=defaults.get('stereoisomers', None),
        help='Stereoisomers enumeration.\nWill be applied by default when gen3d is on (use --no-stereoisomers to disable)'
    )
    smiles_group.add_argument(
        '--max_isomers','-ms',
        type=int,
        default=defaults.get('max_isomers', max_isomers),
        help=f'Maximum number of stereoisomers to consider (default: {max_isomers})'
    )
    smiles_group.add_argument(
        '--stereo_timeout', '-sto',
        type=int,
        default=defaults.get('stereo_timeout', 60),
        help=(
            'Per-molecule timeout in seconds for stereoisomer enumeration (default: 60).\n'
            'If the timeout is reached, any stereoisomers found up to that point are kept; '
            'if none are found, the input SMILES is kept unchanged.'
        )
    )
    smiles_group.add_argument(
        '--tautomers', '-tau',
        action='store_true',
        default=defaults.get('tautomers', False),
        help='Tautomers enumeration.'
    )
    smiles_group.add_argument(
        '--extended-tautomers', '-et',
        action='store_true',
        default=defaults.get('extended_tautomers', False),
        help="Extended enumeration of tautomers using additional tautomerization rules \n" \
              "(will also activate '--tautomers'). Enumerates less probable tautomeric forms\n"
              "for specific chemotypes: aromatic nitrogen compounds, amidine-like structures,\n"
              "and vinylogous acids (up to 20 forms per compound total)."
    )
    smiles_group.add_argument(
        '--protonation', '-prot',
        action='store_true',
        default=defaults.get('protonation', False),
        help='(De)protonate the structures'
    )
    smiles_group.add_argument(
        '--pH', '-p',
        type=int,
        default=defaults.get('pH', pH),
        help='pH for the protonation (default: 7)'
    )
    smiles_group.add_argument(
        '--pH_range', '-r',
        type=int,
        default=defaults.get('pH_range', pH_range),
        help='pH range for the protonation (default: 0)'
    )
    smiles_group.add_argument(
        '--notaurdkit',
        action='store_false',
        dest='taurdkit',
        default=True,
        help='Do not use RDKit to canonicalize the tautomeric form of the input SMILES' if show_advanced_help else argparse.SUPPRESS
    )
    smiles_group.add_argument(
        '--standardize', '-std',
        action='store_true',
        default=defaults.get('standardize', False),
        help='Standardize structures for machine learning using RDKit'
    )

    # Group 4: 3D related options
    gen3d = parser.add_argument_group("Generate 3D conformers options")
    gen3d.add_argument(
        '--gen3d', '-3d',
        action='store_true',
        default=defaults.get('gen3d', False),
        help='Generate 3D conformers')
    gen3d.add_argument(
        '--format', '-f',
        choices=['db2', 'db2.tgz', 'pdbqt', 'sdf', 'mol2', 'oeb', 'oeb.lib'],
        default=defaults.get('format', ['db2.tgz']),
        nargs='*',
        help='Output file format. Multiple formats simultaneously supported.\n(Default: db2.tgz - Options: sdf, db2, db2.tgz, mol2, pdbqt, oeb, oeb.lib.)\noeb: one .oeb.gz per molecule; oeb.lib: one .oeb.gz library per input file.\nBoth need the optional OpenEye toolkits and OE_LICENSE.')
    gen3d.add_argument(
        '--method', '-m',
        choices=['rdkit', 'obabel', 'corina'],
        default=defaults.get('method', embed_method),
        help=f'Embedding method (default: {embed_method} - options: rdkit, obabel, corina)')
    gen3d.add_argument(
        '--corinaPath',
        type=str,
        default=defaults.get('corinaPath', corina_exe),
        help=f'Path to the CORINA executable (default: {corina_exe})')
    gen3d.add_argument(
        '--numconfs', '-nconfs',
        type=int,
        default=defaults.get('numconfs', numconfs),
        help=f'Maximum number of conformers to generate (default: {numconfs})')
    gen3d.add_argument(
        '--randomSeed', '-rs',
        type=int,
        default=defaults.get('randomSeed', 42),
        help=f'Seed for reproducibility (default: 42)' if show_advanced_help else argparse.SUPPRESS)
    gen3d.add_argument(
        '--timeout', '-to',
        type=float,
        default=defaults.get('timeout', timeout),
        help=f'Timeout for the initial embedding for each entry before using OpenBabel\nDefault: {timeout} minutes')
    gen3d.add_argument(
        '--energywindow', '-w',
        type=float,
        default=defaults.get('energywindow', energy_window),
        help=f'Energy window for sampling the conformations (default: {energy_window} kcal/mol)')
    gen3d.add_argument(
        '--rigid',
        type = str,
        default=defaults.get('rigid', None),
        help='Only align the DB2 on this rigid scaffold in SMARTS format. All rings if not provided.' if show_advanced_help else argparse.SUPPRESS)
    gen3d.add_argument(
        '--nringconfs', '-nr',
        type=int,
        default=defaults.get('nringconfs', 1),
        help='Maximum number of ring conformers to generate (default: 1)')
    gen3d.add_argument(
        '--mode', '-mode',
        choices=['fixed', 'random', 'ignoretorlib'],
        default=defaults.get('mode', 'fixed'),
        help='Mode for generating conformers\nDefault: fixed - Options: fixed, random, ignoretorlib')
    gen3d.add_argument(
        '--tolerance', '-tol',
        type=float,
        default=defaults.get('tolerance', 30),
        help='Minimum angle for differentiating two conformers (default: 30)' if show_advanced_help else argparse.SUPPRESS)
    gen3d.add_argument(
        '--nocleanup', 
        action='store_false', 
        dest='cleanup', 
        default = True, 
        help='Do not clean up the temporary files' if show_advanced_help else argparse.SUPPRESS)
    gen3d.add_argument(
        '--allowNonring', 
        action='store_true', 
        default=defaults.get('allowNonring', True), 
        help='Allow the full sampling of non-ring compounds (default undersample to 30 confs).')
    gen3d.add_argument(
        '--eps', 
        type=float, 
        default=defaults.get('eps', 4), 
        help='The dielectric constant for electrostatic calculations (default: 4).' if show_advanced_help else argparse.SUPPRESS)
    gen3d.add_argument(
        '--rmsd', '-rmsd', 
        type=float, 
        default=defaults.get('rmsd', 0.5), 
        help='Minimum RMSD between two conformers (default: 0.5 Å).' )
    gen3d.add_argument(
        '--timeout_conf', '-toc',
        type=float, 
        default=defaults.get('timeout_conf', 1), 
        help='Timeout for conformational sampling stage (default: 1 minute).' if show_advanced_help else argparse.SUPPRESS)
    gen3d.add_argument(
        '--clash_scale', '-clash',
        type=float,
        default=defaults.get('clash_scale', 0.6),
        help='Scale applied to the sum of atomic van der Waals radii for clash detection (default: 0.6).' if show_advanced_help else argparse.SUPPRESS)
    gen3d.add_argument(
        '--forcefield', '-ff',
        type=str,
        default=defaults.get('forcefield', 'MMFF94s'),
        help='Force field to use (default: MMFF94s).' if show_advanced_help else argparse.SUPPRESS)

    # Group 5: Miscellaneous
    misc_group = parser.add_argument_group("Miscellaneous")
    misc_group.add_argument(
        '--config', '-c',
        type=str, 
        help='Path to the YAML configuration file')
    misc_group.add_argument(
        '--create_config', 
        action='store_true', 
        help='Create a template for the configuration file')
    misc_group.add_argument(
        "--debug", "-d", 
        action="store_true", 
        help="Enable debugging mode" if show_advanced_help else argparse.SUPPRESS)
    misc_group.add_argument(
        '--lazy', 
        action='store_true', 
        default=defaults.get('lazy', False), 
        help='Implement all the processing and preparation steps')
    misc_group.add_argument(
        '--numcores', '-j', 
        type=int, 
        default=defaults.get('numcores', 4), 
        help='Number of cores to use for parallel processing (default: 4)')
    misc_group.add_argument(
        "--help", "-h", 
        action="help", 
        help="Show this help message and exit")
    misc_group.add_argument(
        '--help_advanced', '-xh', 
        action='help', 
        help='Show advanced help message with additional options')
    misc_group.add_argument(
        '--timing', 
        action='store_true', 
        help=argparse.SUPPRESS)
    misc_group.add_argument(
        '--test', 
        action='store_true', 
        help=argparse.SUPPRESS)
    misc_group.add_argument(
        '--version', '-v', 
        action='store_true', 
        help='Show the current version of MolSanitizer')
    misc_group.add_argument(
        '--create_protlib', 
        action='store_true', 
        help='Create a template for customized protonation scheme' if show_advanced_help else argparse.SUPPRESS)
    misc_group.add_argument(
        '--create_taulib', 
        action='store_true', 
        help='Create a template for customized tautomerization scheme' if show_advanced_help else argparse.SUPPRESS)
    misc_group.add_argument(
        '--create_torsion', 
        action='store_true', 
        help='Create a template for customized torsion definition' if show_advanced_help else argparse.SUPPRESS)
    misc_group.add_argument(
        '--protlib',  
        type=str, 
        default=defaults.get('protlib', None), 
        help='Path to the protonation library file (default: msani/Data/ionizations_v3.txt).' if show_advanced_help else argparse.SUPPRESS)
    misc_group.add_argument(
        '--taulib', 
        type=str, 
        default=defaults.get('taulib', None), 
        help='Path to the tautomer library file (default:  msani/Data/tautomers_v3.txt).' if show_advanced_help else argparse.SUPPRESS)
    misc_group.add_argument(
        '--torsion', '-tor', 
        type=str, 
        default=defaults.get('torsion', None), 
        help='Path to the customized torsion definitions.' if show_advanced_help else argparse.SUPPRESS)

    if batch_mode:
        # Group 6: Batch mode options
        batch_group = parser.add_argument_group("Batch mode options")
        batch_group.add_argument(
            '--projectName', '-A',
            default=defaults.get('projectName', slurm_account),
            dest='proj_name', 
            type=none_or_str, 
            help=f'Project name for the SLURM script (default: {slurm_account}).\nUse "none" or "null" to not use a project name')
        batch_group.add_argument(
            '--partition',
            default=defaults.get('partition', slurm_partition),
            dest='partition', 
            type=none_or_str, 
            help=f'Partition for the SLURM script (default: {slurm_partition})')
        batch_group.add_argument(
            '--lines_per_job', '-l',
            dest='lines', 
            type=int, 
            default=defaults.get('lines_per_job', lines_per_job), 
            help=f'Number of lines to process per job (default: {lines_per_job})')
        batch_group.add_argument(
            '--timelimit', '-tl', 
            type=int, 
            default=defaults.get('timelimit', time_limit), 
            help=f'Time limit for the SLURM job in hours (default: {time_limit})')
        batch_group.add_argument(
            '--max_jobs','-mj', 
            type=int, 
            default=defaults.get('max_jobs', max_jobs), 
            help=f'Maximum number of jobs to run simultaneously (default: {max_jobs})')
        batch_group.add_argument(
            '--whole_node', 
            action='store_true', 
            default=defaults.get('whole_node', whole_node),
            help='Run the job on a whole node (default: False)')
        batch_group.add_argument(
            '--whole_node_cores', 
            type=int, 
            default=defaults.get('whole_node_cores', whole_node_cores), 
            help=f'Number of CPU cores per node (default: {whole_node_cores}). Only to set when WHOLE_NODE is true.')

    
    # Parse the arguments
    args = parser.parse_args(args)
    if args.input_files and args.smiles:
        parser.error('Please provide either input files or SMILES strings, not both.')

    if args.input_list is not None:
        if not Path(args.input_list).is_file():
            parser.error(f'The input file: {inFile} does not exist.')
        with open(Path(args.input_list).resolve(), 'r') as f:
            temp_list = [Path(line.strip()).resolve() for line in f if line.strip()]
        if not temp_list:
            parser.error(f'The input list file: {args.input_list} is empty or contains no valid paths.')
        if args.input_files: args.input_files += temp_list
        else: args.input_files = temp_list
        args.input_list = None
            
    if args.input_files is not None:
        for inFile in args.input_files:
            if not Path(inFile).is_file():
                parser.error(f'The input file: {inFile} does not exist.')
        args.input_files = [Path(inFile).resolve() for inFile in args.input_files]

    if args.method == 'corina':
        import shutil
        args.corinaPath = shutil.which(str(Path(args.corinaPath or 'corina').expanduser()))
        if not args.corinaPath:
            parser.error('Corina path is not correct or corina not found. Please check the configuration file.')

    if args.custom:
        if not Path(args.custom).is_file():
            parser.error(f'The custom file: {args.custom} does not exist.')
        else:
            args.custom = Path(args.custom).resolve()

    if args.protlib:
        if not Path(args.protlib).is_file():
            parser.error(f'The protonation library file: {args.protlib} does not exist.')
        else:
            args.protlib = Path(args.protlib).resolve()

    if args.taulib:
        if not Path(args.taulib).is_file():
            parser.error(f'The tautomerization library file: {args.taulib} does not exist.')
        else:
            args.taulib = Path(args.taulib).resolve()

    if args.torsion:
        if not Path(args.torsion).is_file():
            parser.error(f'The torsion definition file: {args.torsion} does not exist.')
        else:
            args.torsion = Path(args.torsion).resolve() 
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
    if args.extended_tautomers and not args.tautomers:
        print("It seems like you forget the --tautomers flag. We turned it on for you.")
        args.tautomers = True

    if (args.pH != 7 or args.pH_range != 0) and not args.protonation:
        print("It seems like you forget the --protonation flag. We turned it on for you.")
        args.protonation = True

    if args.neutralize is None:
        if args.removesalts or args.tautomers or args.protonation:
            print("\nNeutralization is turned on by default when removesalts, tautomerization, or protonation is requested.\nTo disable, use the --no-neutralize flag.\n")
            args.neutralize = True
        else:
            args.neutralize = False

    if args.gen3d:
        # Always enumerate stereoisomers for before generating DB2 and PDBQT files
        # Maximum number of stereoisomers is set to in parser
        if args.stereoisomers is None:
            print("Stereoisomers enumeration is turned on by default when generating 3D conformers.\nTo disable, use the --no-stereoisomers flag.\n")
            args.stereoisomers = True
    
    return args
