import argparse

def parseArguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Process some chemical structures.")
    
    # Add the required input files argument
    parser.add_argument('input_files', metavar='input_files', type=str, nargs='+',
                        help='Input files containing chemical structures')

    # Add Boolean options
    parser.add_argument('--protonation', action='store_true',
                        help='Apply protonation to the structures')
    parser.add_argument('--removeSalts', action='store_true',
                        help='Remove salts from the structures')
    parser.add_argument('--neutralize', action='store_true',
                        help='Neutralize the structures')
    parser.add_argument('--tautomers', action='store_true',
                        help='Consider tautomers for the structures')
    
    # Add integer option
    parser.add_argument('--maxTautomers', type=int, default=0,
                        help='Maximum number of tautomers to consider (default: 0)')
    
    # Parse the arguments
    args = parser.parse_args()
    
    return args