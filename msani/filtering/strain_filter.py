"""https://greglandrum.github.io/rdkit-blog/posts/2023-02-04-working-with-conformers.html"""
# Author: Thua-Phong Lam, Jens Carlsson Lab, Uppsala University, July 2024
# This script is used to filter out conformers that do not satisfy the torsion rules in the torlib (last update 2022).
# This is a part of the MolSanitizer project. But could be used as a standalone script.
# strain_filter.py -i mol2_file.mol2 -tol 1 -p prefix

import os
import argparse

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from msani.conformers.torsions import TorsionLibrary
from msani.conformers.utils import within_tolerance
from msani.io.utils import log_error_mol2

def get_atoms_template(pattern):
    pattern_atoms = pattern.GetAtoms()
    enumerated_indices = [-1,-1,-1,-1,-1,-1]
    for idx, atom in enumerate(pattern_atoms):
        if atom.GetAtomMapNum() != 0:
            enumerated_indices[atom.GetAtomMapNum()] = idx
    enumerated_indices.pop(0)
    if enumerated_indices[-1] == -1: enumerated_indices.pop()
    return enumerated_indices

def get_atoms_mol(matches, template_map):
    filtered_matches = []
    
    for match in matches:
        filtered_match = [match[i] for i in template_map]
        if len(filtered_match) == len(template_map):
            filtered_matches.append(filtered_match)
    return filtered_matches





def extract_peaks(match, current_rot_bond):
    for rule in match:
        if set(current_rot_bond) == set(rule[1][1:3]):
            return(rule[1], [(prefered, tolerance1, tolerance2, weight) for prefered, tolerance1, tolerance2, weight in rule[2]])

def parseArguments(args = None):
    info = """StrainFilter - A filtering tool based on TorLib v3
    This tool is used to filter out conformers that do not satisfy the in-house modified TorLib v3.
    The tool will output the filtered conformers in the mol2 format.
    By default, the tool will output the satisfied molecules to filt_{input_name}.mol2, 
    and the strain molecules are saved in strained_{input_name}.mol2.
    
    Ex. run
    strain -i input1.mol2 input2.mol2 -tol 1 
    strain -i input1.mol2 -tol 2 -p strainfiltered
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description= info, formatter_class=argparse.RawTextHelpFormatter)
    
    # Add the required input files argument
    parser.add_argument('-i', '--input_files', type=str, nargs='+', help='Input files containing chemical structures')
    parser.add_argument('-tol', '--tolerance', default=1, choices=[1, 2], type=int, help='Tolerance of the filters (1: only allow relaxed conformers, 2: also allow tolerable conformers)')
    parser.add_argument('-p', '--prefix', default='filt', type=str, help='Prefix for the output files')
    parser.add_argument('-d', '--debug', action='store_true', help='Print debug information')

    # Parse the arguments
    args = parser.parse_args()
    for inFile in args.input_files:
        if not Path(inFile).is_file():
            parser.error(f'The input file: {inFile} does not exist.')

    return args

def write_mol2_file(comments, mol2_block, file_path):
    with open(file_path, 'a') as file:
        file.write(comments + "\n")
        file.write("\n")
        file.write(mol2_block)


def process_one_mol(current_comments_str, current_mol2_str, prefix, input, tol, Torlib, debug):
    """    
    This function processes one molecule by checking if all the dihedral angles are within the tolerance range.
    Args:
        current_comments_str (str): The comments of the current molecule.
        current_mol2_str (str): The mol2 block of the current molecule.
        prefix (str): The prefix for the output files.
        input (str): The path of the current file.
        tol (int): The tolerance of the filters.
        Torlib (list): The list of torsion rules.
    """
    mol = Chem.MolFromMol2Block(current_mol2_str, sanitize=True, removeHs=False)
    name = current_comments_str.split('\n')[0].split()[-1]
    if debug: print(f'Processing {name}')
    if mol:
        if debug: print(f"SMILES: {Chem.MolToSmiles(mol)}")
        match = Torlib.get_match_dihedral(mol, mode = 'random')
        relaxed_angles = [False for _ in match]
        for rule_id, rule in enumerate(match):
            dihedral = rdMolTransforms.GetDihedralDeg(mol.GetConformer(0), *rule[1])
            for (prefered, tol1, tol2, score) in rule[2]:
                tolerance = tol1 if tol == 1 else tol2
                if within_tolerance(dihedral, prefered, tolerance):
                    relaxed_angles[rule_id] = True
                    break
        if all(relaxed_angles) == True:
            write_mol2_file(current_comments_str, current_mol2_str, f'{prefix}_{input}')
        else:
            write_mol2_file(current_comments_str, current_mol2_str, f'strained_{input}')
    else:
        print(f'Error in converting MOL2 {name} to rdkit mol object, saving to strain_error.log')
        log_error_mol2(current_comments_str, current_mol2_str)

def process_mol2_file(input, tol, pre, Torlib, debug = False):
    chunks = []
    current_comments = []
    current_mol2_block = []
    in_molecule_block = False

    with open(input, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith("##########"):
                if in_molecule_block:
                    # If we were in a molecule block, this means we are starting a new molecule, so save the previous one
                    current_comments_str = "\n".join(current_comments)
                    current_mol2_str = "\n".join(current_mol2_block)+"\n"
                    process_one_mol(current_comments_str, current_mol2_str, pre, input, tol, Torlib, debug)
                    current_comments = []
                    current_mol2_block = []
                    in_molecule_block = False
                current_comments.append(line)

            elif line.startswith("@<TRIPOS>MOLECULE"):
                in_molecule_block = True
                current_mol2_block.append(line)

            elif in_molecule_block:
                current_mol2_block.append(line)

        # Add the last molecule block if there is one
        if current_mol2_block:
            current_comments_str = "\n".join(current_comments)
            current_mol2_str = "\n".join(current_mol2_block)+"\n"
            process_one_mol(current_comments_str, current_mol2_str, pre, input, tol, Torlib, debug)
                   
            
    return chunks

def strain_filter(args, Torlib):
    for input_file in args.input_files:
        if os.path.isfile(args.prefix + "_" + input_file): os.remove(args.prefix + "_" + input_file)
        if os.path.isfile("strain_" + input_file): os.remove("strain_" + input_file)
        if os.path.isfile("strain_error.log"): os.remove("strain_error.log")
        process_mol2_file(input_file, args.tolerance, args.prefix, Torlib, args.debug)

def main():
    Torlib = TorsionLibrary()
    args = parseArguments(os.sys.argv[1:])
    strain_filter(args, Torlib)

if __name__ == "__main__":
    main()