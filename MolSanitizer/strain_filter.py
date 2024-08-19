"""https://greglandrum.github.io/rdkit-blog/posts/2023-02-04-working-with-conformers.html"""
# Author: Thua-Phong Lam, Jens Carlsson Lab, Uppsala University, July 2024
# This script is used to filter out conformers that do not satisfy the torsion rules in the torlib (last update 2022).
# This is a part of MolSanitizer project. But could be used as a standalone script.
# strain_filter.py -i mol2_file.mol2 -tol 1 -p prefix
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from pathlib import Path
import os
from os import sys
import numpy as np
import argparse


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


def parse_torlib(xml_file = Path(__file__).parent / 'Data' / 'modified_tor_lib_2020.xml'):
    """This function parse the torlib by the specific class to general class GG, 
    and return a list of tuples with the following format:
    (smarts, rdkit object of the smarts, 4_to_5_atoms_template, [(prefered, tolerance), ...])
    Returns:
        Torlib: list of tuples with the following format:
        (smarts, rdkit object of the smarts, 4_to_5_atoms_template, [(prefered, tolerance), ...])
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    Torlib = []
    for Class in (root.iter(tag='hierarchyClass')):
        if Class.get("name") != "GG": #Not the general class
            for Rule in Class.iter(tag='torsionRule'):
                if  "N_lp" in Rule.get("smarts"): continue
                else:
                    pattern = Chem.MolFromSmarts(Rule.get("smarts"))
                    Torlib.append((Rule.get("smarts"),
                                (pattern),
                                get_atoms_template(pattern),
                                [(((float(angle.get("value")))), float(angle.get("tolerance1")), float(angle.get("tolerance2")), float(angle.get("score"))) for angle in Rule.iter(tag='angle')]))

    for Rule in root.find("hierarchyClass[@name='GG']").iter("torsionRule"):
        if  "N_lp" in Rule.get("smarts"): 
            continue
        else:
            pattern = Chem.MolFromSmarts(Rule.get("smarts"))
            Torlib.append((Rule.get("smarts"),
                        (pattern),
                        get_atoms_template(pattern),
                        [(((float(angle.get("value")))), float(angle.get("tolerance1")), float(angle.get("tolerance2")), float(angle.get("score"))) for angle in Rule.iter(tag='angle')]))
    return Torlib


def normalize_degree(degree):
    """Normalize the degree to the range [-180, 180]."""
    while degree > 180:
        degree -= 360
    while degree < -180:
        degree += 360
    return degree

def handle_lp_rules(mol, rule, match_5_atoms):
    """From the three neighboring atoms to the N atom in the sulfonamide,
    estimate the position of the lone pair and apply the corresponding rule from the TorLib.

    Args:
        mol (_type_): _description_
        Torlib (_type_): _description_
        match_5_atoms (_type_): _description_
    """

    temp_conf = mol.GetConformer(0)
    
    idx_1, idx_2, idx_3, idx_4, idx_5 = match_5_atoms # atom_1, S, N, atom_4, atom_5
    positions = mol.GetConformer().GetPositions()
    N_pos = positions[idx_3]
    neighbor_positions = positions[[idx_2, idx_4, idx_5]] 
    
    # Normalize the distances from nitrogen to each of the neighbors, 
    # form a plane by the three normalized neighbors
    # and calculate the normal vector to the plane.
    # The normal vector will be used to estimate the position of the lone pair.

    normalized_positions = []
    for pos in neighbor_positions:
        distance = np.linalg.norm(N_pos - pos)
        normalized_pos = N_pos + (pos - N_pos) / distance
        normalized_positions.append(normalized_pos)
    normalized_positions = np.array(normalized_positions)

    # Calculate vectors for the plane
    vec1 = normalized_positions[1] - normalized_positions[0]
    vec2 = normalized_positions[2] - normalized_positions[0]

    # Calculate the normal vector to the plane formed by the three neighbors
    normal_vector = np.cross(vec1, vec2)
    normal_vector /= np.linalg.norm(normal_vector)

    # Calculate the vector from nitrogen to the projection point
    projection_length = np.dot(N_pos - normalized_positions[0], normal_vector)
    projection_point = N_pos - projection_length * normal_vector    
    
    vector_N_to_projection = projection_point - N_pos

    # Predict the lone pair position by extending the vector in the opposite direction
    lone_pair_position = N_pos - vector_N_to_projection * 5

    # Create an editable molecule
    editable_mol = Chem.EditableMol(mol)

    # Add a dummy atom (e.g., using '*' for dummy) to represent the lone pair
    dummy_atom = Chem.Atom('*')
    dummy_idx = editable_mol.AddAtom(dummy_atom)

    # Update the molecule
    temp_mol = editable_mol.GetMol()
    conf = temp_mol.GetConformer()
    conf.SetAtomPosition(dummy_idx, lone_pair_position)
    diff_angle1_lp = round(abs(rdMolTransforms.GetDihedralDeg(conf,*(dummy_idx, idx_2, idx_3, idx_4))), 2)
    diff_angle2_lp = round(abs(rdMolTransforms.GetDihedralDeg(conf,*(dummy_idx, idx_2, idx_3, idx_5))), 2)
    #print(diff_angle1_lp, diff_angle2_lp)
    # What we want:                                                In case we choose the wrong reference for rotating, we have:
    #                                 lp                                    
    #                                 |                                              | (not in 0 to -90 deg)
    #                                 N                                              N - lp                          
    #   (first quarter: 0 to 90)    / | \   (fourth quarter: 0 to -90)              /   

    # Idea: randomly rotate the dihedral 1,2,3,4 to +90 degrees (1st quarter), 
    # if the dihedral 1,2,3,5 is within 0 to -90 degrees (fourth quarter), 1,2,3,4 is the correct reference for rotating
    # else, choose 1,2,3,5
    rdMolTransforms.SetDihedralDeg(temp_conf,*(idx_1, idx_2, idx_3, idx_4),value = 90)
    angle_2 = rdMolTransforms.GetDihedralDeg(temp_conf,*(idx_1, idx_2, idx_3, idx_5))
    #print(angle_2)
    new_angles = []
    if within_tolerance(angle_2, -45, 45):
        #print("first assumption correct")
        for angle in rule[3]:
            new_angle = list(angle)
            new_angle[0] = normalize_degree(new_angle[0] - diff_angle1_lp)
            new_angles.append(new_angle)
        return [rule[0], (idx_1, idx_2, idx_3, idx_4), new_angles]
    else:
        for angle in rule[3]:
            new_angle = list(angle)
            new_angle[0] = normalize_degree(new_angle[0] - diff_angle2_lp)
            new_angles.append(new_angle)
        return [rule[0], (idx_1, idx_2, idx_3, idx_5), new_angles]


def get_match_dihedral(mol, Torlib):
    """This function filters the molecule by the torsion rules in the Torlib.
       

    Args:
        mol (_type_): _description_
        Torlib (_type_): _description_
    """
    match = []
    seen = set()
    for idx,rule in enumerate(Torlib):
        matches = mol.GetSubstructMatches(rule[1])
        if len(matches)>0:
            matches_atoms = get_atoms_mol(matches, rule[2])
            if len(matches_atoms[0]) == 5: #Handle lonepair rules
                #print(matches_atoms[0])
                for (a, b, c, d, e) in matches_atoms: 
                    if ((b,c) not in seen) and ((c,b) not in seen):
                        seen.add((b,c))
                        match.append(handle_lp_rules(mol, rule, (a, b, c, d, e)))
            
            else:
                for (a, b, c, d) in matches_atoms:
                    if ((b,c) not in seen) and ((c,b) not in seen):
                        seen.add((b,c))
                        match.append([rule[0], (a, b, c, d), rule[3].copy()])
    return match


def within_tolerance(angle, center, tolerance):
    """
    Check if an angle falls within a specified tolerance range around a central angle in the [-180, 180] range.

    Args:
    angle (float): The angle to check.
    center (float): The central angle.
    tolerance (float): The tolerance range.

    Returns:
    bool: True if the angle is within the tolerance range, False otherwise.
    """
    # Normalize angle and center to the [-180, 180] range
    angle = (angle + 180) % 360 - 180
    center = (center + 180) % 360 - 180
    
    # Calculate the difference considering wrap-around
    diff = (angle - center + 180) % 360 - 180
    return abs(diff) <= tolerance


def extract_peaks(match, current_rot_bond):
    for rule in match:
        if set(current_rot_bond) == set(rule[1][1:3]):
            return(rule[1], [(prefered, tolerance1, tolerance2, weight) for prefered, tolerance1, tolerance2, weight in rule[2]])

def parseArguments(args = None):
    info = """StrainFilter - A filtering tool based on TorLib v3
    This tool is used to filter out conformers that do not satisfy the in-house modified TorLib v3.
    The tool will output the filtered conformers in the mol2 format.

    Ex. run
    strain_filter.py -i input1.mol2 input2.mol2 -tol 1 -p prefix
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description= info, formatter_class=argparse.RawTextHelpFormatter)
    
    # Add the required input files argument
    parser.add_argument('-i', '--input_files', type=str, nargs='+', help='Input files containing chemical structures')
    parser.add_argument('-tol', '--tolerance', default=1, choices=[1, 2], type=int, help='Tolerance of the filters (1: only allow relaxed conformers, 2: also allow tolerable conformers)')
    parser.add_argument('-p', '--prefix', default='filt', type=str, help='Prefix for the output files')

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

def process_one_mol(current_comments_str, current_mol2_str, prefix, input, tol, Torlib):
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
    if mol:
        match = get_match_dihedral(mol, Torlib)
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


def process_mol2_file(input, tol, pre, Torlib):
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
                    process_one_mol(current_comments_str, current_mol2_str, pre, input, tol, Torlib)
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
            process_one_mol(current_comments_str, current_mol2_str, pre, input, tol, Torlib)
                   
            
    return chunks

def strain_filter(args, Torlib):
    for input_file in args.input_files:
        if os.path.isfile(args.prefix + "_" + input_file): os.remove(args.prefix + "_" + input_file)
        process_mol2_file(input_file, args.tolerance, args.prefix, Torlib)

def main():
    Torlib = parse_torlib()
    args = parseArguments(sys.argv[1:])
    strain_filter(args, Torlib)

if __name__ == "__main__":
    main()