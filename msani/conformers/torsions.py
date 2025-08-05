"""https://greglandrum.github.io/rdkit-blog/posts/2023-02-04-working-with-conformers.html"""


import csv
from os import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from numpy.linalg import norm
from numpy import array, cross, dot

from msani.conformers.utils import normalize_angle, within_tolerance

TORLIB_XML = Path(__file__).parent.parent / 'Data' / 'modified_tor_lib_2020.xml'
SMALLRINGLIB_XML = Path(__file__).parent.parent / 'Data/sr_confs.xml'

class TorsionLibrary:
    """
    A class to represent the torsion library.
    It contains methods to parse the torsion rules from an XML file
    and handle the torsion rules.
    """
    
    def __init__(self, xml_file = TORLIB_XML):
        self.xml_file = xml_file
        self.parse_torlib()
    
    def __repr__(self):
        print("Torsion Library from:", self.xml_file)
        print("Number of rules:", len(self.Torlib_specific + self.Torlib_general))

    def add_custom_rule(self, smarts, angles, weights, debug = False):
        """
        Add a custom torsion rule to the library.
        This method allows users to add custom torsion rules one by one.

        Args:
            smarts (str): The SMARTS pattern for the torsion rule.
            angles (list): A list of prefered angles.
            weights (list): A list of weights for the angles.
        """
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            raise ValueError(f"Invalid SMARTS pattern: {smarts}")
        rule, temp = [], []
        for angle, score in zip(angles, weights):
            temp.append((angle, 0, 0, score))
        rule.append((smarts, pattern, TorsionLibrary.get_atoms_template(pattern), temp))
        found = False
        for idx, existing_rule in enumerate(self.Torlib_specific):
            if existing_rule[0] == smarts:
                if debug: print(f"Rule with SMARTS {smarts} already exists. Updating the rule.")
                self.Torlib_specific[idx] = rule
                found = True
                break
        for idx, existing_rule in enumerate(self.Torlib_general):
            if existing_rule[0] == smarts:
                if debug: print(f"Rule with SMARTS {smarts} already exists. Updating the rule.")
                self.Torlib_general[idx] = rule
                found = True
                break
        if not found:
            if debug: print(f"Adding new rule with SMARTS on top priority {smarts}.")
            # Insert the new rule at the beginning of the list
            self.Torlib_specific[:0] = rule  # Insert at the beginning of the list

    def add_custom_rules_from_file(self, file_path, debug = False):
        """
        Add custom torsion rules from a file.
        The template for the file is available at: Data/custom_torlib_template.txt
        
        Args:
            file_path (str): Path to the file containing custom rules.
        """
        custom_rules = dict()
        with open(file_path, 'r') as f:
            reader = csv.reader(f, skipinitialspace=True)
            for row in reader:
                if row and not row[0].startswith("#"):  # Skip empty rows or comments
                    if len(row) == 2:
                        smarts, angles = row[0], row[1].split()
                        weights = [1] * len(angles)
                    elif len(row) == 3:
                        smarts, angles, weights = row[0], row[1].split(), row[2].split()
                    else:
                        print(f"Invalid format in row: {row}. Skipping.")
                        continue
                    if len(weights) != len(angles):
                        print(f"Mismatch in number of angles and weights in row: {row}. Skipping.")
                        continue
                    pattern = Chem.MolFromSmarts(smarts)
                    temp = []
                    for angle, weight in zip(angles, weights):
                        try:
                            angle = float(angle)
                            weight = float(weight)
                            temp.append((angle, 0, 0, weight))  # Assuming tolerance1 and tolerance2 are 0 for custom rules
                        except ValueError:
                            print(f"Invalid angle or weight in row: {row}. Skipping.")
                            continue
                    custom_rules[smarts] = (smarts,
                                              pattern,
                                              TorsionLibrary.get_atoms_template(pattern),
                                              temp)

        if debug: 
            print('\tAdding custom rules:')
            for rule in custom_rules.values():
                print(f'\t\t{rule[0]}: {rule[3]}')
                
        if custom_rules:
            # Check if the custom rules already exist in the library
            for idx, existing_rule in enumerate(self.Torlib_specific):
                if existing_rule[0] in custom_rules.keys():
                    if debug: print(f"\t\tRule with SMARTS {existing_rule[0]} already exists. Updating the rule.")
                    self.Torlib_specific[idx] = custom_rules[existing_rule[0]]
                    del custom_rules[existing_rule[0]]
                    continue
            for idx, existing_rule in enumerate(self.Torlib_general):
                if existing_rule[0] in custom_rules.keys():
                    if debug: print(f"\t\tRule with SMARTS {existing_rule[0]} already exists. Updating the rule.")
                    self.Torlib_general[idx] = custom_rules[existing_rule[0]]
                    del custom_rules[existing_rule[0]]
                    continue
            if custom_rules:
                if debug: 
                    print(f"\tRules that will be placed on top of the priority:")
                    for rule in custom_rules.values():
                        print(f"\t\t{rule[0]}")
                remaining_rules = [rule for rule in custom_rules.values()]
                self.Torlib_specific[:0] = remaining_rules  # Insert at the beginning of the list
        else:
            print("No valid custom rules found in the file.")

    @staticmethod    
    def get_atoms_template(pattern):
            pattern_atoms = pattern.GetAtoms()
            enumerated_indices = [-1,-1,-1,-1,-1,-1]
            for idx, atom in enumerate(pattern_atoms):
                if atom.GetAtomMapNum() != 0:
                    enumerated_indices[atom.GetAtomMapNum()] = idx
            enumerated_indices.pop(0)
            if enumerated_indices[-1] == -1: enumerated_indices.pop()
            return enumerated_indices
    
    def parse_torlib(self):
        """
        This function parse the torlib by the specific class to general class GG, 
        and return a list of tuples with the following format:
        (smarts, rdkit object of the smarts, 4_to_5_atoms_template, [(prefered, tolerance), ...])
        """
        
        
        tree = ET.parse(self.xml_file)
        root = tree.getroot()
        self.Torlib_specific = []
        for Class in (root.iter(tag='hierarchyClass')):
            if Class.get("name") != "GG": #Not the general class
                for Rule in Class.iter(tag='torsionRule'):
                    smarts = Rule.get("smarts")
                    if  "N_lp" in smarts: continue
                    else:
                        pattern = Chem.MolFromSmarts(smarts)
                        self.Torlib_specific.append(
                            (smarts,
                            (pattern),
                            TorsionLibrary.get_atoms_template(pattern),
                            [(((float(angle.get("value")))), float(angle.get("tolerance1")), float(angle.get("tolerance2")), round(float(angle.get("score"))+0.05, 2)) for angle in Rule.iter(tag='angle')])
                            )

        self.Torlib_general = []
        for Rule in root.find("hierarchyClass[@name='GG']").iter("torsionRule"):
            smarts = Rule.get("smarts")
            if  "N_lp" in smarts: 
                continue
            else:
                pattern = Chem.MolFromSmarts(smarts)
                if smarts == "[*:1]~[CX4:2]!@[OX2:3]~[*:4]" or\
                   smarts == "[*:1]~[OX2:2]!@[P:3]~[*:4]" or\
                   smarts == "[*:1]~[CX4:2]!@[SX2:3]~[*:4]":
                    self.Torlib_general.append(
                        (smarts,
                        (pattern),
                        TorsionLibrary.get_atoms_template(pattern),           # Special treatment for aliphatic hydroxyls and phosphates
                        [(((float(angle.get("value")))), float(0), float(0), round(float(angle.get("score"))+0.05, 2)) for angle in Rule.iter(tag='angle')])
                        )
                else:
                    self.Torlib_general.append(
                        (smarts,
                        (pattern),
                        TorsionLibrary.get_atoms_template(pattern),     
                        [(((float(angle.get("value")))), float(angle.get("tolerance1")), float(angle.get("tolerance2")), round(float(angle.get("score"))+0.05, 2)) for angle in Rule.iter(tag='angle')])
                        )
                

    def get_match_dihedral(self, mol, mode = 'fixed'):
        """This function filters the molecule by the torsion rules in the Torlib.
        

        Args:
            mol (Chem.Mol): The RDKit molecule object to be filtered.
            mode (str): The mode of filtering. Can be 'fixed' or 'random
        """

        def get_atoms_mol(matches, template_map):
            '''
            Rematch the order from the template map to define the torsion correctly.
            Several cases in Torlib was defined as [*:1]~[*:2]!@[*:4]~[*:3], so need to be corrected'''
            filtered_matches = []
            
            for match in matches:
                filtered_match = [match[i] for i in template_map]
                if len(filtered_match) == len(template_map):
                    filtered_matches.append(filtered_match)
            return filtered_matches

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
                distance = norm(N_pos - pos)
                normalized_pos = N_pos + (pos - N_pos) / distance
                normalized_positions.append(normalized_pos)
            normalized_positions = array(normalized_positions)

            # Calculate vectors for the plane
            vec1 = normalized_positions[1] - normalized_positions[0]
            vec2 = normalized_positions[2] - normalized_positions[0]

            # Calculate the normal vector to the plane formed by the three neighbors
            normal_vector = cross(vec1, vec2)
            normal_vector /= norm(normal_vector)

            # Calculate the vector from nitrogen to the projection point
            projection_length = dot(N_pos - normalized_positions[0], normal_vector)
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
                    new_angle[0] = round(normalize_angle(new_angle[0] - diff_angle1_lp), 1)
                    new_angles.append(new_angle)
                return [rule[0], (idx_1, idx_2, idx_3, idx_4), new_angles]
            else:
                for angle in rule[3]:
                    new_angle = list(angle)
                    new_angle[0] = round(normalize_angle(new_angle[0] - diff_angle2_lp), 1)
                    new_angles.append(new_angle)
                return [rule[0], (idx_1, idx_2, idx_3, idx_5), new_angles]
            
        match = []
        seen = set()
        for rule in self.Torlib_specific:
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
        
        for rule in self.Torlib_general:
            matches = mol.GetSubstructMatches(rule[1])
            if len(matches)>0:
                matches_atoms = get_atoms_mol(matches, rule[2])
                for (a, b, c, d) in matches_atoms:
                    if ((b,c) not in seen) and ((c,b) not in seen):
                        seen.add((b,c))
                        temp_rule = rule[3].copy()
                        if mode == 'fixed':
                            # Downscale the tolerance 2 to tolerance 1 for the general rules
                            temp_rule = [(peak[0], peak[1], peak[1], peak[3]) for peak in rule[3]] 
                        match.append([rule[0], (a, b, c, d), temp_rule])
        return match

class SmallRingLibrary:
    """
    A class to represent the small ring library.
    It contains methods to parse the small ring conformer library from an XML file.
    """
    
    def __init__(self, xml_file = SMALLRINGLIB_XML):
        self.xml_file = xml_file
        self.parse_sr_confs_library()

    def parse_sr_confs_library(self):
        """
        Parse the XML file containing the SR conformer library and extract the data.
        """

        def parse_dihedral_set(dihedral_str):
        # Remove leading/trailing quotes and split by commas
            dihedral_list = dihedral_str.replace("'", "").split(', ')
            return [int(x) for x in dihedral_list]

        tree = ET.parse(self.xml_file)
        root = tree.getroot()
        self.planar = []
        self.non_planar = []
        for ring in root.findall('ring'):
            name = ring.get('name')  # Extract the ring name
            smarts = ring.find('smarts').get('smarts')  # Extract the smarts string
            mol = Chem.MolFromSmarts(smarts)  # Create RDKit molecule from smarts
            
            dihedral_sets = []
            for dihedral_set in ring.find('dihedral_sets').findall('set'):
                dihedral_str = dihedral_set.get('dihedral')  # Extract dihedral string
                dihedral_list = parse_dihedral_set(dihedral_str)  # Parse dihedral to list of ints
                #value = int(dihedral_set.get('value'))  # Extract the value
                dihedral_sets.append(tuple(dihedral_list))
            
            # Append the tuple of extracted data
            if 'planar' in name: self.planar.append((name, smarts, mol, dihedral_sets))
            else: self.non_planar.append((name, smarts, mol, dihedral_sets))

    def get_planar(self):
        return self.planar
    
    def get_non_planar(self):
        return self.non_planar
    