"""https://greglandrum.github.io/rdkit-blog/posts/2023-02-04-working-with-conformers.html"""
# Author: Thua-Phong Lam, Jens Carlsson Lab, Uppsala University, July 2024
# This script is used to filter out conformers that do not satisfy the torsion rules in the torlib (last update 2022).
# This is a part of MolSanitizer project. But could be used as a standalone script.
# strain_filter.py mol2_file.mol2
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from pathlib import Path

def get_4_atoms_template(pattern):
    pattern_atoms = pattern.GetAtoms()
    enumerated_indices = [idx for idx, atom in enumerate(pattern_atoms) if atom.GetAtomMapNum() != 0]
    return enumerated_indices

def get_4_atoms_mol(matches, template_map):
    filtered_matches = []
    
    for match in matches:
        filtered_match = [match[i] for i in template_map]
        if len(filtered_match) == len(template_map):
            filtered_matches.append(filtered_match)
    return filtered_matches

def Mol2MolSupplier (file=None,sanitize=False):
    mols=[]
    with open(file, 'r') as f:
        doc=[line for line in f.readlines()]

    start=[index for (index,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish=[index-1 for (index,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish.append(len(doc))

    interval=list(zip(start,finish[1:]))
    for i in interval:
        block = ",".join(doc[i[0]:i[1]]).replace(',','')
        m=Chem.MolFromMol2Block(block,sanitize=sanitize)
        mols.append(m)
    return(mols)

def parse_torlib():
    """This function parse the torlib by the specific class to general class GG, 
    and return a list of tuples with the following format:
    (smarts, rdkit object of the smarts, 4_atoms_template, [(prefered, tolerance), ...])
    Returns:
        Torlib: list of tuples
    """
    xml_file = Path(__file__).parent / 'Data' / 'tor_lib_2020.xml'
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
                                get_4_atoms_template(pattern),
                                [(((float(angle.get("value")))), float(angle.get("tolerance1"))) for angle in Rule.iter(tag='angle')]))

    for Rule in root.find("hierarchyClass[@name='GG']").iter("torsionRule"):
        if  "N_lp" in Rule.get("smarts"): 
            continue
        else:
            pattern = Chem.MolFromSmarts(Rule.get("smarts"))
            Torlib.append((Rule.get("smarts"),
                        (pattern),
                        get_4_atoms_template(pattern),
                        [(((float(angle.get("value")))), float(angle.get("tolerance1"))) for angle in Rule.iter(tag='angle')]))
    return Torlib

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
            matches_4_atoms = get_4_atoms_mol(matches, rule[2])
            for (a, b, c, d) in matches_4_atoms:
                if ((b,c) not in seen) and ((c,b) not in seen):
                    seen.add((b,c))
                    match.append([rule[0],(a, b, c, d), rule[3].copy()])
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


def check_strain_angle(conf, match, current_rot_bond):
    # Only consider the angles that is involved in the current rotating bond
    for rule in match:
        if set(current_rot_bond) == set(rule[1][1:3]):
            dihedral = rdMolTransforms.GetDihedralDeg(conf, *rule[1])
            #print(current_rot_bond, dihedral)
            for (prefered, tolerance) in rule[2]:
                if within_tolerance(dihedral, prefered, tolerance):    
                    return True 
    return False
                
def check_strain_conformer(conf, match):
    relaxed_conf = False
    relaxed_angles = [False for _ in match]
    for rule_id, rule in enumerate(match):
        dihedral = rdMolTransforms.GetDihedralDeg(conf, *rule[1])
        for (prefered, tolerance) in rule[2]:
            if within_tolerance(dihedral, prefered, tolerance):
                relaxed_angles[rule_id] = True
                break
    relaxed_conf = all(relaxed_angles)
    return relaxed_conf

def filter(mol, match):
    filtered = Chem.Mol(mol)
    filtered.RemoveAllConformers()
    for conf in mol.GetConformers():
        relaxed_conf = check_strain_conformer(conf, match)
        if relaxed_conf: filtered.AddConformer(conf, assignId=True)
    return filtered

def extract_peaks(match, current_rot_bond):
    for rule in match:
        if set(current_rot_bond) == set(rule[1][1:3]):
            if current_rot_bond[1] != rule[1][2]:
                return(rule[1][::-1], [(prefered, tolerance) for prefered, tolerance in rule[2]])
            else: return(rule[1], [(prefered, tolerance) for prefered, tolerance in rule[2]])



if __name__ == "__main__":
    Torlib = parse_torlib()
    script, mol2_file = sys.argv
    mols = Mol2MolSupplier(mol2_file)
    for mol in mols:
        match = get_match_dihedral(mol, Torlib)
        if check_strain_conformer(mol.GetConformer(), match):
            filtered.append(mol)
    #Need to rewrite the mol2 readin and writeout    