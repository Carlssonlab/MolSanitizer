from rdkit import Chem

from rdkit.Chem import rdMolTransforms, AllChem
from rdkit.Chem.rdchem import Mol, Conformer
from rdkit.Chem.rdchem import HybridizationType
import os, shutil
import itertools
from collections import defaultdict
import numpy as np
import pandas as pd
import yaml, subprocess


six_membered_aliphatic_substructure = Chem.MolFromSmarts("[A;!$(N-*=*)]1-[A;!$(N-*=*)]-[A;!$(N-*=*)]-[A;!$(N-*=*)]-[A;!$(N-*=*)]-[A;!$(N-*=*)]-1")
sulfonamide_like_substructure = Chem.MolFromSmarts("[*:1][S;$(S(=*)=*):2]-!@[N&+0;!$([NH2]):3](-[*,#1:4])-[*,#1:5]")
aliphatic_nitrogen_substructure = Chem.MolFromSmarts("[A:1]@[N&+0;!$(N-*=*):2](@[A:3])!@[*,#1:4]")
conjugated_substituted_nitrogen = Chem.MolFromSmarts('[a:1]:[a:2]:[a:3]:[nX3&+0:4]-*')
additional_substituted_nitrogen = Chem.MolFromSmarts('[!#1]-[nX3&+0:1]:[a:2]:[a:3]')
amide_substructure = Chem.MolFromSmarts('[O:1]=[CX3:2]!@[N&+0:3](-[!#1:4])-[#1:5]')
symmetric_ring = Chem.MolFromSmarts('*!@-a1[cH][cH][a][cH][cH]1')

with open(os.path.join(os.path.dirname(__file__), 'msani_configurations.yaml')) as confFile:
    msani_configurations = yaml.full_load(confFile)
CORINA_EXE = msani_configurations['CORINA']

class Mol2Writer:
    def __init__(self, mol=None):
        """
        Initialize the writer with an optional RDKit Mol object.
        """
        self.mol = mol

        # Define SMARTS patterns for functional groups
        # These are examples and may need refinement.
        self.smarts_patterns = {
            'carboxylate': Chem.MolFromSmarts("[C;X3](=O)[O-]"),    # COO- 
            'carboxylic_acid': Chem.MolFromSmarts("[C;X3](=O)[OH1]"),   # COOH
            'amide': Chem.MolFromSmarts("[C;X3](=O)N"),             # CONH (amide)
            'phosphate': Chem.MolFromSmarts("[P](=O)(O)(O)(O)"), # H3PO4 or H2PO4- or HPO4-- or C-PO4*
            'so2': Chem.MolFromSmarts("[S;X4](=O)(=O)"),             # SO2
            'so': Chem.MolFromSmarts("[S;X3;!$(S(=O)=O)]=O")                     # SO (sulfoxide)
        }

    def set_mol(self, mol):
        """
        Set the RDKit Mol to be written out as MOL2.
        """
        self.mol = mol

    def write_mol2(self, filename=None):
        """
        Write the RDKit Mol to MOL2 format.
        If filename is given, write to file. Otherwise, return as string.
        """
        if self.mol is None:
            raise ValueError("No molecule set for writing.")

        # Ensure we have some 3D coordinates. If not, try to embed.
        if self.mol.GetNumConformers() == 0:
            res = AllChem.EmbedMolecule(self.mol, useRandomCoords=True)
            if res == -1:
                # If embedding fails, try 2D coordinates and then 3D again
                AllChem.Compute2DCoords(self.mol)
                res = AllChem.EmbedMolecule(self.mol, useRandomCoords=True)
                if res == -1:
                    # Fall back to zero coordinates if embedding still fails
                    AllChem.Compute2DCoords(self.mol)

        atom_types = self._assign_sybyl_types_to_mol(self.mol)
        mol2_string = self._generate_mol2_string(self.mol, atom_types)

        if filename:
            with open(filename, 'w') as f:
                f.write(mol2_string)
        else:
            return mol2_string

    def _generate_mol2_string(self, mol, atom_types):
        """
        Generate the MOL2 file content as a string with aligned columns.
        """
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "MOL"
        conf = mol.GetConformer()

        lines = []
        lines.append("@<TRIPOS>MOLECULE")
        lines.append(mol_name)
        lines.append(f"{num_atoms} {num_bonds} 0 0 0")
        lines.append("SMALL")
        lines.append("USER_CHARGES")

        # ATOM section header
        lines.append("@<TRIPOS>ATOM")
        atom_fmt = "{:<6d} {:<8s} {:>10.4f} {:>10.4f} {:>10.4f} {:<8s} {:<4s} {:<8s} {:>8.4f}"
        for i, atom in enumerate(mol.GetAtoms(), start=1):
            x, y, z = conf.GetAtomPosition(atom.GetIdx()).x, conf.GetAtomPosition(atom.GetIdx()).y, conf.GetAtomPosition(atom.GetIdx()).z
            at_type = atom_types[atom.GetIdx()]
            atom_name = f"{atom.GetSymbol()}{i}"
            charge = float(atom.GetFormalCharge())
            line = atom_fmt.format(i, atom_name, x, y, z, at_type, "1", "LIG", charge)
            lines.append(line)

        # BOND section header
        lines.append("@<TRIPOS>BOND")
        bond_fmt = "{:<6d} {:<6d} {:<6d} {:<3s}"
        for i, bond in enumerate(mol.GetBonds(), start=1):
            a1 = bond.GetBeginAtomIdx() + 1
            a2 = bond.GetEndAtomIdx() + 1
            bond_order = str(int(bond.GetBondTypeAsDouble()))
            line = bond_fmt.format(i, a1, a2, bond_order)
            lines.append(line)

        return "\n".join(lines) + "\n"

    def _assign_sybyl_types_to_mol(self, mol):
        """Assign SYBYL-like atom types to all atoms in the given RDKit Mol."""
        atom_types = [None]*mol.GetNumAtoms()

        # Pre-assign atom types based on recognized functional groups
        self._pre_assign_atom_types(mol, atom_types)

        # For any atom not assigned yet, call the guesser
        for idx, atom in enumerate(mol.GetAtoms()):
            if atom_types[idx] is None:
                atom_types[idx] = self._guess_sybyl_atom_type(atom, mol)

        return atom_types

    def _pre_assign_atom_types(self, mol, atom_types):
        """Use SMARTS to identify known functional groups and assign their atom types directly."""
        # Carboxylate COO-
        # Pattern: [C;X3](=O)[O-]
        # Match order: C_idx, O(double) idx, O(-) idx
        for match in mol.GetSubstructMatches(self.smarts_patterns['carboxylate']):
            c_idx, o_double_idx, o_minus_idx = match
            atom_types[c_idx] = 'C.2'
            atom_types[o_double_idx] = 'O.co2'
            atom_types[o_minus_idx] = 'O.co2'

        # Carboxylic acid COOH
        # Pattern: [C;X3](=O)OH
        # Match order: C_idx, O(double) idx, O(H) idx (though O does not ensure H, we assume acid form)
        for match in mol.GetSubstructMatches(self.smarts_patterns['carboxylic_acid']):
            c_idx, o_double_idx, o_oh_idx = match
            atom_types[c_idx] = 'C.2'
            atom_types[o_double_idx] = 'O.2'
            # The OH oxygen in a COOH is neutral -> O.2
            atom_types[o_oh_idx] = 'O.2'

        # Amide CONH
        # Pattern: [C;X3](=O)N
        # Match order: C_idx, O_idx, N_idx
        for match in mol.GetSubstructMatches(self.smarts_patterns['amide']):
            c_idx, o_idx, n_idx = match
            atom_types[c_idx] = 'C.2'
            atom_types[o_idx] = 'O.2'
            atom_types[n_idx] = 'N.am'


        # Phosphate
        # Pattern: [P](=O)(O)(O)(O)
        # Match order: P_idx, O(double) idx, O idx, O idx, O idx (all neutral)
        for match in mol.GetSubstructMatches(self.smarts_patterns['phosphate']):
            P_idx, O_double_idx, O1_idx, O2_idx, O3_idx = match
            atom_types[P_idx] = 'P.3'
            atom_types[O_double_idx] = 'O.2'
            # The rest are neutral hydroxyl oxygens => O.2
            atom_types[O1_idx] = 'O.3'
            atom_types[O2_idx] = 'O.3'
            atom_types[O3_idx] = 'O.3'

        # SO2 (Sulfone)
        # Pattern: S(=O)(=O)
        for match in mol.GetSubstructMatches(self.smarts_patterns['so2']):
            # match gives S_idx, O1_idx, O2_idx
            # Assign S.o2 for double SO2
            S_idx, O1_idx, O2_idx = match
            atom_types[S_idx] = 'S.o2'
            atom_types[O1_idx] = 'O.2'
            atom_types[O2_idx] = 'O.2'

        # SO (Sulfoxide)
        # Pattern: S=O
        # Note: This might also match sulfone S(=O)(=O), but we handled sulfone already.
        # If an atom is already assigned from the sulfone pattern, skip reassigning.
        for match in mol.GetSubstructMatches(self.smarts_patterns['so']):
            S_idx, O_idx = match
            # Only assign if not already assigned
            if atom_types[S_idx] is None and atom_types[O_idx] is None:
                atom_types[S_idx] = 'S.o'
                atom_types[O_idx] = 'O.2'

    def _guess_sybyl_atom_type(self, atom, mol):
        """Fallback SYBYL-like atom type assignment if SMARTS didn't assign one."""
        symbol = atom.GetSymbol()
        charge = atom.GetFormalCharge()
        hyb = atom.GetHybridization()
        nbrs = atom.GetNeighbors()
        bonds = [mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()) for nbr in nbrs]


        # General fallback logic
        if symbol == 'P':
            return 'P.3'

        elif symbol in ['F', 'Cl', 'Br', 'I', 'H']:
            return symbol

        elif symbol == 'C':
            if charge == 1:
                return 'C.cat'
            if hyb == HybridizationType.SP3:
                return 'C.3'
            elif hyb == HybridizationType.SP2:
                return 'C.2'
            elif hyb == HybridizationType.SP:
                return 'C.1'
            else:
                return 'C.3'

        elif symbol == 'N':
            if charge == 1 and len(nbrs) == 4:
                return 'N.4'
            if hyb == HybridizationType.SP:
                return 'N.1'
            #elif is_amide_nitrogen(): # This is already handled by the SMARTS
            #    return 'N.am'
            elif hyb == HybridizationType.SP2:
                has_double_bond = any(b.GetBondTypeAsDouble() >= 2.0 for b in bonds)
                if has_double_bond:
                    return 'N.2'
                else:
                    return 'N.pl3'
            else:
                return 'N.3'

        elif symbol == 'O':
            # If we reach here, it's not a recognized functional group O.
            # Default logic:
            if hyb == HybridizationType.SP2:
                return 'O.2'
            else:
                return 'O.3'

        elif symbol == 'S':
            # double_o_count = count_double_bonded_oxygens(atom)
            # if double_o_count == 2:
            #     return 'S.o2'
            # elif double_o_count == 1:
            #     return 'S.o' # This is already handled by the SMARTS
            if hyb == HybridizationType.SP2:
                return 'S.2'
            else:
                return 'S.3'



        return symbol

def embed_smiles_corina(smiles, name, VERBOSE):
    with open('temp.smi', 'w') as f:
        f.write(f"{smiles} {name}")
    subprocess.run([CORINA_EXE, '-i', 't=smiles,scn=1,ncn=2', 
                    '-o', 't=mol2', '-d', 'rc,flapn,de=6,mc=1,wh', 'temp.smi', f'{name}.mol2'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if VERBOSE: print(f"\tGenerated mol2 file for {name} using CORINA")
    mol2block = ''
    with open(f"{name}.mol2", 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            else: mol2block += line
    mol_rdkit = Chem.MolFromMol2Block(mol2block, removeHs=False, sanitize=True)
    net_charge = sum([atom.GetFormalCharge() for atom in mol_rdkit.GetAtoms()])
    rigid_scaffolds = [Chem.Mol(mol_rdkit)] # Replicate the output from embed_rdkit
    return mol_rdkit, net_charge, rigid_scaffolds, list()

def find_symmetric_rings(mol_H: Mol):
    '''
    Find the symmetric rings in the molecule. *!@-a1[cH][cH][a][cH][cH]1
    '''
    return mol_H.GetSubstructMatches(symmetric_ring)

def find_flipped_nitrogen(mol_H: Mol):
    '''
    Find the flippable nitrogen in the molecule. n(R):c:c:c:c.
    '''
    return mol_H.GetSubstructMatches(aliphatic_nitrogen_substructure)

def find_conjugated_substituted_nitrogen1(mol_H: Mol):
    '''
        Find the conjugated substituted nitrogen in the molecule. c:c:n(R):c:c. 

        Only match to the atoms within the same ring as n.
    '''
    ring_info = mol_H.GetRingInfo()
    matches = [tuple(match[:4]) for match in mol_H.GetSubstructMatches(conjugated_substituted_nitrogen)]
    filtered_matches = []
    for match in matches:
        for ring in ring_info.AtomRings():
            if set(match).issubset(ring):
                filtered_matches.append(match)
    return filtered_matches

def find_conjugated_substituted_nitrogen2(mol_H: Mol):
    '''
        Another function to find conjugated substituted nitrogen in the molecule fo fix the dihedral

        Find two *-n:a:a matches for each Ns, then fix them to 180 to make them planar
    '''
    matches = mol_H.GetSubstructMatches(additional_substituted_nitrogen)
    ring_info = mol_H.GetRingInfo()
    filtered_matches = []
    for match in matches:
        for ring in ring_info.AtomRings():
            if set(match[1:]).issubset(ring):
                filtered_matches.append(match)
    return filtered_matches

def find_amide(mol_H: Mol):
    '''
        Find the amide in the molecule. O=C-N(-R)-H
    '''
    return mol_H.GetSubstructMatches(amide_substructure)

def calculate_dihedrals_for_rings(conf: Conformer, ring_atoms):
    """Calculate dihedral angles for all torsions involving four consecutive atoms in the ring."""
    ring_atoms = list(ring_atoms)
    dihedrals = [
        (idx1, idx2, idx3, idx4, rdMolTransforms.GetDihedralDeg(conf, idx1, idx2, idx3, idx4))
        for idx1, idx2, idx3, idx4 in [
            (ring_atoms[i], ring_atoms[(i + 1) % len(ring_atoms)], ring_atoms[(i + 2) % len(ring_atoms)], ring_atoms[(i + 3) % len(ring_atoms)])
            for i in range(len(ring_atoms))
        ]
    ]
    return dihedrals

def normalize_angle(angle):
    """Normalize the angle to the range -180 to 180 degrees."""
    return (angle + 180) % 360 - 180


def ring_conf_clusters(current_conf_ring_descriptors, remaining_confs):
    clusters = []
    for remaining_conf in remaining_confs:
        if not current_conf_ring_descriptors == remaining_conf[2:-1]:
            clusters.append(remaining_conf)
    return clusters

def move_and_rename_mol2_files(name: str, num_rigid_scaffolds: int, numPossibleRigidAlignments: int, VERBOSE = False):
    # Loop over the range of scaffold indices
    for i in range(num_rigid_scaffolds):
        for j in range(numPossibleRigidAlignments):
            # Define the source file path for the current scaffold index
            src_file = f"3d/{name}/{name}_mol{i}_align{j}.mol2"
            
            # Define the new file path and destination
            dest_file = src_file.replace('3d', 'db2')
            
            # Ensure the destination directory exists
            os.makedirs(os.path.dirname(dest_file), exist_ok=True)
            
            # Move and rename the file
            shutil.move(src_file, dest_file)
            if VERBOSE: print(f"\tMoved and renamed: {src_file} -> {dest_file}")

def normalize_scores(peaks):
    """Normalize the scores for each peak in a list of peaks."""
    total_score = sum(peak[-1] for peak in peaks)
    return [(angles[:-1], (peak[-1] / total_score) * 100) for angles, peak in zip(peaks, peaks)]

def generate_combinations(data):
    """Generate all combinations of peaks from different dihedrals and calculate the combined score."""
    normalized_data = []
    combination_indices = []
    
    for smarts, atoms, peaks in data:
        normalized_peaks = normalize_scores(peaks)
        normalized_data.append(normalized_peaks)
        combination_indices.append(range(len(peaks)))

    # Generate all possible combinations of normalized peaks
    all_combinations = itertools.product(*combination_indices)
    
    # Create a dictionary to store the results
    combinations_dict = defaultdict(float)
    
    for combination in all_combinations:
        combined_probability = 1
        for dihedral_idx, peak_idx in enumerate(combination):
            combined_probability *= normalized_data[dihedral_idx][peak_idx][-1]
        combinations_dict[combination] = combined_probability
    
    return combinations_dict

def get_flexible_ring(mol: Mol, ssr: list, planar_lib:list, non_planar_lib:list):
    """Extract the flexible rings from the molecule based on the given rules 
    (planar_lib and non_planar_lib). Only consider the rings from SSR to exclude fused rings.

    Args:
        mol (Mol): Rdkit molecule object
        ssr (list): List containing the indices of the smallest set of smallest rings
        planar_lib and non_planar_lib (list): Lists containing 
            (name, smarts, rdkit mol from smarts, reference set of dihedral)

    Returns:
        planar_rings (set): Set containing the indices of the planar rings
        non_planar_rings (list): List containing the indices of the non-planar rings
    """
    planar_rings = set()
    for rule in planar_lib:
        matches = mol.GetSubstructMatches(rule[2])
        if len(matches) > 0:
            for match in matches:
                planar_rings.add((match))
    # Pre-sort the planar_rings and ssr sets
    sorted_planar_rings = {tuple(sorted(t)) for t in planar_rings}
    sorted_ssr = {tuple(sorted(t)) for t in ssr}         
    non_planar_rings = []
    set_non_planar_rings= set()
    for rule in non_planar_lib:

        matches = mol.GetSubstructMatches(rule[2])
        
        if len(matches) > 0:
            for match in matches:
                sorted_match = tuple(sorted(match))
                if sorted_match in sorted_planar_rings: continue 
                if sorted_match not in sorted_ssr: continue
                if sorted_match in set_non_planar_rings: continue
                set_non_planar_rings.add(sorted_match)
                non_planar_rings.append((rule[0], match, [reference_set for reference_set in rule[3]]))
    return planar_rings, non_planar_rings

def calculate_dihedrals_for_rings(conf: Conformer, ring_atoms):
    """Calculate dihedral angles for all torsions involving four consecutive atoms in the ring."""
    ring_atoms = list(ring_atoms)
    dihedrals = []
    for i in range(len(ring_atoms)):
        dihedrals.append(rdMolTransforms.GetDihedralDeg(conf, ring_atoms[i], ring_atoms[(i + 1) % len(ring_atoms)], ring_atoms[(i + 2) % len(ring_atoms)], ring_atoms[(i + 3) % len(ring_atoms)]))
    return dihedrals

def dfs(mol, current_idx, visited, path, forbidden_idxs):
    visited.add(current_idx)
    path.append(current_idx)
    
    for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        # Check if the atom is not visited, not forbidden, and not a hydrogen atom
        if neighbor_idx not in visited and neighbor_idx not in forbidden_idxs:
            if mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() != 1:  # Avoid hydrogen atoms
                dfs(mol, neighbor_idx, visited, path, forbidden_idxs)

def get_dfs_path(mol, start_idx, forbidden_idxs):
    visited = set()
    path = []
    dfs(mol, start_idx, visited, path, forbidden_idxs)
    return path

def identical_substituents(mol, idx2, idx3, idx4, idx5):
    # Create forbidden indices set
    forbidden_idxs = {idx3, idx2}

    # Get paths using DFS
    path4 = get_dfs_path(mol, idx4, forbidden_idxs)
    path5 = get_dfs_path(mol, idx5, forbidden_idxs)
    # Generate SMILES for the paths
    smiles4 = Chem.MolFragmentToSmiles(mol, atomsToUse=path4, rootedAtAtom=idx4)
    smiles5 = Chem.MolFragmentToSmiles(mol, atomsToUse=path5, rootedAtAtom=idx5)

    # Compare the SMILES strings
    return smiles4 == smiles5
    
def find_sulfonamide_like_scaffolds(mol_H: Mol):
    """Find all Sulfonamide-like scaffolds (S(O2)-N(R1)R2 or (S(O)(N)-N(R1)(R2))."""
    preliminary_sulfonamide = mol_H.GetSubstructMatches(sulfonamide_like_substructure)
    matches_sulfonamide = []
    for (a, b, c, d, e) in preliminary_sulfonamide:
        if identical_substituents(mol_H, b, c, d, e): continue
        else: matches_sulfonamide.append((a, b, c, d, e))
    return matches_sulfonamide

def classify_confs(conf, energy, non_planar_rings, flippable_Ns, sulfo_matches, conf_ring_descriptors_df, tolerance=20):
    temp_dict = {
        'Conformer': conf,
        'Energy': energy
    }
    
    # Process non-planar rings
    for ring_idx, (ring_name, ring_atoms, reference_sets) in enumerate(non_planar_rings):
        dihedral_one_ring = calculate_dihedrals_for_rings(conf, ring_atoms)
        sign_vector = np.sign(dihedral_one_ring)
        
        ring_conf_id = -1
        for idx, reference_set in enumerate(reference_sets):
            if np.any(np.sign(reference_set) * sign_vector < 0):
                continue
            abs_diff = np.abs(np.array(dihedral_one_ring) - np.array(reference_set))
            if np.any(abs_diff > tolerance):
                continue
            ring_conf_id = idx
            break
        
        temp_dict[f"{ring_name}{ring_idx}"] = ring_conf_id

    # Process flippable Nitrogens
    flippable_N_descriptors = [1 if rdMolTransforms.GetDihedralDeg(conf, *flippable_N) > 0 else 0 for flippable_N in flippable_Ns]
    temp_dict['flippable_N_descriptors'] = flippable_N_descriptors if flippable_N_descriptors else [-1]

    # Process sulfo matches
    sulfo_descriptors = tuple([1 if rdMolTransforms.GetDihedralDeg(conf, d, b, c, e) > 0 else 0 for (a, b, c, d, e) in sulfo_matches])
    temp_dict['sulfo_descriptors'] = sulfo_descriptors if sulfo_descriptors else [-1]

    # Create dataframe and append to conf_ring_descriptors_df
    temp_df = pd.DataFrame([temp_dict])
    conf_ring_descriptors_df = pd.concat([conf_ring_descriptors_df, temp_df], ignore_index=True)

    return conf_ring_descriptors_df

def remove_unfavorable_confs(conf_ring_descriptors_df: pd.DataFrame, name: str)-> pd.DataFrame:
    for column in conf_ring_descriptors_df.columns[2:-2]:
        if (conf_ring_descriptors_df[column] == -1).all():
            conf_ring_descriptors_df.drop(columns=[column], inplace=True)
            #print(f"\t While handling {name}, found no good conformation of ring: {column}")

    # Step 2: Drop rows where any value in remaining columns (2 onwards) is -1
    for column in conf_ring_descriptors_df.columns[2:-2]:
        conf_ring_descriptors_df = conf_ring_descriptors_df[conf_ring_descriptors_df[column] != -1]
    return conf_ring_descriptors_df

def get_sdf_mol2_filename(name: str, rigid_scaffold_idx: int, align_copy: int):
    return f"{name}_mol{rigid_scaffold_idx}_align{align_copy}.sdf", f"{name}_mol{rigid_scaffold_idx}_align{align_copy}.mol2"

def is_slurm_job():
    '''Check if SLURM_JOB_ID is present in environment variables'''
    return 'SLURM_JOB_ID' in os.environ

def remove_folders(folders_to_remove: list):
    ''' Remove specified folders and all their contents if they exist'''
    for folder in folders_to_remove:
        if os.path.exists(folder) and os.path.isdir(folder):
            try:
                shutil.rmtree(folder)
            except Exception as e:
                print(f"Error removing folder {folder}: {e}")

