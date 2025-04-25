from rdkit import Chem
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem import AllChem

class Mol2Writer:
    ''' 
    A class to convert the RDKit Mol Object to Mol2 format. 
    Only the first conformer of the Mol object is written out.
    '''
    def __init__(self, mol = None, mol2_template:str = None):
        """
        Initialize the writer with an optional RDKit Mol object.
        Args:
            mol (rdkit.Chem.Mol): An RDKit Mol object with 3D coordinates.
            mol2_template (str): A string containing the MOL2 template.
        """
        self.mol = mol

        
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
        if mol2_template is not None:
            self.atom_types, self.bond = self._read_mol2_template(mol2_template)
        else:
            Chem.Kekulize(self.mol, clearAromaticFlags=True)  
            # Define SMARTS patterns for functional groups
            # These are examples and may need refinement.
            # Still not supporting guanidine,...
            self.smarts_patterns = {
                'carboxylate': Chem.MolFromSmarts("[C;X3](=O)[O-]"),    # COO- 
                'carboxylic_acid': Chem.MolFromSmarts("[C;X3](=O)[OH1]"),   # COOH
                'amide': Chem.MolFromSmarts("[C;X3](=O)N"),             # CONH (amide)
                'phosphate': Chem.MolFromSmarts("[P](=O)(O)(O)(O)"), # H3PO4 or H2PO4- or HPO4-- or C-PO4*
                'so2': Chem.MolFromSmarts("[S;X4](=O)(=O)"),             # SO2
                'so': Chem.MolFromSmarts("[S;X3;!$(S(=O)=O)]=O"),   # SO (sulfoxide)
                'aro6': Chem.MolFromSmarts("A1=A-A=A-A=A1"), # Aromatic 6-membered ring
            }
            self.amide_bonds = set()
            self.atom_types = self._assign_sybyl_types_to_mol(self.mol)
            self.bond = self._bond_section(self.mol)

    def write_mol2(self, filename=None):
        """
        Write the RDKit Mol to MOL2 format.
        If filename is given, write to file. Otherwise, return as string.
        """
        if self.mol is None:
            raise ValueError("No molecule set for writing.")

        # Initialize a string to store the MOL2 content if no filename is provided
        mol2_string = ""
        if filename:
            with open(filename, 'w') as f:        
                for confId in range(self.mol.GetNumConformers()):
                    # Generate the MOL2 string for the current conformer
                    conf_mol2_string = self._generate_mol2_string(self.mol, confId=confId)
                    # Write the MOL2 string to a file or return it
                    f.write(conf_mol2_string)
        else:
            for confId in range(self.mol.GetNumConformers()):
                conf_mol2_string = self._generate_mol2_string(self.mol, confId=confId)
                mol2_string += conf_mol2_string
            return mol2_string
        
    def _read_mol2_template(self, mol2_template):
        """
        Read a MOL2 template file and extract atom types and bonds.
        """
        atom_types = []
        bonds = []
        lines = mol2_template.split("\n")
        # Read the ATOM section
        atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                break
            if atom_section:
                parts = line.split()
                if len(parts) > 5:
                    atom_types.append(parts[5])
        # Read the BOND section
        bond_section = False
        bond_idx = 0
        bond_fmt = "{:>4d} {:>4d} {:>4d} {:<2s}"
        for line in lines:
            if line.startswith("@<TRIPOS>BOND"):
                bond_section = True
                bonds.append(line)
                continue
            if bond_section:
                bond_idx += 1
                parts = line.split()
                if len(parts) > 3:
                    bonds.append(bond_fmt.format(bond_idx, int(parts[1]), int(parts[2]), parts[3]))
        return atom_types, bonds
    
    # ATOM section
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
            self.amide_bonds.add(tuple(sorted((c_idx, n_idx))))


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

        # Aromatic 6-membered ring
        # Pattern: A1=A-A=A-A=A1
        for match in mol.GetSubstructMatches(self.smarts_patterns['aro6']):
            # match gives A1_idx, A2_idx, A3_idx, A4_idx, A5_idx, A6_idx
            # Assign C.ar and N.ar matching the 6-membered aromatic ring.
            for idx in match:
                if mol.GetAtomWithIdx(idx).GetSymbol() == 'O': atom_types[idx] = 'O.3' # pyrrilium
                else: atom_types[idx] = mol.GetAtomWithIdx(idx).GetSymbol() + '.ar'

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
                if atom.GetFormalCharge() == 1:
                    return 'O.3'
                else:
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
    
    # BOND section   
    def _bond_section(self, mol) -> str:
        """
        Generate the BOND section of the MOL2 file.
        """
        lines = []
        lines.append("@<TRIPOS>BOND")
        bond_fmt = "{:>4d} {:>4d} {:>4d} {:>2s}"
        
        for i, bond in enumerate(mol.GetBonds(), start=1):
            a1 = bond.GetBeginAtomIdx() + 1
            a2 = bond.GetEndAtomIdx() + 1
            if all(self.atom_types[i].endswith('.ar') for i in [a1-1, a2-1]):
                bond_order = 'ar'
            elif set([self.atom_types[a1-1], self.atom_types[a2-1]]) == {'C.2', 'O.co2'}:
                bond_order = 'ar'
            elif tuple(sorted((a1-1, a2-1))) in self.amide_bonds:
                bond_order = 'am'
            else:
                bond_order = str(int(bond.GetBondTypeAsDouble()))
            line = bond_fmt.format(i, a1, a2, bond_order)
            lines.append(line)
        
        return lines
    
    def _generate_mol2_string(self, mol, confId=0):
        """
        Generate the MOL2 file content as a string with aligned columns.
        """
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "MOL"
        conf = mol.GetConformer(id=confId)

        lines = []
        lines.append("@<TRIPOS>MOLECULE")
        lines.append(mol_name)
        lines.append(f"{num_atoms}\t{num_bonds}\t0\t0\t0")
        lines.append("SMALL")
        lines.append("NO_CHARGES")

        # ATOM section header
        lines.append("@<TRIPOS>ATOM")
        atom_fmt = "{:>4d} {:<11s} {:>9.4f} {:>10.4f} {:>10.4f} {:<8s}"
        for i, atom in enumerate(mol.GetAtoms(), start=1):
            x, y, z = conf.GetAtomPosition(atom.GetIdx()).x, conf.GetAtomPosition(atom.GetIdx()).y, conf.GetAtomPosition(atom.GetIdx()).z
            at_type = self.atom_types[atom.GetIdx()]
            atom_name = f"{atom.GetSymbol()}{i}"
            #charge = float(atom.GetFormalCharge())
            line = atom_fmt.format(i, atom_name, x, y, z, at_type)#, "1", "LIG", charge)
            lines.append(line)

        # BOND section header
        lines += self.bond

        return "\n".join(lines) + "\n"

    