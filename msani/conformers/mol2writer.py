from rdkit import Chem
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem import AllChem

class Mol2Writer:
    ''' 
    A class to convert the RDKit Mol Object to Mol2 format. 
    Only the first conformer of the Mol object is written out.

    Args:
        mol (rdkit.Chem.Mol): RDKit Mol object with 3D coordinates.
        mol2_template (str): A string containing the MOL2 template from CORINA. If provided, atom types and bonds will be taken from this template.
        atom_attributes (bool): If True, include atom attributes (eg. formal charge) in the output.
        If the Mol object does not have 3D coordinates, an attempt will be made to generate them.

    '''
    def __init__(self,
                 mol = None,
                 mol2_template: str = None,
                 atom_attributes:bool = False):
        """
        Initialize the writer with an optional RDKit Mol object.
        Args:
            mol (rdkit.Chem.Mol): An RDKit Mol object with 3D coordinates.
            mol2_template (str): A string containing the MOL2 template.
            atom_attributes (bool): If True, include atom attributes in the output.
        """
        self.mol = mol

        if atom_attributes:
            self.atom_attributes = self._get_atom_attributes()
        else:
            self.atom_attributes = []
        
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
            self.aromatic_ring_bonds = set()
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
        conformer_ids = [conf.GetId() for conf in self.mol.GetConformers()]
        if filename:
            with open(filename, 'w') as f:        
                for confId in conformer_ids:
                    # Generate the MOL2 string for the current conformer
                    conf_mol2_string = self._generate_mol2_string(self.mol, confId=confId)
                    # Write the MOL2 string to a file or return it
                    f.write(conf_mol2_string)
        else:
            for confId in conformer_ids:
                conf_mol2_string = self._generate_mol2_string(self.mol, confId=confId)
                mol2_string += conf_mol2_string
            return mol2_string

    def to_db2_topology(self, name=None, smiles="fake", longname="fake"):
        """Build reusable DB2 topology without conformer-specific data."""
        from msani.db2.molecule import MoleculeData, prepare_molecule_for_db2

        topology = MoleculeData()
        topology.name = name or (
            self.mol.GetProp("_Name") if self.mol.HasProp("_Name") else "fake"
        )
        topology.protein_name = "fake"
        topology.smiles = smiles
        topology.long_name = longname or "fake"

        num_atoms = self.mol.GetNumAtoms()
        if len(self.atom_types) != num_atoms:
            raise ValueError(
                f"MOL2 atom-type count ({len(self.atom_types)}) does not match "
                f"RDKit atom count ({num_atoms})"
            )

        topology.atom_numbers = list(range(1, num_atoms + 1))
        topology.atom_names = [
            f"{atom.GetSymbol()}{atom.GetIdx() + 1}" for atom in self.mol.GetAtoms()
        ]
        topology.atom_types = list(self.atom_types)
        topology.atom_bonds = [[] for _ in range(num_atoms)]

        for bond_line in self.bond:
            if bond_line.startswith("@<TRIPOS>BOND"):
                continue
            tokens = bond_line.split()
            if len(tokens) < 4:
                continue

            bond_num, start, end = map(int, tokens[:3])
            bond_type = tokens[3]
            topology.bond_numbers.append(bond_num)
            topology.bond_starts.append(start)
            topology.bond_ends.append(end)
            topology.bond_types.append(bond_type)
            topology.atom_bonds[start - 1].append((end - 1, bond_type))
            topology.atom_bonds[end - 1].append((start - 1, bond_type))

        return prepare_molecule_for_db2(topology)

    @staticmethod
    def with_db2_conformers(topology, mol):
        """Copy cached topology and attach coordinates from an RDKit molecule."""
        num_atoms = mol.GetNumAtoms()
        if len(topology.atom_numbers) != num_atoms:
            raise ValueError(
                f"Cached topology atom count ({len(topology.atom_numbers)}) does not match "
                f"RDKit atom count ({num_atoms})"
            )

        molecule_data = topology.copy_topology()

        for conf in mol.GetConformers():
            # Transfer the whole conformer once while preserving float tuples.
            molecule_data.conformers.append(
                list(map(tuple, conf.GetPositions().tolist())))

        return molecule_data

    def to_db2_molecule(self, name=None, smiles="fake", longname="fake"):
        """Build complete DB2 molecule data directly from the RDKit molecule."""
        topology = self.to_db2_topology(name=name, smiles=smiles, longname=longname)
        return self.with_db2_conformers(topology, self.mol)
        
    def _read_mol2_template(self, mol2_template):
        """
        Read a MOL2 template file and extract atom types and bonds.
        mol2_template: str - The Mol2 atom types from CORINA.
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
        # Amide CONH
        # Pattern: [C;X3](=O)N
        # Match order: C_idx, O_idx, N_idx
        for match in mol.GetSubstructMatches(self.smarts_patterns['amide']):
            c_idx, o_idx, n_idx = match
            atom_types[c_idx] = 'C.2'
            atom_types[o_idx] = 'O.2'
            atom_types[n_idx] = 'N.am'
            self.amide_bonds.add(tuple(sorted((c_idx, n_idx))))

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
            # Neutral hypervalent and anionic sulfur rings are Kekule rings
            # in CORINA, despite matching the alternating-bond pattern.
            if any(mol.GetAtomWithIdx(idx).GetSymbol() == 'S'
                   and mol.GetAtomWithIdx(idx).GetHybridization() != HybridizationType.SP2
                   for idx in match):
                continue
            for begin, end in zip(match, match[1:] + match[:1]):
                self.aromatic_ring_bonds.add(tuple(sorted((begin, end))))
            # match gives A1_idx, A2_idx, A3_idx, A4_idx, A5_idx, A6_idx
            # Assign C.ar and N.ar matching the 6-membered aromatic ring.
            # Resolve issue 58 related to wrongly perceived S.ar
            for idx in match:
                symbol = mol.GetAtomWithIdx(idx).GetSymbol()
                if symbol in ('O', 'S'):
                    atom_types[idx] = symbol + '.3'
                elif symbol in ('C', 'N'):
                    atom_types[idx] = symbol + '.ar'

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
            # CORINA uses S.o for neutral tetravalent sulfur even without
            # oxygen. Sulfoxides and sulfones were already assigned above.
            if charge == 0 and atom.GetTotalValence() == 4:
                return 'S.o'
            # Charged trivalent sulfur and thiophene-like sulfur use S.3;
            # RDKit's SP2 hybridization alone does not imply SYBYL S.2.
            if charge != 0 and atom.GetTotalValence() == 3:
                return 'S.3'
            if any(b.GetBondTypeAsDouble() == 2.0 for b in bonds):
                return 'S.2'
            else:
                return 'S.3'
        return symbol
    
    # Atom attributes
    def _get_atom_attributes(self):
        """
        Extract atom attributes from the RDKit Mol object.
        """
        atom_attributes = {}
        for atom in self.mol.GetAtoms():
            atom_idx = atom.GetIdx()
            charge = atom.GetFormalCharge()
            if charge != 0:
                atom_attributes[atom_idx] = charge
        return atom_attributes
    
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
            if (tuple(sorted((a1-1, a2-1))) in self.aromatic_ring_bonds
                    or all(self.atom_types[i].endswith('.ar') for i in [a1-1, a2-1])):
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
    
    def _generate_mol2_string(self, 
                              mol, 
                              confId=0):
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
            line = atom_fmt.format(i, atom_name, x, y, z, at_type)
            lines.append(line)
            
        # ATOM attributes if request:
        if self.atom_attributes:
            lines.append("@<TRIPOS>UNITY_ATOM_ATTR")
            for atom_idx, attr in self.atom_attributes.items():
                lines.append(f"{atom_idx+1} 1")
                lines.append(f"charge {attr}")

        # BOND section header
        lines += self.bond

        return "\n".join(lines) + "\n"
