# -*- coding: utf-8 -*-
"""
molvs.tautomer
~~~~~~~~~~~~~~

This module contains tools for enumerating tautomers and determining a canonical tautomer.
This is a part of the molvs package, which could be found here: https://github.com/mcs07/MolVS/blob/master/molvs/tautomer.py
We modified the original code to match with the improvements made by the RDKit team.
molvs is under MIT license.
RDKit is under BSD-3 license.
Date modified May 1st, 2025
"""

import copy
import logging

from rdkit import Chem
from rdkit.Chem import SanitizeFlags
from rdkit.Chem.rdchem import BondDir, BondStereo, BondType, HybridizationType

from .utils import memoized_property, pairwise

log = logging.getLogger('eirvs')


class TautomerTransform(object):
    """Rules to transform one tautomer to another.

    Each TautomerTransform is defined by a SMARTS pattern where the transform involves moving a hydrogen from the first
    atom in the pattern to the last atom in the pattern. By default, alternating single and double bonds along the
    pattern are swapped accordingly to account for the hydrogen movement. If necessary, the transform can instead define
    custom resulting bond orders and also resulting atom charges.
    """

    BONDMAP = {'-': BondType.SINGLE, '=': BondType.DOUBLE, '#': BondType.TRIPLE, ':': BondType.AROMATIC}
    CHARGEMAP = {'+': 1, '0': 0, '-': -1}

    def __init__(self, name, smarts, bonds=(), charges=(), radicals=()):
        """Initialize a TautomerTransform with a name, SMARTS pattern and optional bonds and charges.

        The SMARTS pattern match is applied to a Kekule form of the molecule, so use explicit single and double bonds
        rather than aromatic.

        Specify custom bonds as a string of ``-``, ``=``, ``#``, ``:`` for single, double, triple and aromatic bonds
        respectively. Specify custom charges as ``+``, ``0``, ``-`` for +1, 0 and -1 charges respectively.

        :param string name: A name for this TautomerTransform.
        :param string smarts: SMARTS pattern to match for the transform.
        :param string bonds: Optional specification for the resulting bonds.
        :param string charges: Optional specification for the resulting charges on the atoms.
        """
        self.name = name
        self.tautomer_str = smarts
        self.bonds = [self.BONDMAP[b] for b in bonds]
        self.charges = [self.CHARGEMAP[b] for b in charges]
        # TODO: Raise error (ValueError?) if bonds and charges lists are not the correct length

    @memoized_property
    def tautomer(self):
        return Chem.MolFromSmarts(self.tautomer_str)

    def __repr__(self):
        return 'TautomerTransform({!r}, {!r}, {!r}, {!r})'.format(self.name, self.tautomer_str, self.bonds, self.charges)

    def __str__(self):
        return self.name


class TautomerScore(object):
    """A substructure defined by SMARTS and its score contribution to determine the canonical tautomer."""

    def __init__(self, name, smarts, score):
        """Initialize a TautomerScore with a name, SMARTS pattern and score.

        :param name: A name for this TautomerScore.
        :param smarts: SMARTS pattern to match a substructure.
        :param score: The score to assign for this substructure.
        """
        self.name = name
        self.smarts_str = smarts
        self.score = score

    @memoized_property
    def smarts(self):
        return Chem.MolFromSmarts(self.smarts_str)

    def __repr__(self):
        return 'TautomerScore({!r}, {!r}, {!r})'.format(self.name, self.smarts_str, self.score)

    def __str__(self):
        return self.name


#: The default list of TautomerTransforms.
# Also keep pace with RDKit (updated extended SMARTS for many rules)
TAUTOMER_TRANSFORMS = (
    TautomerTransform('1,3 (thio)keto/enol f', '[CX4!H0R{0-2}]-[C;z{1-2}]=[O,S,Se,Te;X1]'),
    TautomerTransform('1,3 (thio)keto/enol r', '[O,S,Se,Te;X2!H0]-[#6;z{1-2}]=[C,cz{0-1}R{0-1}]'),
    TautomerTransform('1,5 (thio)keto/enol f', '[CX4z0,NX3;!H0]-[C]=[C][Cz1H0]=[O,S,Se,Te;X1]'),
    TautomerTransform('1,5 (thio)keto/enol r', '[O,S,Se,Te;X2!H0]-[Cz1H0]=[C]-[C]=[Cz0,N]'),
    TautomerTransform('aliphatic imine f', '[CX4R{0-2}!H0]-[Cz1]=[NX2]'),
    TautomerTransform('aliphatic imine r', '[NX3!H0]-[C;z{1-2}]=[CX3]'),
    TautomerTransform('special imine f', '[Nz0!H0]-[C]=[Cz0X3R0]'),
    TautomerTransform('special imine r', '[Cz0R0X4!H0]-[c]=[nz0]'),
    TautomerTransform('1,3 aromatic heteroatom H shift f', '[#7+0!H0]-[#6R1]=[O,#7X2+0]'),
    TautomerTransform('1,3 aromatic heteroatom H shift r', '[O,#7+0;!H0]-[#6R1]=[#7+0X2]'),
    TautomerTransform('1,3 heteroatom H shift', '[#7+0,S,O,Se,Te;!H0]-[#7X2,#6,#15]=[#7+0,#16,#8,Se,Te]'),
    TautomerTransform('1,5 aromatic heteroatom H shift', '[#7+0,#16,#8;!H0]-[#6,#7]=[#6]-[#6,#7]=[#7+0,#16,#8;H0]'),
    TautomerTransform('1,5 aromatic heteroatom H shift f', '[#7+0,#16,#8,Se,Te;!H0]-[#6,nX2]=[#6,nX2]-[#6,#7X2]=[#7X2+0,S,O,Se,Te]'),
    TautomerTransform('1,5 aromatic heteroatom H shift r', '[#7+0,S,O,Se,Te;!H0]-[#6,#7X2]=[#6,nX2]-[#6,nX2]=[#7+0,#16,#8,Se,Te]'),
    TautomerTransform('1,7 aromatic heteroatom H shift f', '[#7+0,#8,#16,Se,Te;!H0]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6]-[#6,#7X2]=[#7X2+0,S,O,Se,Te,Cz0X3]'),
    TautomerTransform('1,7 aromatic heteroatom H shift r', '[#7+0,S,O,Se,Te,Cz0X4;!H0]-[#6,#7X2]=[#6]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[NX2,S,O,Se,Te]'),
    TautomerTransform('1,9 aromatic heteroatom H shift f', '[#7+0,O;!H0]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#7+0,O]'),
    TautomerTransform('1,11 aromatic heteroatom H shift f', '[#7+0,O;!H0]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#7X2+0,O]'),
    TautomerTransform('furanone f', '[O,S,N;!H0]-[#6z2r5]=,:[#6X3r5;$([#6]([#6;r5])=,:[#6X3r5])]'),
    TautomerTransform('furanone r', '[#6r5!H0;$([#6]([#6r5])[#6r5])][#6z2r5]=[O,S,N]'),
    TautomerTransform('keten/ynol f', '[C!H0]=[C]=[O,S,Se,Te;X1]', bonds='#-'),
    TautomerTransform('keten/ynol r', '[O,S,Se,Te;!H0X2]-[C]#[C]', bonds='=='),
    TautomerTransform('ionic nitro/aci-nitro f', '[C!H0]-[N+;$([N][O-])]=[O]'),
    TautomerTransform('ionic nitro/aci-nitro r', '[O!H0]-[N+;$([N][O-])]=[C]'),
    TautomerTransform('oxim/nitroso f', '[O!H0]-[Nz1]=[C]'),
    TautomerTransform('oxim/nitroso r', '[C!H0]-[Nz1]=[O]'),
    TautomerTransform('oxim/nitroso via phenol f', '[O!H0]-[N]=[C]-[C]=[C]-[C]=[OH0]'),
    TautomerTransform('oxim/nitroso via phenol r', '[O!H0]-[c]=,:[c][c]=,:[c]-[N]=[OH0]'),
    TautomerTransform('cyano/iso-cyanic acid f', '[O!H0]-[C]#[N]', bonds='=='),
    TautomerTransform('cyano/iso-cyanic acid r', '[N!H0]=[C]=[O]', bonds='#-'),
    TautomerTransform('formamidinesulfinic acid f', '[O,N;!H0]-[C]=[S,Se,Te;v6]=[O]', bonds='=--'),  # ACTIVATED
    TautomerTransform('formamidinesulfinic acid r', '[O!H0]-[S,Se,Te;v4]-[C]=[O,N]', bonds='==-'),
    TautomerTransform('isocyanide f', '[C-0!H0]#[N+0]', bonds='#', charges='-+'),
    TautomerTransform('isocyanide r', '[N+!H0]#[C-]', bonds='#', charges='-+'),
    TautomerTransform('phosphonic acid f', '[OH]-[PH0]', bonds='='),
    TautomerTransform('phosphonic acid r', '[PH]=[O]', bonds='-'),
)

#: The default list of TautomerScores.
# Original MolVS has 10, RDKit has 12
TAUTOMER_SCORES = (
    TautomerScore('benzoquinone', '[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]', 25),
    TautomerScore('oxim', '[#6]=[N][OH]', 4),
    TautomerScore('C=O', '[#6]=,:[#8]', 2),
    TautomerScore('N=O', '[#7]=,:[#8]', 2),
    TautomerScore('P=O', '[#15]=,:[#8]', 2),
    TautomerScore('C=hetero', '[C]=[!#1;!#6]', 1), #Modified according to RDKit
    TautomerScore("C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2), #Modified according to RDKit
    TautomerScore("aromatic C = exocyclic N", "[c]=!@[N]", -1), #Modified according to RDKit
    TautomerScore('methyl', '[CX4H3]', 1),
    TautomerScore('guanidine terminal=N', '[#7]C(=[NR0])[#7H0]', 1),
    TautomerScore('guanidine endocyclic=N', '[#7;R][#6;R]([N])=[#7;R]', 2),
    TautomerScore('aci-nitro', '[#6]=[N+]([O-])[OH]', -4),
    TautomerScore('amide', '[NH1]-[#6]=[#8]', 1), # I added here
)

#: The default value for the maximum number of tautomers to enumerate, a limit to prevent combinatorial explosion.
# --- Constants ---
DEFAULT_MAX_TAUTOMERS = 1000
DEFAULT_MAX_TRANSFORMS = 1000
_CIPCode = "_CIPCode"
_isotopicHs = "_isotopicHs"

   
# --- RDKit C++ style Sanitization Flags ---
SANITIZE_OPTIONS = (
    SanitizeFlags.SANITIZE_KEKULIZE |
    SanitizeFlags.SANITIZE_SETAROMATICITY |
    SanitizeFlags.SANITIZE_SETCONJUGATION |
    SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
    SanitizeFlags.SANITIZE_ADJUSTHS
)

class TautomerEnumerator(object):
    """
    Handles tautomer enumeration, scoring, and canonicalization.

    Provides methods:
    - enumerate: Generates tautomers and tracks modifications.
    - ScoreTautomer: Scores a single molecule based on defined rules.
    - Canonicalize: Enumerates, scores, selects the best tautomer.
    """

    def __init__(self,
                 transforms=TAUTOMER_TRANSFORMS,
                 scores=TAUTOMER_SCORES,
                 score_func=None,          # Optional custom scoring function for canonicalize
                 max_tautomers=DEFAULT_MAX_TAUTOMERS,
                 max_transforms=DEFAULT_MAX_TRANSFORMS,
                 remove_sp3_stereo=True,   # Corresponds to C++ tautomerRemoveSp3Stereo
                 remove_bond_stereo=True,  # Corresponds to C++ tautomerRemoveBondStereo
                 remove_isotopic_hs=True, # Corresponds to C++ tautomerRemoveIsotopicHs
                 reassign_stereo=True,     # Controls stereo reassignment *during* enumeration
                 callback=None,            # C++ callback equivalent for enumerate
                 debug=False):
        """
        Initialize the TautomerEnumerator.

        :param transforms: A list of TautomerTransforms to use.
        :param scores: A list of TautomerScores for default scoring.
        :param score_func: A custom scoring function (Mol -> int) for canonicalize.
                           Overrides default scoring logic if provided.
        :param max_tautomers: Max number of unique tautomers to generate.
        :param max_transforms: Max number of transform applications allowed.
        :param remove_sp3_stereo: Clear chiral tags on sp3 atoms involved in tautomerism.
        :param remove_bond_stereo: Clear bond stereo on bonds involved in tautomerism.
        :param remove_isotopic_hs: Remove isotopic H information from atoms involved.
        :param reassign_stereo: If True, run AssignStereochemistry on generated tautomers
                                 after specific adjustments during enumeration.
        :param callback: A function called during enumeration. If it returns False,
                         enumeration stops. Signature: callback(original_mol,
                         current_tautomers_dict, modified_atoms_set, modified_bonds_set) -> bool
        :param debug: Enable debug logging.
        """
        self.transforms = [t for t in transforms if t.tautomer is not None]
        self.scores = [s for s in scores if s.smarts is not None] # Also filter invalid scores
        self.score_func = score_func
        self.max_tautomers = max_tautomers
        self.max_transforms = max_transforms
        self.remove_sp3_stereo = remove_sp3_stereo
        self.remove_bond_stereo = remove_bond_stereo
        self.remove_isotopic_hs = remove_isotopic_hs
        self.reassign_stereo = reassign_stereo # Used by enumerate & _adjust_stereo
        self.callback = callback
        self.debug = debug

        self._num_atoms = 0 # Internal state, set during enumerate

    # --- Private Helper for Stereo/Iso Adjustments ---
    def _adjust_stereo_and_iso_hs(self, original_mol, taut_mol, modified_atoms, modified_bonds):
        """Internal method to adjust stereochemistry and isotopic Hs based on flags."""
        # --- (Code for _adjust_stereo_and_iso_hs is identical to the previous version) ---
        # --- It uses self.remove_sp3_stereo, self.remove_bond_stereo, etc. ---
        if not modified_atoms and not modified_bonds:
            return # Nothing to do if no modifications yet

        # --- Adjust Atom Stereo and Isotopes ---
        for atom_idx in range(self._num_atoms):
            if atom_idx not in modified_atoms:
                continue

            orig_atom = original_mol.GetAtomWithIdx(atom_idx)
            taut_atom = taut_mol.GetAtomWithIdx(atom_idx)

            # Clear chiral tag on sp2 atoms OR if remove_sp3_stereo is True
            if taut_atom.GetHybridization() == HybridizationType.SP2 or self.remove_sp3_stereo:
                if taut_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    taut_atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
                    if taut_atom.HasProp(_CIPCode):
                        taut_atom.ClearProp(_CIPCode)
            else: # Preserve original stereo if sp3 and remove_sp3_stereo is False
                 if taut_atom.GetChiralTag() != orig_atom.GetChiralTag():
                    taut_atom.SetChiralTag(orig_atom.GetChiralTag())
                    if orig_atom.HasProp(_CIPCode):
                         taut_atom.SetProp(_CIPCode, orig_atom.GetProp(_CIPCode))
                    elif taut_atom.HasProp(_CIPCode): # Clear if original didn't have it
                        taut_atom.ClearProp(_CIPCode)


            # Remove isotopic Hs if requested OR if atom now has no hydrogens
            if taut_atom.HasProp(_isotopicHs):
                if self.remove_isotopic_hs or taut_atom.GetTotalNumHs() == 0:
                   taut_atom.ClearProp(_isotopicHs)

        # --- Adjust Bond Stereo ---
        bonds_to_clear_dirs = set()
        for bond_idx in range(original_mol.GetNumBonds()):
            if bond_idx not in modified_bonds:
                continue

            orig_bond = original_mol.GetBondWithIdx(bond_idx)
            taut_bond = taut_mol.GetBondWithIdx(bond_idx)
            clear_dirs_for_this_bond = False
            if orig_bond.GetBondType() == BondType.DOUBLE and orig_bond.GetStereo() > BondStereo.STEREOANY:
                 clear_dirs_for_this_bond = True

            # Clear bond stereo if bond type changes OR remove_bond_stereo is True
            if taut_bond.GetBondType() != BondType.DOUBLE or self.remove_bond_stereo:
                if taut_bond.GetStereo() != BondStereo.STEREONONE:
                    taut_bond.SetStereo(BondStereo.STEREONONE)
                if clear_dirs_for_this_bond:
                     for atom_idx in (orig_bond.GetBeginAtomIdx(), orig_bond.GetEndAtomIdx()):
                        atom = original_mol.GetAtomWithIdx(atom_idx)
                        for nbr_bond in atom.GetBonds():
                            if nbr_bond.GetBondDir() in {BondDir.ENDDOWNRIGHT, BondDir.ENDUPRIGHT}:
                                bonds_to_clear_dirs.add(nbr_bond.GetIdx())
            else: # Preserve original bond stereo if remains DOUBLE and remove_bond_stereo is False
                 if taut_bond.GetStereo() != orig_bond.GetStereo():
                    taut_bond.SetStereo(orig_bond.GetStereo())
                 orig_stereo_atoms = orig_bond.GetStereoAtoms()
                 if len(orig_stereo_atoms) == 2: pass # Rely on SetStereo for now

                 if clear_dirs_for_this_bond:
                     for atom_idx in (orig_bond.GetBeginAtomIdx(), orig_bond.GetEndAtomIdx()):
                        atom = original_mol.GetAtomWithIdx(atom_idx)
                        for nbr_bond in atom.GetBonds():
                             if nbr_bond.GetBondDir() in {BondDir.ENDDOWNRIGHT, BondDir.ENDUPRIGHT}:
                                 taut_nbr_bond = taut_mol.GetBondWithIdx(nbr_bond.GetIdx())
                                 if taut_nbr_bond.GetBondDir() != nbr_bond.GetBondDir():
                                     taut_nbr_bond.SetBondDir(nbr_bond.GetBondDir())

        # Clear directions on associated bonds if the double bond stereo was cleared
        for bond_idx in bonds_to_clear_dirs:
             taut_bond = taut_mol.GetBondWithIdx(bond_idx)
             if taut_bond.GetBondDir() != BondDir.NONE:
                 taut_bond.SetBondDir(BondDir.NONE)

        # Final Reassignment based on instance flag (self.reassign_stereo)
        if self.reassign_stereo:
            try:
                Chem.AssignStereochemistry(taut_mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
            except TypeError: # Fallback for older RDKit
                 Chem.AssignStereochemistry(taut_mol, cleanIt=True, force=True)
        else:
            taut_mol.SetProp("_StereochemDone", "1")

    # --- Public Methods ---
    def enumerate(self, mol):
        """
        Enumerate tautomers for the input molecule.

        Applies transformations and handles stereo/isotopes according to
        the settings provided during initialization (e.g., self.reassign_stereo).

        :param mol: The input molecule (rdkit.Chem.rdchem.Mol).
        :return: A tuple containing:
                 (list of tautomer molecules,
                  set of modified atom indices,
                  set of modified bond indices)
        :rtype: tuple(list[rdkit.Chem.rdchem.Mol], set[int], set[int])
        """
        # --- (Code for enumerate is identical to the previous version) ---
        # --- It uses self.reassign_stereo directly via _adjust_stereo_and_iso_hs ---
        if not mol:
            return ([], set(), set())

        self._num_atoms = mol.GetNumAtoms()
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        tautomers = {smiles: copy.deepcopy(mol)}
        modified_atoms = set()
        modified_bonds = set()
        kekulized_mols = {}
        try:
            kekulized_mol = copy.deepcopy(mol)
            Chem.Kekulize(kekulized_mol)
            kekulized_mols[smiles] = kekulized_mol
        except Exception as e:
            log.warning(f"Kekulization failed for input molecule {smiles}: {e}. Skipping.")
            # Ensure original molecule respects reassign_stereo flag if returned directly
            # Although _adjust_stereo... expects modifications, let's apply AssignStereo if needed
            mol_copy = copy.deepcopy(mol)
            if self.reassign_stereo:
                 try: Chem.AssignStereochemistry(mol_copy, cleanIt=True, force=True, flagPossibleStereoCenters=True)
                 except TypeError: Chem.AssignStereochemistry(mol_copy, cleanIt=True, force=True)
            else: mol_copy.SetProp("_StereochemDone", "1")
            return ([mol_copy], set(), set())

        processed_smiles = set()
        n_transforms = 0
        status = "completed"

        while len(processed_smiles) < len(tautomers):
            if len(tautomers) >= self.max_tautomers: status = "max_tautomers"; break
            if n_transforms >= self.max_transforms: status = "max_transforms"; break
            if self.callback and not self.callback(mol, tautomers, modified_atoms, modified_bonds): status = "canceled"; break

            current_smiles = sorted([s for s in tautomers if s not in processed_smiles])[0]
            current_kekulized = kekulized_mols[current_smiles]
            processed_smiles.add(current_smiles)
            if self.debug: log.debug(f"Processing tautomer: {current_smiles}")

            for transform in self.transforms:
                if not transform.tautomer: continue
                matches = current_kekulized.GetSubstructMatches(transform.tautomer, uniquify=True)
                if not matches: continue

                # Check limits/callback before processing matches
                if len(tautomers) >= self.max_tautomers: status = "max_tautomers"; break
                if n_transforms >= self.max_transforms: status = "max_transforms"; break
                if self.callback and not self.callback(mol, tautomers, modified_atoms, modified_bonds): status = "canceled"; break

                for match in matches:
                    # Check limits/callback again inside match loop
                    if len(tautomers) >= self.max_tautomers: status = "max_tautomers"; break
                    if n_transforms >= self.max_transforms: status = "max_transforms"; break
                    if self.callback and not self.callback(mol, tautomers, modified_atoms, modified_bonds): status = "canceled"; break

                    n_transforms += 1
                    if self.debug: log.debug(f"Applying rule '{transform.name}' to {current_smiles} (match: {match})")
                    product = copy.deepcopy(current_kekulized)

                    # --- (Apply H shifts, bond changes, charge changes - Code identical) ---
                    try:
                        first_idx, last_idx = match[0], match[-1]
                        first, last = product.GetAtomWithIdx(first_idx), product.GetAtomWithIdx(last_idx)
                        first.SetNumExplicitHs(max(0, first.GetTotalNumHs() - 1))
                        last.SetNumExplicitHs(last.GetTotalNumHs() + 1)
                        first.SetNoImplicit(True); last.SetNoImplicit(True)
                        modified_atoms.update([first_idx, last_idx])

                        for bi, pair in enumerate(pairwise(match)):
                            bond = product.GetBondBetweenAtoms(pair[0], pair[1])
                            if not bond: raise ValueError(f"Bond not found {pair}")
                            bond_idx = bond.GetIdx()
                            modified_bonds.add(bond_idx)
                            if transform.bonds: bond.SetBondType(transform.bonds[bi])
                            else: bond.SetBondType(BondType.DOUBLE if bond.GetBondType() == BondType.SINGLE else BondType.SINGLE)

                        if transform.charges:
                            if len(transform.charges) != len(match): raise ValueError("Charge length mismatch")
                            for ci, atom_idx in enumerate(match):
                                atom = product.GetAtomWithIdx(atom_idx)
                                atom.SetFormalCharge(atom.GetFormalCharge() + transform.charges[ci])
                                modified_atoms.add(atom_idx)

                    except (IndexError, ValueError, RuntimeError) as e:
                         log.warning(f"Error applying transform {transform.name}: {e}. Match: {match}. Skipping.")
                         continue # Skip this match

                    # --- Sanitize and Finalize Product ---
                    try:
                        Chem.SanitizeMol(product, sanitizeOps=SANITIZE_OPTIONS)
                        for bond_idx in range(product.GetNumBonds()):
                             if bond_idx not in modified_bonds:
                                 orig_bond = mol.GetBondWithIdx(bond_idx)
                                 prod_bond = product.GetBondWithIdx(bond_idx)
                                 if (orig_bond.GetBondType() != prod_bond.GetBondType() or
                                     orig_bond.GetIsAromatic() != prod_bond.GetIsAromatic()):
                                     modified_bonds.add(bond_idx)
                                     if self.debug: log.debug(f"Bond {bond_idx} modified by sanitization.")

                        # Adjust stereo using the *instance's* reassign_stereo flag
                        self._adjust_stereo_and_iso_hs(mol, product, modified_atoms, modified_bonds)

                        prod_smiles = Chem.MolToSmiles(product, isomericSmiles=True)
                        if prod_smiles not in tautomers:
                            if self.debug: log.debug(f"New tautomer: {prod_smiles} (from rule {transform.name})")
                            tautomers[prod_smiles] = product
                            try:
                                kekulized_prod = copy.deepcopy(product)
                                Chem.Kekulize(kekulized_prod)
                                kekulized_mols[prod_smiles] = kekulized_prod
                            except Exception as e:
                                log.warning(f"Kekulization failed for generated tautomer {prod_smiles}: {e}. Removing.")
                                if prod_smiles in tautomers: del tautomers[prod_smiles]

                    except Exception as e:
                        log.warning(f"Error sanitizing/finalizing product from rule {transform.name}: {e}. Skipping product.")

                if status != "completed": break # Break from transform loop
            if status != "completed": break # Break from while loop

        # --- Final Return ---
        if status != "completed" and self.debug:
            log.warning(f"Tautomer enumeration stopped: {status}")

        return list(tautomers.values()), modified_atoms, modified_bonds

    
    def Canonicalize(self, mol):
        """
        Return the canonical tautomer for the input molecule.

        Internally calls self.enumerate (temporarily overriding reassign_stereo=False),
        scores the results using self.ScoreTautomer or self.score_func, selects
        the best, and performs final stereo assignment.

        :param mol: The input molecule (rdkit.Chem.rdchem.Mol).
        :return: The canonical tautomer molecule (rdkit.Chem.rdchem.Mol) or None if input is invalid.
        """
        if not mol:
            return None

        # --- Enumerate Tautomers with reassign_stereo temporarily off ---
        original_reassign_setting = self.reassign_stereo
        self.reassign_stereo = False # Override for internal enumeration step
        canonical_mol = None
        try:
            results = self.enumerate(mol)
            tautomers_list = results[0]
        except Exception as e:
            log.error(f"Exception during tautomer enumeration for canonicalization: {e}")
            tautomers_list = [] # Ensure list exists even on error
        finally:
            self.reassign_stereo = original_reassign_setting # Restore original setting

        # --- Select or Score ---
        if not tautomers_list:
            log.warning("No tautomers found during canonicalization, returning original molecule.")
            canonical_mol = copy.deepcopy(mol) # Use copy
        elif len(tautomers_list) == 1:
            if self.debug: log.debug("Only one tautomer found, selecting it as canonical.")
            canonical_mol = copy.deepcopy(tautomers_list[0]) # Use copy
        else:
            # --- Score Tautomers ---
            best_score = -float('inf')
            best_smiles = ""
            best_mol_ref = None # Reference to the best mol in the list

            # Choose scoring function: custom provided or default instance method
            scoring_function = self.score_func if self.score_func else self.ScoreTautomer

            if self.debug: log.debug(f"Scoring {len(tautomers_list)} tautomers...")

            for t in tautomers_list:
                try:
                    smiles = Chem.MolToSmiles(t, isomericSmiles=True)
                    score = scoring_function(t) # Call the chosen function
                    if self.debug: log.debug(f"  Tautomer: {smiles}, Score: {score}")

                    if score > best_score:
                        best_score = score
                        best_smiles = smiles
                        best_mol_ref = t
                    elif score == best_score and smiles < best_smiles:
                         best_smiles = smiles
                         best_mol_ref = t
                except Exception as e:
                    try: err_smiles = Chem.MolToSmiles(t, isomericSmiles=True)
                    except: err_smiles = "[Failed to generate SMILES]"
                    log.warning(f"Error scoring tautomer {err_smiles}: {e}")

            if best_mol_ref is None:
                log.error("Failed to select a best tautomer after scoring. Returning original.")
                canonical_mol = copy.deepcopy(mol)
            else:
                if self.debug: log.debug(f"Selected canonical tautomer: {best_smiles} (Score: {best_score})")
                canonical_mol = copy.deepcopy(best_mol_ref) # Make copy before final processing


        # --- Final Stereo Assignment on Canonical Choice ---
        if canonical_mol:
            try:
                 # Clear _StereochemDone if it was set by the enumeration step
                 if canonical_mol.HasProp("_StereochemDone"):
                    canonical_mol.ClearProp("_StereochemDone")

                 # Perform final assignment
                 Chem.AssignStereochemistry(canonical_mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
            except TypeError: # Fallback for older RDKit
                 if canonical_mol.HasProp("_StereochemDone"): canonical_mol.ClearProp("_StereochemDone")
                 Chem.AssignStereochemistry(canonical_mol, cleanIt=True, force=True)
            except Exception as e:
                log.error(f"Error during final AssignStereochemistry on canonical tautomer: {e}")
        else: # Handle case where even original molecule failed processing earlier
            log.error("Canonicalization resulted in no valid molecule.")
            return None


        return canonical_mol
    
    def ScoreTautomer(self, mol):
        """
        Calculate the score for a single tautomer molecule.

        Uses the default scoring logic based on rings, substructures (self.scores),
        and heteroatom penalties.

        :param mol: The tautomer molecule to score (rdkit.Chem.rdchem.Mol).
        :return: The calculated score (int).
        """
        score = 0
        if not mol: return -9999 # Handle None input

        # Ring scoring
        try:
            ring_info = mol.GetRingInfo()
            if ring_info and ring_info.NumRings() > 0:
                bond_rings = ring_info.BondRings()
                for ring_bond_indices in bond_rings:
                    is_all_aromatic = True
                    is_all_carbon = True
                    for b_idx in ring_bond_indices:
                        bond = mol.GetBondWithIdx(b_idx)
                        if not bond.GetIsAromatic(): is_all_aromatic = False; break
                        if bond.GetBeginAtom().GetAtomicNum() != 6 or bond.GetEndAtom().GetAtomicNum() != 6: is_all_carbon = False
                    if is_all_aromatic:
                        score += 100
                        if self.debug: log.debug("Score +100 (aromatic ring)")
                        if is_all_carbon:
                            score += 150
                            if self.debug: log.debug("Score +150 (carbocyclic aromatic ring)")
        except Exception as e: log.warning(f"Error during ring scoring: {e}")

        # Substructure scoring
        for tscore in self.scores: # Use self.scores
            if not tscore.smarts: continue
            try:
                matches = mol.GetSubstructMatches(tscore.smarts, uniquify=True)
                if matches:
                    match_count = len(matches)
                    score_change = match_count * tscore.score
                    score += score_change
                    if self.debug: log.debug(f"Score {score_change:+d} ({match_count}x {tscore.name})")
            except Exception as e: log.warning(f"Error during substructure scoring for '{tscore.name}': {e}")

        # Heteroatom H scoring
        for atom in mol.GetAtoms():
            anum = atom.GetAtomicNum()
            if anum in {15, 16, 34, 52}: # P, S, Se, Te
                try:
                    hs = atom.GetTotalNumHs()
                    if hs > 0:
                        score -= hs
                        if self.debug: log.debug(f"Score {-hs:+d} ({atom.GetSymbol()}-H bonds)")
                except Exception as e: log.warning(f"Error during heteroatom H scoring for atom {atom.GetIdx()}: {e}")

        return score