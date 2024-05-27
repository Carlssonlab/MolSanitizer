#!/usr/bin/env python3

"""
This is a script to standardize a SMILES file without considering organometallics.

[Flavio Ballante, CBCS-Karolinska Institutet - 2023]
"""

import argparse
import logging
import sys
import os
import pprint
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolfiles import SmilesWriter

def ParseArgs():
    """
    Parse the arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='source .smi file')
    parser.add_argument('-o', '--output', required=True, help='output .smi file')
    parser.add_argument('-f', '--field', required=True, help='use: _Name')
    parser.parse_args()
    parser.set_defaults(verbose=False)
    args = parser.parse_args()
    return args

def standardize(smiles):
    # follows the steps in
    # https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb
    # https://www.youtube.com/watch?v=eWTApNX8dJQ
    mol = Chem.MolFromSmiles(smiles)
     
    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.tautomerRemoveIsotopicHs = False
    clean_mol = rdMolStandardize.Cleanup(mol, params) 

     
    # if many fragments, get the "parent"
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol, params)
         
    # try to neutralize molecule
    uncharger = rdMolStandardize.Uncharger()
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
    
    # tautomer enumerator
    te = rdMolStandardize.TautomerEnumerator(params) 
    taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)
     
    return taut_uncharged_parent_clean_mol

def main():
    """
    MAIN
    """
    args = ParseArgs()
    path = os.getcwd()

    #Sc, Y, In, Sn, W, Ac are not included for now (need to check how they can hit other structures).
    organometallics=['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Ga', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',\
'Cd', 'La', 'Hf ', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Ce', 'Pr', 'Nd',\
'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf ', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Ge', 'Sb']
    print(rdBase.rdkitVersion)

    mols = [mol for mol in Chem.SmilesMolSupplier(args.input) if mol != None]
    with open(args.output, 'w') as outfile, open(args.input,'r'): as infile

        outfile.write('SMILES'+ ' ' + 'ID\n')
        
        for line in infile:

                mol = Chem.SmilesToMol(line.split()[0])

                if mol == None: 

			print(f'Error with SMILES:{line}')
                        continue

                #for mol in mols:
                match = True in [item in Chem.MolToSmiles(mol) for item in organometallics]

                if match == False:
                    outfile.write("{} {}\n".format(Chem.MolToSmiles(standardize(Chem.MolToSmiles(mol))), mol.GetProp(args.field)))
                elif match == True:
                    outfile.write("{} {}\n".format(Chem.MolToSmiles(mol), mol.GetProp(args.field)))

if __name__ == "__main__":
    main()
