import logging
import os
import multiprocessing
import subprocess
import shutil
import random
import itertools
import tarfile, io
import time
import argparse


from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolAlign, rdMolTransforms, PropertyPickleOptions

#from msani.io.parsers import CustomHelpFormatter
from msani.conformers import utils, torsions

# Import C++ implementation - using direct transfer module (no SMILES conversion)
try:
    import stochastic_sampling_combined as cpp_stochastic
    CPP_AVAILABLE = True
    print(dir(cpp_stochastic))
except ImportError as e:
    print(f"C++ module not available, falling back to Python implementation. Error: {e}")
    CPP_AVAILABLE = False

# Global flag to control which implementation to use
USE_CPP = True  # Set to False to force Python implementation


def stochastic_sampling_v2(
                               mol,
                               angle_map,
                               score_map,
                               numConfs,
                               possible_numConfs,
                               importance_order,
                               window = 25,
                               max_attempts=50_000,
                               hetero_H_bonds = [],
                               timeout_conf = 2,
                               rmsd = 0.5, 
                               verbose = False,
                               ):
    """
    Perform stochastic sampling of conformers based on a given angle map and score map.
    
    This version automatically selects between C++ and Python implementations.
    The C++ implementation includes energy filtering and conformer management.

    Args:
        mol (Chem.Mol): The molecule to sample conformers for.
        angle_map (dict): A dictionary mapping bond indices to tuples of (bond, dihedral atoms, possible angles).
        score_map (dict): A dictionary mapping bond indices to scores for each possible angle.
        numConfs (int): The number of conformers to generate.
        possible_numConfs (int): The total number of possible conformations.
        importance_order (list): A list of weights for the importance of each bond.
        window (float): The energy window for accepting conformers.
        max_attempts (int): The maximum number of attempts to generate conformers.
        hetero_H_bonds (list): A list of tuples representing heteroatom-hydrogen bonds to consider.
        timeout_conf (int): The timeout in minutes for generating conformers (converted to seconds internally).
        product (list): A list to store the generated conformers and their energies.

    Returns:
        product (list): A list of tuples containing the generated conformers and their energies.
    """
    
    # Try to use C++ implementation first
    if CPP_AVAILABLE and USE_CPP:
        try:
            # Call C++ direct transfer implementation
            result = cpp_stochastic.stochastic_sampling_discrete(
                mol=mol,
                angle_map=angle_map,
                score_map=score_map,
                possible_numConfs=possible_numConfs,
                importance_order=importance_order.tolist(),
                window=window,
                max_attempts=max_attempts,
                hetero_H_bonds=hetero_H_bonds,
                timeout_conf=int(timeout_conf * 60),  # Convert minutes to seconds (int)
                rmsd=rmsd,
                numConfs=numConfs,
                
                randomSeed=42,
                clash_threshold=1.6,  # Default clash threshold
                verbose=verbose,
                mmff_variant="MMFF94s",
                eps=1.0,

            )

            
            # Direct transfer C++ implementation returns molecule with conformers
            if result is not None and hasattr(result, 'GetNumConformers'):
                return(result)
            
        except Exception as e:
            print(f"C++ implementation failed with error: {e}")
            pass
    
    # Fallback to original Python implementation
   # print("🐍 Using Python implementation")
    return []

def main():
    parser = argparse.ArgumentParser(description="Generate conformers for a given SMILES string/file.\nTwo-column files are required.",
                                     #formatter_class=CustomHelpFormatter,
                                     add_help = False)  # Suppress default -h/--help)
    parser.add_argument('--smiles', '-s', type = str, default = None, help='Input SMILES string.')
    parser.add_argument('--numconfs', '-nconfs', type=int, default=600, help='Number of conformers to generate (default: 2000).')
    parser.add_argument('--debug', '-d', action='store_true', help='Enable verbose output for debugging.')
    parser.add_argument('--timeout', '-to',type=int, default=2, help='Timeout in minutes for RDKit-based conformation generation.')
    parser.add_argument('--timeout_conf', '-toc', type=float, default=2, help='Timeout in minutes for conformational sampling (default: 2 minutes).')
    parser.add_argument('--energywindow', '-w',type=float, default=25.0, help='Energy window for conformer generation.')
    parser.add_argument('--numcores', '-j', type=int, default=4, help='Number of CPU cores to use (default: 4).')
    parser.add_argument('--timing', action='store_true', help='Log timing information for each step.')
    parser.add_argument('--nocleanup', action='store_false', dest='cleanup', help='Do not remove intermediate files after processing.')
    parser.add_argument('--rigid', '-r', type=str, default=None, help='SMILES/SMARTS for conformers to be aligned to.')
    parser.add_argument("--help", "-h", action="help", help="Show this help message and exit")
    parser.add_argument('--randomSeed', '-rs',type=int, default=42, help=argparse.SUPPRESS)
    parser.add_argument('--test', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--synthon', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--torsion', '-t', type=str, default=None, help='File containing custom torsion rules.')
    parser.add_argument('--rmsd', '-rmsd', type=float, default=0.5, help=argparse.SUPPRESS)
    args = parser.parse_args()

    mol = Chem.MolFromSmiles(args.smiles)
    mol_H = Chem.AddHs(mol)
    rdDistGeom.EmbedMolecule(mol_H, randomSeed=args.randomSeed)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol_H, maxIters=2000, mmffVariant='MMFF94s')
    Torlib = torsions.TorsionLibrary()
    amide_linkages = utils.find_amide(mol_H)
    ignoreTorlib = False
    rot_bonds = utils.getDihedralMatches(mol_H)
    possible_numConfs, angle_map, score_map, rot_bonds, hetero_H_bonds = utils.count_confs_by_rotbonds_v2(mol = mol_H,
                                                                                rot_bonds = rot_bonds,
                                                                                amide_bonds = amide_linkages,
                                                                                ignoretorlib = ignoreTorlib,
                                                                                torlib = Torlib,
                                                                                VERBOSE= args.debug,
                                                                                )
    importance_order = utils.get_importance_order(mol_H, rot_bonds)
    requested_num_confs = args.numconfs
    processing_mol = Chem.Mol(mol_H)
    product = stochastic_sampling_v2(mol = processing_mol,
                                    angle_map = angle_map,
                                    score_map = score_map,
                                    numConfs = requested_num_confs,
                                    possible_numConfs = possible_numConfs,
                                    importance_order =  importance_order, 
                                    window = args.energywindow, 
                                    max_attempts = 50_000,
                                    hetero_H_bonds = hetero_H_bonds, 
                                    timeout_conf = args.timeout_conf,  # Convert minutes to seconds (int)
                                    rmsd = args.rmsd,
                                    verbose = args.debug,
                                    )
    print(f"Generated {product.GetNumConformers()} conformers for {args.smiles}")
    with Chem.SDWriter("output.sdf") as writer:
        for conf in product.GetConformers():
            idx = conf.GetId()
            rdMolAlign.AlignMol(product, mol_H, prbCid=idx, refCid=0, atomMap=[(i,i) for i in [0, 1, 2]])
            writer.write(product, confId=idx)

        
if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(f"Total time: {end - start:.2f} seconds")