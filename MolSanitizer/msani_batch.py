"""
MolSanitizer.
"""

__author__ = "Israel Cabeza de Vaca Lopez, Thua-Phong Lam, Szymon Pach"
__place__ = "Jens Carlsson lab, Uppsala University, Sweden"
__license__ = "MIT"

import logging
logger = logging.getLogger('molsani')

import pandas as pd

import pathlib
import os
import time
import sys

from . import parsers
import subprocess

from rdkit import Chem
from rdkit import rdBase

slurm_header = '''#!/bin/bash
#SBATCH -A PROJECT_NAME
#SBATCH -n 1
#SBATCH -J msani_db2
#SBATCH -t 23:59:59
#SBATCH --mail-type=FAIL
#SBATCH --mem=64G
# '''

slurm_script='''
dirs=( $(cat dirlista) )
TASK_ID=${SLURM_ARRAY_TASK_ID}
smiles_file=${dirs[$TASK_ID]}

~/.conda/envs/msani -i $smiles_file'''

def parse_flags_single_job(args: dict):
    """Parse the flags for a single job

    Args:
        args (dict): Arguments from the command line

    Returns:
        str: The flags for a single job
    """
    flags = ''
    if args.enamine: flags += ' --enamine'
    if args.lazy: flags += ' --lazy'
    if args.removesalts: flags += ' --removesalts'
    if args.tautomers: flags += ' --tautomers'
    if args.pains: flags += ' --pains'
    if args.unwanted is not None: flags += f' --unwanted {" ".join(args.unwanted)}'
    if args.stereoisomers: flags += ' --stereoisomers'
    if args.protonation: flags += ' --protonation'
    if args.db2: flags += ' --db2'
    if not(args.cleanup): flags += ' --nocleanup'
    if args.custom is not None: flags += f' --custom {args.custom}'

    if args.max_isomers != 0: flags += f' --max_isomers {args.max_isomers}'
    if args.numconfs != 2000: flags += f' --numconfs {args.numconfs}'
    if args.randomSeed != 42: flags += f' --randomSeed {args.randomSeed}'
    return flags
def write_single_job_script(slurm_header: str, slurm_script: str):
    """Write the script for a single job

    Args:
        slurm_header (str): The header for the SLURM script
        slurm_script (str): The script for the SLURM job

    Returns:
        None
    """
    with open('submit_msani.sh', 'w') as f:
        f.write(slurm_header)
        f.write(slurm_script)

def Split_Submit_jobs(args: dict):
    """Split the input files into chunks and submit jobs to the cluster

    Args:
        args (dict): Arguments from the command line

    Returns:
        None
    """
    # Replace the PROJECT_NAME with the project name for SLURM
    global slurm_header
    slurm_header = slurm_header.replace('PROJECT_NAME', args.proj_name)
    flags = parse_flags_single_job(args)
    global slurm_script
    slurm_script = slurm_script + flags
    print(slurm_script)
    for file in args.input_files:
        prefix = file.split('.')[0]
        subprocess.run(f"mkdir -p {prefix}", shell=True)
        subprocess.run(f"split -l {args.lines} -d {file} -a 3 {prefix}/in", shell=True)
        os.chdir(prefix)
        n_jobs = len(os.listdir())
        subprocess.run(f"ls in* > dirlista", shell=True)
        write_single_job_script(slurm_header, slurm_script)
        subprocess.run(f"sbatch --array=0-{n_jobs-1}%100 submit_msani.sh", shell=True)
        os.chdir('..')
        
def main():

    args = parsers.parseArguments_batch(sys.argv[1:])
    Split_Submit_jobs(args)




if __name__=="__main__":
    main()