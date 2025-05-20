"""
EirVS in the batch mode.
"""

__author__ = "Thua-Phong Lam, Israel Cabeza de Vaca Lopez, Szymon Pach"
__place__ = "Jens Carlsson lab, Uppsala University, Sweden"
__license__ = "GPLv2"
__version__ = "0.2.3"



import os
import sys
import time
import math
import yaml
import subprocess

from pathlib import Path
from .io import parsers
from rdkit import rdBase

slurm_header = '''#!/bin/bash
#SBATCH -A PROJECT_NAME
#SBATCH -n 1
#SBATCH -J eirvs_3d
#SBATCH -t TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mem=MEMORY
'''

slurm_script='''
dirs=( $(cat dirlista) )
TASK_ID=${SLURM_ARRAY_TASK_ID}
smiles_file=${dirs[$TASK_ID]}
ARRAY_ID=${SLURM_ARRAY_JOB_ID}

log_prefix=$(basename "$smiles_file")
log_prefix="${log_prefix%.*}"  # Remove the extension
log_file="${log_prefix}.log"
EIRVS_PATH -i $smiles_file -j 2'''

cleanup_script ="""
task_count=$(ls *.lock 2>/dev/null | wc -l)
echo "$task_count remaining jobs in the queue."

# Remove the SLURM output file for this array task
slurm_out_file="slurm-${ARRAY_ID}_${TASK_ID}.out"
cat $slurm_out_file >> "$log_file"
rm -f "$slurm_out_file"

# If the count is equal to 0 (last job in the array), perform cleanup and check for failed tasks
if [ "$task_count" -eq 0 ]; then
    echo "Proceeding with cleanup..."

    # Find and delete empty directories in the 3d directory
    find 3d -type d -empty -delete

    # Find and delete empty directories in the solv directory
    find solv -type d -empty -delete
    mkdir -p in/processed in/removed log
    mv *.log log
    mv in*_rejected* in/removed -f 2>/dev/null
    mv in*_clean* in/processed
    find in/removed -type d -empty -delete

    # Check for failed tasks using sacct
    failed_tasks=$(sacct -j "${ARRAY_ID}" --format='JobID%30,State' --noheader | grep 'NODE_FAIL' | awk -F_ '{print $2}' | awk '{print $1}' | tr '\n' ',' | sed 's/,$//')

    if [ -n "$failed_tasks" ]; then
        echo "Failed tasks detected: $failed_tasks"
        echo "sbatch --array=${failed_tasks} submit_eirvs.sh" > RESUBMIT_FAILED_JOBS.txt
        echo "Instructions for resubmitting failed jobs written to RESUBMIT_FAILED_JOBS.txt"

        # Create a pattern to exclude failed task files with zero-padded IDs
        exclude_pattern=$(echo $failed_tasks | tr ',' '\n' | awk '{printf "in%04d.smi ", $1}' | tr '\n' ' ')
        echo "Excluding files: $exclude_pattern"

        # Move all `in*.smi` files except those corresponding to failed tasks
        for file in in*.smi; do
            if [[ ! $exclude_pattern =~ $(basename "$file") ]]; then
                mv "$file" in
            fi
        done
    else
        echo "No failed tasks detected."
        mv in*.smi in
    fi

    echo "Cleanup complete."
fi
"""

remove_lock_files = '''
rm -f "${log_prefix}.lock"
'''
with open(os.path.join(os.path.dirname(__file__), 'eirvs_configurations.yaml')) as confFile:
    configurations = yaml.full_load(confFile)
    slurm_account = configurations['SLURM_ACCOUNT']
    time_limit = configurations['TIME_LIMIT']
    lines_per_job = configurations['LINES_PER_JOB']
    max_jobs = configurations['MAX_JOBS']
    max_array_size = configurations['MAX_ARRAY_SIZE']
    max_limit_project = configurations['MAX_LIMIT_PROJECT']
    timeout = configurations['TIMEOUT']
    corina_exe = configurations['CORINA']
    energy_window = configurations['ENERGY_WINDOW']
    numconfs = configurations['NUMCONFS']
    max_stereoisomers = configurations['MAX_STEREOISOMERS']
    pH = configurations['PH']
    pH_range = configurations['PH_RANGE']

def count_lines_bash(file_path):
    result = subprocess.run(['wc', '-l', file_path], stdout=subprocess.PIPE)
    return int(result.stdout.split()[0])

def parse_flags_single_job(args: dict, parser):
    """Parse the flags for a single job

    Args:
        args (dict): Arguments from the command line

    Returns:
        str: The flags for a single job
    """
    flags = []
    omitted_args = ["input_files", "smiles", "proj_name", "timelimit", "lines", "max_jobs", "help"]

    for action in parser._actions:
        arg = action.dest  # Argument name
        if not hasattr(args, arg) or arg in omitted_args:  # Skip help and internal arguments
            continue
        current_value = getattr(args, arg)  # Current value in the Namespace
        default_value = action.default      # Default value from the parser
        
        # Skip if the value is the same as the default or if the argument isn't specified
        if current_value == default_value:
            continue

        if isinstance(current_value, bool):
            # Boolean flags
            # Handle the counterintuitive flags
            if arg in ["taurdkit", "cleanup", "neutralize"] and not current_value:
                flags.append(f"--no{arg}")
            # Handle the rest of the flags
            elif current_value:
                flags.append(f"--{arg}")
        elif isinstance(current_value, list):
            # List arguments
            value_str = " ".join(map(str, current_value))
            flags.append(f"--{arg} {value_str}")
        else:
            # Other arguments
            flags.append(f"--{arg} '{current_value}'")

    
    flags = " " + " ".join(flags)
    return flags

def write_single_job_script(slurm_header: str, slurm_script: str):
    """Write the script for a single job

    Args:
        slurm_header (str): The header for the SLURM script
        slurm_script (str): The script for the SLURM job

    Returns:
        None
    """
    result = subprocess.run(['which', 'eirvs'], stdout=subprocess.PIPE, text=True)

    # Get the stdout from the result and strip any extra whitespace
    eirvs_path = result.stdout.strip()

    slurm_script = slurm_script.replace('EIRVS_PATH', eirvs_path)

    with open('submit_eirvs.sh', 'w') as f:
        f.write(slurm_header)
        f.write(slurm_script)

def test_batch_mode(args: dict):
    global slurm_header
    global slurm_script
    file = Path(args.prefix) / args.input_files[0]
    prefix = Path(args.prefix) / file.stem  # Ensure prefix is within temp_dir

    subprocess.run(f"mkdir -p {prefix}", shell=True)
    subprocess.run(f"split -l {args.lines} -d -a 4 --additional-suffix=.smi {file} {prefix}/in", shell=True)
    os.chdir(prefix)
    subprocess.run(f"ls in* > dirlista", shell=True)
    write_single_job_script(slurm_header, slurm_script)

def Split_Submit_jobs(args: dict, parser):
    """Split the input files into chunks and submit jobs to the cluster

    Args:
        args (dict): Arguments from the command line

    Returns:
        None
    """
    # Replace the PROJECT_NAME with the project name and time limit for SLURM
    
    global slurm_header
    slurm_header = slurm_header.replace('PROJECT_NAME', args.proj_name)
    slurm_header = slurm_header.replace('TIME_LIMIT', f'{args.timelimit}:00:00')
    if not(args.gen3d): slurm_header = slurm_header.replace('eirvs_3d', 'eirvs_2d')
    if args.lines >= 250_000: slurm_header = slurm_header.replace('MEMORY', '12G')
    elif args.lines >= 100_000: slurm_header = slurm_header.replace('MEMORY', '6G')
    else: slurm_header = slurm_header.replace('MEMORY', '4G')
    flags = parse_flags_single_job(args, parser)
    global slurm_script
    slurm_script = slurm_script + flags + remove_lock_files

    if args.cleanup: slurm_script += cleanup_script
    if args.test: 
        test_batch_mode(args)
    else:
        print(f"\nStarting EirVS in batch mode\n")
        print(f"Using project name (-A): {args.proj_name}")
        print(f"Time limit for each job (-tl): {args.timelimit} hours")
        print(f"Maximum number of jobs in an array: {max_array_size} jobs")
        print(f"Maximum number of jobs running parallelly (-mj): {args.max_jobs} jobs")
        print(f"Number of compounds per job (-l): {args.lines} lines\n")

        result = subprocess.run(f'squeue -A {args.proj_name} -r | wc -l', shell=True, stdout=subprocess.PIPE, text=True)
        try: 
            current_running_jobs = int(result.stdout.strip())
        except:
            current_running_jobs = 0
            pass
        
        n_jobs = 0
        for file in args.input_files:
            if not os.path.exists(file):
                print(f"File {file} does not exist. Please check the path and try again.")
                print(f"Exitting EirVS...")
                return
            line_count = count_lines_bash(file)
            n_jobs += math.ceil(line_count/args.lines)
        print(f"Total number of jobs to submit: {n_jobs}\n")
        if n_jobs > max_array_size:
            print(f"Too many jobs to submit ({n_jobs}). Please increase the number of lines per job or decrease the number of input files")
            print(f"Exitting EirVS...")
            return
        
        if current_running_jobs + n_jobs > max_limit_project:
            print(f"Current number of jobs running in the project {args.proj_name}: {current_running_jobs}")
            print(f"Total number of jobs to submit: {n_jobs}")
            print(f"Total number of jobs will exceed the limit of {max_limit_project} jobs in the project {args.proj_name}")
            print(f"Please wait for the current jobs to finish before submitting new jobs.")
            print(f"Exitting EirVS...")
            return

        # Wait for 5 seconds before proceeding
        print("Waiting 5 seconds to review the configurations...")
        time.sleep(5)
        print('Submitting jobs...\n')
        for file in args.input_files:
            prefix = file.stem
            if os.path.exists(prefix):
                remove_folder = input(f"Folder {prefix} already exists. Do you want to remove it? (y/n): ")
                if remove_folder.lower() in ['y','yes']:
                    print(f"Removing folder {prefix}...\n")
                    subprocess.run(f"rm -rf {prefix}", shell=True)
                else:
                    print(f"Exitting EirVS...\n")
                    return
            subprocess.run(f"mkdir -p {prefix}", shell=True)
            subprocess.run(f"split -l {args.lines} -d -a 4 --additional-suffix=.smi {file} {prefix}/in", shell=True)
            os.chdir(prefix)
            subprocess.run(f"ls in* > dirlista", shell=True)
            n_jobs = sum(1 for line in open('dirlista'))
            with open('dirlista') as f:
                for line in f:
                    jobname = line.strip().split('.')[0]
                    with open(f'{jobname}.lock', 'w') as lock:
                        lock.write('')
            print(f"Submitting {n_jobs} jobs\n")
            write_single_job_script(slurm_header, slurm_script)
            subprocess.run(f"sbatch --array=0-{n_jobs-1}%{args.max_jobs} submit_eirvs.sh", shell=True)
            os.chdir('..')
        
def main():

    args, parser = parsers.parseArguments(sys.argv[1:], batch_mode=True)
    rdkit_version = rdBase.rdkitVersion
    if args.version:
        print(f"EirVS version: {__version__}")
        print(f"RDKit version: {rdkit_version}")
        return
    Split_Submit_jobs(args, parser)


if __name__=="__main__":
    main()
