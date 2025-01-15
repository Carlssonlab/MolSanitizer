"""
MolSanitizer in the batch mode.
"""

__author__ = "Thua-Phong Lam, Israel Cabeza de Vaca Lopez, Szymon Pach"
__place__ = "Jens Carlsson lab, Uppsala University, Sweden"
__license__ = "GPLv2"
__version__ = "0.2.2"



import os
import sys
import time
import math
import yaml
from . import parsers
import subprocess
from rdkit import rdBase

slurm_header = '''#!/bin/bash
#SBATCH -A PROJECT_NAME
#SBATCH -n 1
#SBATCH -J msani_3d
#SBATCH -t TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mem=64G
'''

slurm_script='''
dirs=( $(cat dirlista) )
TASK_ID=${SLURM_ARRAY_TASK_ID}
smiles_file=${dirs[$TASK_ID]}
ARRAY_ID=${SLURM_ARRAY_JOB_ID}

log_prefix=$(basename "$smiles_file")
log_prefix="${log_prefix%.*}"  # Remove the extension
log_file="${log_prefix}.log"
MSANI_PATH -i $smiles_file -j 2'''

cleanup_script ="""
# Get the number of tasks with the name msani_3d from the user's squeue
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
        echo "sbatch --array=${failed_tasks} submit_msani.sh" > RESUBMIT_FAILED_JOBS.txt
        echo "Instructions for resubmitting failed jobs written to RESUBMIT_FAILED_JOBS.txt"

        # Create a pattern to exclude failed task files with zero-padded IDs
        exclude_pattern=$(echo $failed_tasks | tr ',' '\n' | awk '{printf "in%03d.smi ", $1}' | tr '\n' ' ')
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
with open(os.path.join(os.path.dirname(__file__), 'msani_configurations.yaml')) as confFile:
    configurations = yaml.full_load(confFile)
    slurm_account = configurations['SLURM_ACCOUNT']
    time_limit = configurations['TIME_LIMIT']
    lines_per_job = configurations['LINES_PER_JOB']
    max_jobs = configurations['MAX_JOBS']
    max_array_size = configurations['MAX_ARRAY_SIZE']
    max_limit_project = configurations['MAX_LIMIT_PROJECT']
    timeout = configurations['TIMEOUT']
    use_corina = configurations['USE_CORINA']
    corina_exe = configurations['CORINA']
    energy_window = configurations['ENERGY_WINDOW']
    numconfs = configurations['NUMCONFS']
    max_stereoisomers = configurations['MAX_STEREOISOMERS']
    pH = configurations['PH']
    pH_range = configurations['PH_RANGE']

def count_lines_bash(file_path):
    result = subprocess.run(['wc', '-l', file_path], stdout=subprocess.PIPE)
    return int(result.stdout.split()[0])

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
    if not(args.taurdkit): flags += ' --notaurdkit'
    if args.pains: flags += ' --pains'
    if args.unwanted is not None: flags += f' --unwanted {" ".join(args.unwanted)}'
    if args.ha is not None: flags += f' --ha {args.ha}'
    if args.logp is not None: flags += f' --logp {args.logp}'
    if args.stereoisomers: 
        flags += ' --stereoisomers'
        if args.max_stereoisomers != max_stereoisomers: flags += f' --max_stereoisomers {args.max_stereoisomers}'

    if args.protonation: 
        flags += ' --protonation'
        if args.pH != pH: flags += f' --pH {args.pH}'
        if args.pH_range != pH_range: flags += f' --pH_range {args.pH_range}'

    if args.db2: flags += ' --db2'
    if not(args.cleanup): flags += ' --nocleanup'
    if args.debug: flags += ' --debug'
    if args.timing: flags += ' --timing'
    if args.custom is not None: flags += f' --custom ../{args.custom}'
    if args.energywindow != energy_window: flags += f' --energywindow {args.energywindow}'
    if args.numconfs != numconfs: flags += f' --numconfs {args.numconfs}'
    if args.randomSeed != 42: flags += f' --randomSeed {args.randomSeed}'
    if args.timeout != timeout: flags += f' --timeout {args.timeout}'
    if args.enrichment: flags += ' --enrichment'
    if args.corina: flags += f' --corina'
    if args.ignoretorlib: flags += ' -igtor'
    if not(args.neutralize): flags += ' --noneutralize'
    if args.synthon: flags += ' --synthon'
    return flags

def write_single_job_script(slurm_header: str, slurm_script: str):
    """Write the script for a single job

    Args:
        slurm_header (str): The header for the SLURM script
        slurm_script (str): The script for the SLURM job

    Returns:
        None
    """
    result = subprocess.run(['which', 'msani'], stdout=subprocess.PIPE, text=True)

    # Get the stdout from the result and strip any extra whitespace
    msani_path = result.stdout.strip()

    #print(f"msani_path: {msani_path}")
    slurm_script = slurm_script.replace('MSANI_PATH', msani_path)

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
    # Replace the PROJECT_NAME with the project name and time limit for SLURM
    global slurm_header
    slurm_header = slurm_header.replace('PROJECT_NAME', args.proj_name)
    slurm_header = slurm_header.replace('TIME_LIMIT', f'{args.timelimit}:00:00')
    flags = parse_flags_single_job(args)
    global slurm_script
    slurm_script = slurm_script + flags + remove_lock_files

    if args.cleanup: slurm_script += cleanup_script
    
    print(f"\nStarting MolSanitizer in batch mode\n")
    print(f"Using project name (-p): {args.proj_name}")
    print(f"Time limit for each job (-t): {args.timelimit} hours")
    print(f"Maximum number of jobs running parallelly (--max_jobs): {args.max_jobs} jobs")
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
            print(f"Exitting MolSanitizer...")
            return
        line_count = count_lines_bash(file)
        n_jobs += math.ceil(line_count/args.lines)
    print(f"Total number of jobs to submit: {n_jobs}\n")
    if n_jobs > max_array_size:
        print(f"Too many jobs to submit ({n_jobs}). Please increase the number of lines per job or decrease the number of input files")
        print(f"Exitting MolSanitizer...")
        return
    
    if current_running_jobs + n_jobs > max_limit_project:
        print(f"Current number of jobs running in the project {args.proj_name}: {current_running_jobs}")
        print(f"Total number of jobs to submit: {n_jobs}")
        print(f"Total number of jobs will exceed the limit of {max_limit_project} jobs in the project {args.proj_name}")
        print(f"Please wait for the current jobs to finish before submitting new jobs.")
        print(f"Exitting MolSanitizer...")
        return

    # Wait for 5 seconds before proceeding
    print("Waiting 5 seconds to review the configurations...")
    time.sleep(5)
    
    for file in args.input_files:
        prefix = file.split('.')[0]
        if os.path.exists(prefix):
            remove_folder = input(f"Folder {prefix} already exists. Do you want to remove it? (y/n): ")
            if remove_folder.lower() == 'y' or remove_folder.lower() == 'yes':
                print(f"Removing folder {prefix}...\n")
                subprocess.run(f"rm -rf {prefix}", shell=True)
            else:
                print(f"Exitting MolSanitizer...\n")
                return
        subprocess.run(f"mkdir -p {prefix}", shell=True)
        subprocess.run(f"split -l {args.lines} -d -a 3 --additional-suffix=.smi {file} {prefix}/in", shell=True)
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
        subprocess.run(f"sbatch --array=0-{n_jobs-1}%{args.max_jobs} submit_msani.sh", shell=True)
        os.chdir('..')
        
def main():

    args = parsers.parseArguments_batch(sys.argv[1:])
    rdkit_version = rdBase.rdkitVersion
    if rdkit_version != '2024.09.1':
        print('\n###########################################################')
        print('RDKit version 2024.09.1 is recommended for MolSanitizer.')
        print("Use 'conda install rdkit==2024.9.1' to avoid potential issues.")
        print('##############################################################\n')
        time.sleep(2)
    Split_Submit_jobs(args)


if __name__=="__main__":
    main()
