"""
MolSanitizer in the batch mode.
"""

from msani import __version__


import os
import sys
import time
import math
import shutil
import subprocess
import tempfile
import gzip
import bz2
import lzma

from yaml import safe_load
from pathlib import Path
from rdkit import rdBase

from msani.io import parsers
from msani.io.readers import detect_input_compression


_DECOMPRESSOR_COMMANDS = {
    'gzip': 'gzip',
    'bz2': 'bzip2',
    'xz': 'xz',
}

_DECOMPRESSOR_OPEN = {'gzip': gzip.open, 'bz2': bz2.open, 'xz': lzma.open}

logo=r""" __  __         _  _____                _  _    _                 
|  \/  |       | |/  ___|              (_)| |  (_)                
| .  . |  ___  | |\ `--.   __ _  _ __   _ | |_  _  ____ ___  _ __ 
| |\/| | / _ \ | | `--. \ / _` || '_ \ | || __|| ||_  // _ \| '__|
| |  | || (_) || |/\__/ /| (_| || | | || || |_ | | / /|  __/| |   
\_|  |_/ \___/ |_|\____/  \__,_||_| |_||_| \__||_|/___|\___||_|   

                                                In the batch mode
"""

slurm_header = '''#!/bin/bash
PROJECT_NAME
PARTITION_NAME
#SBATCH -n 1
#SBATCH -J msani_3d
#SBATCH -t TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mem=MEMORY
'''

slurm_header_whole_node = '''#!/bin/bash
PROJECT_NAME
PARTITION_NAME
#SBATCH --nodes=1
#SBATCH --ntasks=NODE_CORES
#SBATCH -J msani_3d
#SBATCH -t TIME_LIMIT
#SBATCH --mail-type=FAIL
'''

slurm_script='''
TASK_ID=${SLURM_ARRAY_TASK_ID}
ARRAY_ID=${SLURM_ARRAY_JOB_ID}

offsets=( $(cat input.offsets) )
start_offset=${offsets[$TASK_ID]}
end_offset=${offsets[$(( TASK_ID + 1 ))]}
byte_count=$(( end_offset - start_offset ))

printf -v log_prefix "in%04d" "$TASK_ID"
log_file="${log_prefix}.log"
scratch_parent=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
task_dir=$(mktemp -d "${scratch_parent}/msani-${ARRAY_ID}-${TASK_ID}-XXXXXX")
smiles_file="${task_dir}/${log_prefix}.smi"
trap 'rm -rf "$task_dir"' EXIT

dd if=input.data of="$smiles_file" iflag=skip_bytes,count_bytes \\
   skip="$start_offset" count="$byte_count" status=none
MSANI_PATH -i "$smiles_file" -j 1'''

slurm_script_whole_node='''
offsets=( $(cat input.offsets) )
n_total=$(( ${#offsets[@]} - 1 ))
ARRAY_ID=${SLURM_ARRAY_JOB_ID}
TASK_ID=${SLURM_ARRAY_TASK_ID}

START_IDX=$(( TASK_ID * NODE_CORES ))
END_IDX=$(( START_IDX + NODE_CORES - 1 ))
if [ "$END_IDX" -ge "$n_total" ]; then
    END_IDX=$(( n_total - 1 ))
fi

run_chunk() {
    local i=$1
    local start_offset=${offsets[$i]}
    local end_offset=${offsets[$(( i + 1 ))]}
    local byte_count=$(( end_offset - start_offset ))
    local scratch_parent=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
    local task_dir
    local smiles_file
    local log_prefix
    local log_file
    local chunk_status

    printf -v log_prefix "in%04d" "$i"
    log_file="${log_prefix}.log"
    task_dir=$(mktemp -d "${scratch_parent}/msani-${ARRAY_ID}-${i}-XXXXXX")
    smiles_file="${task_dir}/${log_prefix}.smi"
    dd if=input.data of="$smiles_file" iflag=skip_bytes,count_bytes \\
       skip="$start_offset" count="$byte_count" status=none
    MSANI_PATH -i "$smiles_file" -j 1'''

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
    mv in*_clean* in/processed -f 2>/dev/null

    merge_shards() {
        local source_dir=$1
        local pattern=$2
        local destination=$3
        local temporary="${destination}.tmp"
        local shard
        local shards=()

        while IFS= read -r -d '' shard; do
            shards+=("$shard")
        done < <(find "$source_dir" -maxdepth 1 -type f -name "$pattern" -print0 | sort -zV)

        # Always produce a deterministic output, including an empty file when
        # no molecules belong to this category.
        : > "$temporary" || return 1
        for shard in "${shards[@]}"; do
            cat -- "$shard" >> "$temporary" || {
                rm -f "$temporary"
                return 1
            }
        done
        mv -f "$temporary" "$destination" || return 1

        if [ "${#shards[@]}" -gt 0 ]; then
            rm -f -- "${shards[@]}" || return 1
        fi
    }

    if ! merge_shards in/processed 'in*_clean*' in/processed/processed.smi; then
        echo "Failed to merge processed output shards; retaining staging data for recovery."
        exit 1
    fi
    if ! merge_shards in/removed 'in*_rejected*' in/removed/removed.smi; then
        echo "Failed to merge rejected output shards; retaining staging data for recovery."
        exit 1
    fi

    # Remove the shared staging data only after both ordered merges succeed.
    rm -f input.data input.offsets

    # Check for failed tasks using sacct
    failed_tasks=$(sacct -j "${ARRAY_ID}" --format='JobID%30,State' --noheader | grep 'NODE_FAIL' | awk -F_ '{print $2}' | awk '{print $1}' | tr '\n' ',' | sed 's/,$//')

    if [ -n "$failed_tasks" ]; then
        echo "Failed tasks detected: $failed_tasks"
        echo "sbatch --array=${failed_tasks} submit_msani.sh" > RESUBMIT_FAILED_JOBS.txt
        echo "Instructions for resubmitting failed jobs written to RESUBMIT_FAILED_JOBS.txt"

    else
        echo "No failed tasks detected."
    fi

    echo "Cleanup complete."
fi
"""

remove_lock_files = '''
msani_status=$?
if [ "$msani_status" -eq 0 ]; then
    rm -f "${log_prefix}.lock"
fi
'''

remove_lock_files_whole_node = '''
pids=()
for i in $(seq $START_IDX $END_IDX); do
    run_chunk "$i" &
    pids+=("$!")
done

overall_status=0
for pid in "${pids[@]}"; do
    wait "$pid" || overall_status=1
done
'''

def count_input_lines(file_path):
    """Count input lines with native tools, streaming decompression if needed."""
    compression = detect_input_compression(file_path)
    if compression and not shutil.which(_DECOMPRESSOR_COMMANDS[compression]):
        with _DECOMPRESSOR_OPEN[compression](file_path, 'rb') as source:
            return sum(block.count(b'\n') for block in iter(
                lambda: source.read(1024 * 1024), b''))
    if compression is None:
        result = subprocess.run(
            ['wc', '-l', str(file_path)],
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )
    else:
        decompressor_command = [
            _DECOMPRESSOR_COMMANDS[compression],
            '-cd',
            str(file_path),
        ]
        decompressor = subprocess.Popen(
            decompressor_command,
            stdout=subprocess.PIPE,
        )
        try:
            result = subprocess.run(
                ['wc', '-l'],
                stdin=decompressor.stdout,
                stdout=subprocess.PIPE,
                text=True,
            )
        finally:
            if decompressor.stdout is not None:
                decompressor.stdout.close()
        decompressor_returncode = decompressor.wait()
        result.check_returncode()
        if decompressor_returncode != 0:
            raise subprocess.CalledProcessError(
                decompressor_returncode,
                decompressor_command,
            )

    return int(result.stdout.split()[0])


def _write_offset_index(source, index_path, lines_per_job, destination=None):
    """Stream input data and record byte boundaries for each batch job."""
    boundaries = [0]
    total_bytes = 0
    line_count = 0
    last_newline_offset = 0

    while True:
        block = source.read(1024 * 1024)
        if not block:
            break
        if destination is not None:
            destination.write(block)

        block_start = total_bytes
        cursor = 0
        while True:
            newline = block.find(b'\n', cursor)
            if newline < 0:
                break
            line_count += 1
            last_newline_offset = block_start + newline + 1
            if line_count % lines_per_job == 0:
                boundaries.append(last_newline_offset)
            cursor = newline + 1
        total_bytes += len(block)

    # A final line without a newline still belongs to the last job.
    if total_bytes > last_newline_offset:
        line_count += 1
    if total_bytes > boundaries[-1]:
        boundaries.append(total_bytes)

    with open(index_path, 'w', encoding='ascii') as index_stream:
        for boundary in boundaries:
            index_stream.write(f'{boundary}\n')

    return len(boundaries) - 1


def prepare_batch_input(file_path, output_dir, lines_per_job):
    """Create one shared data file (or symlink) and its byte-offset index."""
    if lines_per_job <= 0:
        raise ValueError('lines_per_job must be greater than zero')

    file_path = Path(file_path).resolve()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_path = output_dir / 'input.data'
    index_path = output_dir / 'input.offsets'
    compression = detect_input_compression(file_path)

    if compression is None:
        with open(file_path, 'rb') as source:
            n_jobs = _write_offset_index(
                source,
                index_path,
                lines_per_job,
            )
        # A symlink gives every job the same folder-local path without copying
        # an already uncompressed input.
        data_path.symlink_to(file_path)
        return n_jobs

    if not shutil.which(_DECOMPRESSOR_COMMANDS[compression]):
        with _DECOMPRESSOR_OPEN[compression](file_path, 'rb') as source:
            with open(data_path, 'wb') as destination:
                return _write_offset_index(
                    source, index_path, lines_per_job, destination)

    decompressor_command = [
        _DECOMPRESSOR_COMMANDS[compression],
        '-cd',
        str(file_path),
    ]
    decompressor = subprocess.Popen(
        decompressor_command,
        stdout=subprocess.PIPE,
    )
    try:
        with open(data_path, 'wb') as destination:
            n_jobs = _write_offset_index(
                decompressor.stdout,
                index_path,
                lines_per_job,
                destination,
            )
    finally:
        if decompressor.stdout is not None:
            decompressor.stdout.close()
    decompressor_returncode = decompressor.wait()
    if decompressor_returncode != 0:
        raise subprocess.CalledProcessError(
            decompressor_returncode,
            decompressor_command,
        )
    return n_jobs


def split_input_file(file_path, output_dir, lines_per_file):
    """Split an input using native tools, streaming decompression if needed."""
    if lines_per_file <= 0:
        raise ValueError('lines_per_file must be greater than zero')

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / 'in'
    split_command = [
        'split',
        '-l', str(lines_per_file),
        '-d',
        '-a', '4',
    ]

    compression = detect_input_compression(file_path)
    if compression and not shutil.which(_DECOMPRESSOR_COMMANDS[compression]):
        chunk_paths = []
        with _DECOMPRESSOR_OPEN[compression](file_path, 'rb') as source:
            destination = None
            try:
                for index, line in enumerate(source):
                    if index % lines_per_file == 0:
                        if destination is not None:
                            destination.close()
                        chunk_path = output_dir / f'in{len(chunk_paths):04d}.smi'
                        destination = open(chunk_path, 'wb')
                        chunk_paths.append(chunk_path)
                    destination.write(line)
            finally:
                if destination is not None:
                    destination.close()
        return chunk_paths
    if compression is None:
        subprocess.run(
            split_command + [str(file_path), str(output_prefix)],
            check=True,
        )
    else:
        decompressor_command = [
            _DECOMPRESSOR_COMMANDS[compression],
            '-cd',
            str(file_path),
        ]
        decompressor = subprocess.Popen(
            decompressor_command,
            stdout=subprocess.PIPE,
        )
        try:
            split_result = subprocess.run(
                split_command + ['-', str(output_prefix)],
                stdin=decompressor.stdout,
            )
        finally:
            if decompressor.stdout is not None:
                decompressor.stdout.close()
        decompressor_returncode = decompressor.wait()
        split_result.check_returncode()
        if decompressor_returncode != 0:
            raise subprocess.CalledProcessError(
                decompressor_returncode,
                decompressor_command,
            )

    # BSD split lacks GNU's --additional-suffix option. Renaming the relatively
    # small number of chunks keeps the data path native and works on both.
    chunk_paths = []
    for path in sorted(output_dir.glob('in[0-9]*')):
        if not path.name[2:].isdigit():
            continue
        chunk_path = path.with_name(f'{path.name}.smi')
        path.rename(chunk_path)
        chunk_paths.append(chunk_path)
    return chunk_paths


def prepare_batch_chunks(file_path, output_dir, lines_per_file):
    """Create job chunks and a deterministic directory listing."""
    chunk_paths = split_input_file(file_path, output_dir, lines_per_file)
    with open(Path(output_dir) / 'dirlista', 'w', encoding='utf-8') as stream:
        for chunk_path in chunk_paths:
            stream.write(f'{chunk_path.name}\n')
    return chunk_paths

def parse_flags_single_job(args: dict, parser):
    """Parse the flags for a single job

    Args:
        args (dict): Arguments from the command line

    Returns:
        str: The flags for a single job
    """
    flags = []
    omitted_args = ["input_files", "input_list", "config", "smiles", "proj_name", "timelimit", "lines", "max_jobs", "help", "whole_node", "whole_node_cores", "partition"]
    # Load the config file to see if the defaults are really the defaults by intention or already set
    # by the config file
    config_defaults = {}
    if args.config:
        try:
            with open(args.config, 'r') as f:
                config_defaults = safe_load(f) or {}
        except Exception as e:
            print(f"Error reading YAML config: {e}")
            sys.exit(1)

    for action in parser._actions:
        arg = action.dest  # Argument name
        if not hasattr(args, arg) or arg in omitted_args:  # Skip help and internal arguments
            continue
        current_value = getattr(args, arg)  # Current value in the Namespace
        default_value = action.default      # Default value from the parser
        
        # Only skip the value is the same as the default or if the argument isn't specified
        if arg not in config_defaults and current_value == default_value:
            continue

        if isinstance(current_value, bool):
            # Boolean flags
            # Handle the counterintuitive flags
            if arg in ("taurdkit", "cleanup") and not current_value:
                flags.append(f"--no{arg}")
            # Handle the rest of the flags
            elif arg in ("stereoisomers", "neutralize") and not current_value:
                flags.append(f"--no-{arg}")
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

def write_single_job_script(slurm_header: str, slurm_script: str, defaults=None):
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

    slurm_script = slurm_script.replace('MSANI_PATH', msani_path)

    from msani.config import load_defaults
    from yaml import safe_dump
    snapshot = 'msani_defaults.yaml'
    if defaults is None:
        defaults, _ = load_defaults()
    Path(snapshot).write_text(safe_dump(defaults), encoding='utf-8')
    # SLURM copies scripts to its spool directory; the snapshot stays at submission.
    slurm_script = '\nexport MSANI_CONFIG="${SLURM_SUBMIT_DIR:-$PWD}/msani_defaults.yaml"\n' + slurm_script
    with open('submit_msani.sh', 'w') as f:
        f.write(slurm_header)
        f.write(slurm_script)

def test_batch_mode(args: dict, header: str, script: str, defaults=None):
    file = Path(args.prefix) / args.input_files[0]
    prefix = Path(args.prefix) / file.stem  # Ensure prefix is within temp_dir

    prepare_batch_input(file, prefix, args.lines)
    os.chdir(prefix)
    write_single_job_script(header, script, defaults)

def Split_Submit_jobs(args: dict, parser):
    """Split the input files into chunks and submit jobs to the cluster

    Args:
        args (dict): Arguments from the command line

    Returns:
        None
    """
    from msani.config import load_defaults
    configurations, _ = load_defaults()
    max_array_size = configurations['MAX_ARRAY_SIZE']
    max_limit_project = configurations['MAX_LIMIT_PROJECT']

    # Replace the PROJECT_NAME with the project name and time limit for SLURM
    
    cores = args.whole_node_cores

    if args.whole_node:
        # --- Whole-node header (work on a local copy) ---
        active_header = slurm_header_whole_node
        if args.proj_name is not None:
            active_header = active_header.replace('PROJECT_NAME', f'#SBATCH -A {args.proj_name}')
        else:
            active_header = active_header.replace('PROJECT_NAME', '')
        if args.partition is not None:
            active_header = active_header.replace('PARTITION_NAME', f'#SBATCH --partition={args.partition}')
        else:
            active_header = active_header.replace('PARTITION_NAME', '')
        active_header = active_header.replace('NODE_CORES', str(cores))
        active_header = active_header.replace('TIME_LIMIT', f'{args.timelimit}:00:00')
        if not args.gen3d:
            active_header = active_header.replace('msani_3d', 'msani_2d')

        # --- Whole-node execution script (local copy) ---
        # Inject NODE_CORES and append flags + background launch + wait
        flags = parse_flags_single_job(args, parser)
        active_script = slurm_script_whole_node.replace('NODE_CORES', str(cores))
        active_script += flags + ''' >> "$log_file" 2>&1
    chunk_status=$?
    rm -rf "$task_dir"
    if [ "$chunk_status" -eq 0 ]; then
        rm -f "${log_prefix}.lock"
    fi
    return "$chunk_status"
}
'''
        active_script += remove_lock_files_whole_node
        if args.cleanup:
            active_script += '\nlog_file="node_slurm_${ARRAY_ID}_${TASK_ID}.log"\n' + cleanup_script
        active_script += '\nexit "$overall_status"\n'
    else:
        # --- Standard (1-core-per-task) header (local copy) ---
        active_header = slurm_header
        if args.proj_name is not None:
            active_header = active_header.replace('PROJECT_NAME', f'#SBATCH -A {args.proj_name}')
        else:
            active_header = active_header.replace('PROJECT_NAME', '')

        if args.partition is not None:
            active_header = active_header.replace('PARTITION_NAME', f'#SBATCH --partition={args.partition}')
        else:
            active_header = active_header.replace('PARTITION_NAME', '')
        active_header = active_header.replace('TIME_LIMIT', f'{args.timelimit}:00:00')
        if not args.gen3d:
            active_header = active_header.replace('msani_3d', 'msani_2d')
        if args.lines >= 100_000: active_header = active_header.replace('MEMORY', '8G')
        elif args.lines >= 50_000: active_header = active_header.replace('MEMORY', '6G')
        elif args.lines >= 25_000: active_header = active_header.replace('MEMORY', '4G')
        else: active_header = active_header.replace('MEMORY', '2G')

        # --- Standard execution script (local copy) ---
        flags = parse_flags_single_job(args, parser)
        active_script = slurm_script + flags + ' >> "$log_file" 2>&1\n' + remove_lock_files
        if args.cleanup:
            active_script += cleanup_script
        active_script += '\nexit "$msani_status"\n'


    if args.test: 
        test_batch_mode(args, active_header, active_script, configurations)
    else:
        print(logo)
        if args.proj_name is not None: print(f"Using project name (-A): {args.proj_name}")
        if args.whole_node:
            print(f"Using whole-node mode:   partition = {args.partition}, {cores} cores per node")
        else:
            if args.partition is not None: print(f"Using partition = {args.partition}")
        print(f"Time limit for each job (-tl): {args.timelimit} hours")
        print(f"Maximum number of jobs in an array: {max_array_size} jobs")
        print(f"Maximum number of jobs running parallelly (-mj): {args.max_jobs} jobs")
        print(f"Number of compounds per job (-l): {args.lines} lines\n")

        squeue_command = ['squeue', '--noheader', '--array']
        if args.proj_name is not None:
            squeue_command.extend(['--account', str(args.proj_name)])
        result = subprocess.run(
            squeue_command,
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )
        current_running_jobs = len(result.stdout.splitlines())
        
        print("Checking input files...")
        for file in args.input_files:
            if not os.path.exists(file):
                print(f"File {file} does not exist. Please check the path and try again.")
                print(f"Exitting MolSanitizer...")
                return

        # Wait for 5 seconds before proceeding
        print("Waiting 5 seconds to review the configurations...")
        time.sleep(5)
        print('Submitting jobs...\n')
        for file in args.input_files:
            file = Path(file)
            prefix = Path(file.stem)
            replace_prefix = False
            if os.path.exists(prefix):
                remove_folder = input(f"Folder {prefix} already exists. Do you want to remove it? (y/n): ")
                if remove_folder.lower() in ['y','yes']:
                    replace_prefix = True
                else:
                    print(f"Skipping input file {file}.\n")
                    continue

            staging_dir = Path(tempfile.mkdtemp(
                prefix=f'.{prefix.name}.msani-',
                dir=prefix.parent,
            ))
            try:
                # Staging and indexing are also the counting pass. Compressed
                # inputs are therefore decompressed exactly once.
                n_jobs = prepare_batch_input(
                    file,
                    staging_dir,
                    args.lines,
                )
                n_array_tasks = (
                    math.ceil(n_jobs / cores) if args.whole_node else n_jobs
                )

                if n_jobs == 0:
                    print(f"Input file {file} is empty; skipping it.\n")
                    continue

                print(f"File {file} was indexed into {n_jobs} jobs.")
                if args.whole_node:
                    print(
                        f"This input requires {n_array_tasks} array tasks "
                        f"({cores} jobs/node)."
                    )

                # Each input is submitted as its own array, so the scheduler's
                # array-size limit applies to this input independently.
                if n_array_tasks > max_array_size:
                    print(
                        f"Skipping {file}: its array requires {n_array_tasks} "
                        f"tasks, exceeding the limit of {max_array_size}. "
                        "Please increase the number of lines per job.\n"
                    )
                    continue

                # Project capacity is shared by all arrays. Include arrays
                # submitted earlier in this invocation because they may not be
                # visible in squeue immediately.
                if current_running_jobs + n_array_tasks > max_limit_project:
                    available_tasks = max(
                        0,
                        max_limit_project - current_running_jobs,
                    )
                    print(
                        f"Skipping {file}: it requires {n_array_tasks} array "
                        f"tasks, but the project currently has capacity for "
                        f"{available_tasks}.\n"
                    )
                    continue

                if replace_prefix:
                    print(f"Removing folder {prefix}...\n")
                    shutil.rmtree(prefix)
                staging_dir.rename(prefix)

                previous_dir = Path.cwd()
                try:
                    os.chdir(prefix)
                    for job_index in range(n_jobs):
                        with open(f'in{job_index:04d}.lock', 'w') as lock:
                            lock.write('')
                    if args.whole_node:
                        print(f"Submitting {n_array_tasks} array tasks ({n_jobs} jobs, {cores} per node)\n")
                    else:
                        print(f"Submitting {n_jobs} jobs\n")
                    write_single_job_script(active_header, active_script, configurations)
                    submission = subprocess.run(
                        [
                            'sbatch',
                            f'--array=0-{n_array_tasks-1}%{args.max_jobs}',
                            'submit_msani.sh',
                        ],
                    )
                finally:
                    os.chdir(previous_dir)

                if submission.returncode == 0:
                    current_running_jobs += n_array_tasks
                else:
                    print(f"Submission failed for {file}; continuing with the next input.\n")
            finally:
                # Once renamed, staging_dir no longer exists. On validation or
                # splitting failures this removes only our private temporary
                # directory and leaves any existing result directory intact.
                if staging_dir.exists():
                    shutil.rmtree(staging_dir)
        
def main():
    from msani.config import dispatch
    if dispatch(sys.argv[1:]):
        return

    args, parser = parsers.parseArguments(sys.argv[1:], batch_mode=True)
    rdkit_version = rdBase.rdkitVersion
    if args.version:
        print(logo)
        print(f"MolSanitizer version: {__version__}")
        print(f"RDKit version: {rdkit_version}")
        return
    if (args.input_files is None or len(args.input_files) == 0):
        print("No input files provided. Please provide input files using the -i or --input_files option.")
        print("Exiting MolSanitizer...")
        return
    Split_Submit_jobs(args, parser)


if __name__=="__main__":
    main()
