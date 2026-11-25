"""
MolSanitizer in the batch mode.
"""

from msani import __version__


import argparse
import os
import re
import shlex
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
MEMORY_LIMIT
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
    # A missing lock means this chunk already succeeded (resubmitted task).
    [ -f "${log_prefix}.lock" ] || return 0
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

    mkdir -p in/processed in/removed log
    mv *.log log
    mv in*_rejected* in/removed -f 2>/dev/null
    mv in*_metal_error* in/removed -f 2>/dev/null
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
    if ! merge_shards in/removed 'in*_metal_error*' in/removed/metal_error.smi; then
        echo "Failed to merge metal error shards; retaining staging data for recovery."
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

# Below this default, one-core jobs are at risk of being killed for memory.
LOW_MEMORY_MB = 2000


def _default_memory_setting(text):
    """Read DefMemPerCPU/DefMemPerNode from scontrol output.

    Returns the value in MB, math.inf for UNLIMITED or 0 (which SLURM also
    treats as no limit), or None when absent. A one-core job gets the same
    amount under either setting.
    """
    match = re.search(r'\bDefMemPer(?:CPU|Node)\s*=\s*(\S+)', text)
    if match is None:
        return None
    value = match.group(1)
    if value.upper() == 'UNLIMITED':
        return math.inf
    if not value.isdigit():
        return None
    return int(value) or math.inf


def _scontrol(*arguments):
    return subprocess.run(
        ['scontrol', *arguments],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=True,
        timeout=10,
    ).stdout


def default_job_memory_mb(partition=None):
    """Memory a one-core job gets when submitted without --mem.

    Returns MB, math.inf when the cluster sets no limit, or None when it
    cannot be determined. This is advisory, so it never raises.
    """
    try:
        if partition is None:
            # One line per partition; the default one is used when none is named.
            lines = _scontrol('show', 'partition', '--oneliner').splitlines()
            partitions = [line for line in lines if re.search(r'\bDefault=YES\b', line)]
        else:
            partitions = [_scontrol('show', 'partition', name, '--oneliner')
                          for name in partition.split(',') if name]
        # A partition without its own default reports UNLIMITED and inherits
        # the cluster-wide setting.
        limits = [_default_memory_setting(text) for text in partitions]
        if any(limit in (None, math.inf) for limit in limits) or not limits:
            limits.append(_default_memory_setting(_scontrol('show', 'config')))
        limits = [limit for limit in limits if limit is not None]
        # With several partitions the job may land on the smallest default.
        finite = [limit for limit in limits if limit != math.inf]
        if finite:
            return min(finite)
        return math.inf if limits else None
    except (OSError, subprocess.SubprocessError):
        return None


# msani reads its input in chunks of this many records (cli.READ_CHUNK_SIZE).
LARGE_3D_JOB_LINES = 10_000
# Formats collected into one archive per job, which is rewritten after every chunk.
ARCHIVE_FORMATS = ('db2.tgz', 'oeb.lib')


def confirm_large_3d_jobs(args):
    """Ask before submitting 3D archive jobs that span more than one read chunk.

    Only the db2.tgz and oeb.lib formats are affected; the other formats
    write one file per molecule. Returns True when the submission should go
    ahead.
    """
    if not args.gen3d or args.lines <= LARGE_3D_JOB_LINES:
        return True
    formats = [args.format] if isinstance(args.format, str) else (args.format or [])
    archives = [name for name in ARCHIVE_FORMATS if name in formats]
    if not archives:
        return True
    archive_text = (f"the {archives[0]} archive" if len(archives) == 1
                    else f"the {' and '.join(archives)} archives")
    print(f"Warning: -l {args.lines} puts more than {LARGE_3D_JOB_LINES:,} molecules in each 3D job. "
          f"Each job works through its input in chunks of {LARGE_3D_JOB_LINES:,} molecules and rewrites "
          f"{archive_text} after every chunk, which slows the job down, "
          "and a failed or timed-out job loses more work. "
          f"Prefer many small jobs: keep -l at or below {LARGE_3D_JOB_LINES:,}; "
          "we recommend fewer than 5,000 lines per job.")
    if getattr(args, 'yes', False):
        print("Submitting anyway because --yes was given.\n")
        return True
    try:
        answer = input("Do you want to submit anyway? (y/n): ")
    except EOFError:
        # No terminal to answer on, e.g. when called from a script.
        answer = ''
    return answer.strip().lower() in ['y', 'yes']


def parse_flags_single_job(args: dict, parser):
    """Parse the flags for a single job

    Args:
        args (dict): Arguments from the command line

    Returns:
        str: The flags for a single job
    """
    flags = []
    omitted_args = ["input_files", "input_list", "config", "smiles", "proj_name", "timelimit", "lines", "max_jobs", "help", "whole_node", "whole_node_cores", "partition", "mem", "metal_error_file", "retry", "yes"]
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
                # The option string, not the dest: --embed-only is stored as embed_only.
                flags.append(action.option_strings[0] if action.option_strings else f"--{arg}")
        elif isinstance(current_value, list):
            # List arguments
            value_str = " ".join(shlex.quote(str(value)) for value in current_value)
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
    slurm_script = ('\nexport MSANI_CONFIG="${SLURM_SUBMIT_DIR:-$PWD}/msani_defaults.yaml"\n'
                    'export MSANI_BATCH_WORKER=1\n' + slurm_script)
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
        # Without --mem SLURM grants the per-core default, which is the most
        # memory a one-core job can use without being billed for extra cores.
        if args.mem is not None:
            active_header = active_header.replace('MEMORY_LIMIT', f'#SBATCH --mem={args.mem}')
        else:
            active_header = active_header.replace('MEMORY_LIMIT\n', '')

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
            if args.mem is not None: print("Warning: --mem is ignored with --whole_node.")
        else:
            if args.partition is not None: print(f"Using partition = {args.partition}")
            if args.mem is not None:
                print(f"Memory per job (--mem): {args.mem}")
            else:
                default_memory = default_job_memory_mb(args.partition)
                if default_memory is not None and default_memory < LOW_MEMORY_MB:
                    print(f"Warning: this cluster gives one-core jobs only {default_memory} MB of memory by default. "
                          "Set --mem (or SLURM_MEMORY in the msani configuration) if jobs are killed for memory.")
        print(f"Time limit for each job (-tl): {args.timelimit} hours")
        print(f"Maximum number of jobs in an array: {max_array_size} jobs")
        print(f"Maximum number of jobs running parallelly (-mj): {args.max_jobs} jobs")
        print(f"Number of compounds per job (-l): {args.lines} lines\n")

        if not confirm_large_3d_jobs(args):
            print("Nothing was submitted. Exiting MolSanitizer...")
            return

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
        
def _array_spec(task_ids):
    """Compress task IDs into a SLURM array specification, e.g. 1,4-6,9."""
    ranges = []
    for task_id in sorted(set(task_ids)):
        if ranges and task_id == ranges[-1][1] + 1:
            ranges[-1][1] = task_id
        else:
            ranges.append([task_id, task_id])
    return ','.join(str(first) if first == last else f'{first}-{last}'
                    for first, last in ranges)


def explicit_options(parser, argv, config_path=None):
    """Destinations of the options the user set on this command line.

    An option counts when it was typed or appears in the per-job YAML given
    with --config. A value that merely equals the current default cannot be
    told from "not given" by looking at the parsed arguments, hence this.
    """
    saved_actions = [(action, action.default) for action in parser._actions]
    saved_defaults = dict(parser._defaults)
    try:
        for action in parser._actions:
            action.default = argparse.SUPPRESS
        parser._defaults.clear()
        typed, _ = parser.parse_known_args(argv)
    finally:
        for action, default in saved_actions:
            action.default = default
        parser._defaults.update(saved_defaults)
    explicit = set(vars(typed))

    if config_path:
        try:
            with open(config_path, 'r') as stream:
                config_keys = set(safe_load(stream) or {})
        except Exception:
            config_keys = set()
        for action in parser._actions:
            if action.dest in config_keys or any(
                    option.lstrip('-') in config_keys for option in action.option_strings):
                explicit.add(action.dest)
    return explicit


# Options that change how --retry resubmits; every other one is baked into
# the folder's submit_msani.sh.
RETRY_OPTIONS = {'retry', 'timelimit', 'mem', 'max_jobs', 'proj_name',
                 'partition', 'test', 'config', 'yes'}


def active_array_tasks(batch_dir):
    """Array task IDs that SLURM still has queued or running for a batch folder.

    Only the user's own jobs are listed. They are matched by their working
    directory, which is the batch folder for everything msani_batch submits.
    Raises OSError or
    subprocess.SubprocessError when squeue cannot be queried.
    """
    batch_dir = os.path.realpath(batch_dir)
    # --array lists every array element on its own line, pending ones included.
    result = subprocess.run(
        ['squeue', '--me', '--noheader', '--array', '--format=%i|%Z'],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=True,
    )
    active = set()
    for line in result.stdout.splitlines():
        job_id, separator, work_dir = line.strip().partition('|')
        if not separator or os.path.realpath(work_dir) != batch_dir:
            continue
        match = re.fullmatch(r'\d+_(\d+)', job_id)
        if match is None:
            # An array element that cannot be told apart, or a job that is not
            # part of an array: nothing in this folder is safe to resubmit.
            raise subprocess.SubprocessError(
                f'cannot tell which tasks job {job_id} covers')
        active.add(int(match.group(1)))
    return active


def retry_batch_dir(batch_dir, args, explicit=frozenset()):
    """Resubmit the jobs of one batch folder that still hold a lock file.

    A lock file is removed only when its job succeeds, so the remaining ones
    mark jobs that failed, timed out, or never ran. Jobs SLURM still has
    queued or running are left alone. Returns the number of resubmitted jobs.
    """
    batch_dir = Path(batch_dir)
    script = batch_dir / 'submit_msani.sh'
    offsets = batch_dir / 'input.offsets'
    if not batch_dir.is_dir():
        print(f"Skipping {batch_dir}: it is not a folder.\n")
        return 0
    if not script.is_file():
        print(f"Skipping {batch_dir}: no submit_msani.sh, so it was not created by msani_batch.\n")
        return 0

    locked = sorted(
        int(match.group(1))
        for path in batch_dir.iterdir()
        if (match := re.fullmatch(r'in(\d+)\.lock', path.name))
    )
    if not locked:
        print(f"{batch_dir}: no lock files left, all jobs finished. Nothing to retry.\n")
        return 0
    if not offsets.is_file() or not (batch_dir / 'input.data').exists():
        print(f"Skipping {batch_dir}: the staged input (input.data, input.offsets) is missing. "
              "Prepare a new batch from the original input instead.\n")
        return 0
    n_jobs = len(offsets.read_text().split()) - 1
    unknown = [job for job in locked if job >= n_jobs]
    if unknown:
        print(f"Skipping {batch_dir}: lock files in{unknown[0]:04d}.lock and above "
              f"do not match the {n_jobs} jobs of this batch.\n")
        return 0

    # In whole-node mode one array task runs several jobs.
    whole_node = re.search(r'START_IDX=\$\(\( TASK_ID \* (\d+) \)\)', script.read_text())
    jobs_per_task = int(whole_node.group(1)) if whole_node else 1

    try:
        active = active_array_tasks(batch_dir)
    except (OSError, subprocess.SubprocessError) as error:
        print(f"Skipping {batch_dir}: could not check with squeue which jobs are still "
              f"queued or running ({error}).\n")
        return 0

    tasks = sorted({job // jobs_per_task for job in locked})
    still_active = [task for task in tasks if task in active]
    to_submit = [task for task in tasks if task not in active]
    n_retry = sum(1 for job in locked if job // jobs_per_task not in active)

    print(f"{batch_dir}: {len(locked)} of {n_jobs} jobs are unfinished.")
    if still_active:
        print(f"  {len(locked) - n_retry} of them are still queued or running in SLURM "
              f"(array tasks {_array_spec(still_active)}) and are left alone.")
    if not to_submit:
        print("  Nothing to retry.\n")
        return 0

    command = ['sbatch', f'--array={_array_spec(to_submit)}%{args.max_jobs}']
    # The folder's script keeps its original settings. Options given on this
    # command line replace them for the retried jobs, as sbatch options take
    # precedence over the #SBATCH lines of the script.
    if 'proj_name' in explicit and args.proj_name is not None:
        command.append(f'--account={args.proj_name}')
    if 'partition' in explicit and args.partition is not None:
        command.append(f'--partition={args.partition}')
    if 'timelimit' in explicit:
        command.append(f'--time={args.timelimit}:00:00')
    # Memory also applies when it comes from SLURM_MEMORY in the configuration:
    # setting it there is a request for every one-core job.
    if args.mem is not None:
        if not whole_node:
            command.append(f'--mem={args.mem}')
        elif 'mem' in explicit:
            print("  Warning: --mem is ignored for whole-node batches.")
    command.append(script.name)

    if whole_node:
        print(f"  Resubmitting {n_retry} jobs as {len(to_submit)} array tasks "
              f"({jobs_per_task} jobs per node): {' '.join(command)}")
    else:
        print(f"  Resubmitting {n_retry} jobs: {' '.join(command)}")
    if args.test:
        print("  Test mode: nothing was submitted.\n")
        return 0
    submission = subprocess.run(command, cwd=batch_dir)
    if submission.returncode != 0:
        print(f"  Submission failed for {batch_dir}.\n")
        return 0
    print()
    return n_retry


def retry_jobs(args, parser, argv=None):
    """Resubmit the unfinished jobs of batch folders created by msani_batch.

    argv is the command line the arguments were parsed from (default:
    sys.argv[1:]); it tells which options were actually given.
    """
    print(logo)
    if argv is None:
        argv = sys.argv[1:]
    explicit = explicit_options(parser, argv, args.config)
    ignored = [action.option_strings[-1] if len(action.option_strings[0]) == 2
               else action.option_strings[0]
               for action in parser._actions
               if action.dest in explicit - RETRY_OPTIONS and action.option_strings]
    if ignored:
        print("Note: retried jobs keep the options of their folder's submit_msani.sh, "
              f"so these are ignored: {', '.join(ignored)}.\n")
    for option, dest in (('--projectName', 'proj_name'), ('--partition', 'partition')):
        if dest in explicit and getattr(args, dest) is None:
            print(f"Note: {option} cannot be unset for a retry; the folder's own setting is kept.\n")
    total = sum(retry_batch_dir(batch_dir, args, explicit) for batch_dir in args.retry)
    print(f"Resubmitted {total} jobs in total.")


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
    if args.retry:
        retry_jobs(args, parser, sys.argv[1:])
        return
    if (args.input_files is None or len(args.input_files) == 0):
        print("No input files provided. Please provide input files using the -i or --input_files option.")
        print("Exiting MolSanitizer...")
        return
    from msani.io.loggers import protonation_notice
    parsers.validate_protonation_options(args)
    protonation_notice(args, batch_submission=True)
    Split_Submit_jobs(args, parser)


if __name__=="__main__":
    main()
