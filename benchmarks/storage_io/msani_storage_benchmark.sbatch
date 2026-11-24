#!/bin/bash
#SBATCH -J msani_storage_bench
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH -o msani_storage_bench-%A_%a.out
#
# Benchmark where msani should write while it generates conformers.
#
# For every storage location in TARGETS, on the same node and with the same
# msani flags, the job:
#   1. stages the input SMILES into a fresh working directory on that storage,
#   2. runs msani with that directory as its working directory (msani writes
#      its outputs to the current directory),
#   3. copies the outputs back to RESULTS_DIR,
# and records the wall time of each phase in RESULTS_DIR/timings.tsv.
#
# "stage" + "run" + "copyback" is the end-to-end time until the outputs sit in
# RESULTS_DIR, which is the fair comparison between running on node-local tmp
# (and copying back) and writing directly to a shared filesystem.
#
# Everything is configured through environment variables, so nothing in this
# file needs editing. Examples:
#
#   # One target
#   INPUT=/path/to/DUDEZ/AA2AR.smi \
#   SCRATCH_ROOT=/scratch/$USER MASC_ROOT=/masc/$USER \
#   sbatch -A <account> -p <partition> msani_storage_benchmark.sbatch
#
#   # Several DUDE-Z targets as an array, one target per task
#   ls /path/to/DUDEZ/*.smi > targets.txt
#   INPUT_LIST=targets.txt SCRATCH_ROOT=/scratch/$USER MASC_ROOT=/masc/$USER \
#   sbatch --array=0-$(( $(wc -l < targets.txt) - 1 )) \
#          -A <account> -p <partition> msani_storage_benchmark.sbatch
#
# Variables (defaults in brackets):
#   INPUT         SMILES file to process (required unless INPUT_LIST is set)
#   INPUT_LIST    file with one SMILES path per line; array task N uses line N+1
#   TARGETS       storage labels to test, in order ["tmp scratch masc"]
#   TMP_ROOT      node-local directory [$SLURM_TMPDIR, $TMPDIR or /tmp]
#   SCRATCH_ROOT  directory on scratch (required if "scratch" is in TARGETS)
#   MASC_ROOT     directory on masc (required if "masc" is in TARGETS)
#   RESULTS_DIR   where outputs and timings end up
#                 [$SLURM_SUBMIT_DIR/storage_bench_<jobid>]
#   MSANI_CMD     msani executable ["msani"]
#   NUMCORES      value for -j [$SLURM_CPUS_PER_TASK or 1]
#   MSANI_ARGS    msani flags other than -i and -j
#                 ["--gen3d --method rdkit --protonation --pH 7 -r 1
#                   --numconfs 600 --rmsd 0.5 --timeout 1
#                   --format db2 oeb.lib sdf"]
#   REPEATS       times to run every target [1]; the order of TARGETS is
#                 rotated on each repeat so no storage always runs first
#   KEEP_OUTPUTS  1 keeps the copied-back outputs, 0 deletes them after the
#                 size is recorded [1]
#
# Any additional label can be tested by exporting <LABEL>_ROOT, for example
# TARGETS="tmp scratch home" HOME_ROOT=/home/$USER/bench.

set -uo pipefail

job_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local$$}}
task_id=${SLURM_ARRAY_TASK_ID:-0}

if [[ -n "${INPUT_LIST:-}" ]]; then
    INPUT=$(sed -n "$(( task_id + 1 ))p" "$INPUT_LIST")
fi
if [[ -z "${INPUT:-}" || ! -f "$INPUT" ]]; then
    echo "Input SMILES file not found: '${INPUT:-}'. Set INPUT or INPUT_LIST." >&2
    exit 1
fi
INPUT=$(readlink -f "$INPUT")

TARGETS=${TARGETS:-"tmp scratch masc"}
TMP_ROOT=${TMP_ROOT:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}
MSANI_CMD=${MSANI_CMD:-msani}
NUMCORES=${NUMCORES:-${SLURM_CPUS_PER_TASK:-1}}
MSANI_ARGS=${MSANI_ARGS:-"--gen3d --method rdkit --protonation --pH 7 -r 1 --numconfs 600 --rmsd 0.5 --timeout 1 --format db2 oeb.lib sdf"}
REPEATS=${REPEATS:-1}
KEEP_OUTPUTS=${KEEP_OUTPUTS:-1}
RESULTS_DIR=${RESULTS_DIR:-${SLURM_SUBMIT_DIR:-$PWD}/storage_bench_${job_id}}

input_name=$(basename "$INPUT")
input_stem=${input_name%%.*}
n_mols=$(grep -cv '^[[:space:]]*$' "$INPUT")

mkdir -p "$RESULTS_DIR"
timings="$RESULTS_DIR/timings.tsv"
# Several array tasks append to the same file; take a lock per line.
exec 9>>"$timings.lock"
flock 9
if [[ ! -s "$timings" ]]; then
    printf 'job\ttask\tnode\tinput\tn_mols\tstorage\troot\trepeat\tstage_s\trun_s\tcopyback_s\ttotal_s\texit_code\toutput_files\toutput_bytes\n' > "$timings"
fi
flock -u 9

# Resolve the root directory of a storage label: tmp -> TMP_ROOT,
# scratch -> SCRATCH_ROOT, masc -> MASC_ROOT, foo -> FOO_ROOT.
root_for() {
    local var
    var="$(echo "$1" | tr '[:lower:]' '[:upper:]')_ROOT"
    echo "${!var:-}"
}

for label in $TARGETS; do
    root=$(root_for "$label")
    if [[ -z "$root" ]]; then
        echo "No directory for storage '$label': set $(echo "$label" | tr '[:lower:]' '[:upper:]')_ROOT." >&2
        exit 1
    fi
    mkdir -p "$root" || { echo "Cannot create $root" >&2; exit 1; }
done

now() { date +%s.%N; }
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.2f", b - a }'; }

echo "Node:        $(hostname)"
echo "Input:       $INPUT ($n_mols molecules)"
echo "Targets:     $TARGETS"
echo "Command:     $MSANI_CMD -i <staged input> -j $NUMCORES $MSANI_ARGS"
echo "Results dir: $RESULTS_DIR"
"$MSANI_CMD" --version 2>/dev/null | tail -2

read -r -a target_list <<< "$TARGETS"
read -r -a msani_args <<< "$MSANI_ARGS"
n_targets=${#target_list[@]}

workdir=""
trap '[[ -n "$workdir" ]] && rm -rf "$workdir"' EXIT

for (( rep = 1; rep <= REPEATS; rep++ )); do
    for (( k = 0; k < n_targets; k++ )); do
        label=${target_list[$(( (k + rep - 1) % n_targets ))]}
        root=$(root_for "$label")
        dest="$RESULTS_DIR/${input_stem}/${label}_rep${rep}"
        mkdir -p "$dest"

        echo
        echo "== $input_stem | $label ($root) | repeat $rep =="

        t0=$(now)
        workdir=$(mktemp -d "$root/msani-bench-${job_id}-${task_id}-${label}-XXXXXX")
        cp "$INPUT" "$workdir/$input_name"
        sync
        t1=$(now)

        ( cd "$workdir" && "$MSANI_CMD" -i "$input_name" -j "$NUMCORES" "${msani_args[@]}" ) \
            > "$dest/msani.log" 2>&1
        rc=$?
        sync
        t2=$(now)

        # Copy everything msani produced, but not the staged input itself.
        ( cd "$workdir" && rm -f "$input_name" && cp -a . "$dest/" )
        sync
        t3=$(now)

        rm -rf "$workdir"
        workdir=""

        n_files=$(find "$dest" -type f ! -name msani.log | wc -l)
        n_bytes=$(du -sb "$dest" | cut -f1)
        [[ "$KEEP_OUTPUTS" == "1" ]] || find "$dest" -mindepth 1 ! -name msani.log -delete

        stage=$(elapsed "$t0" "$t1")
        run=$(elapsed "$t1" "$t2")
        copy=$(elapsed "$t2" "$t3")
        total=$(elapsed "$t0" "$t3")
        echo "stage ${stage}s  run ${run}s  copyback ${copy}s  total ${total}s  exit $rc  files $n_files"

        flock 9
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$job_id" "$task_id" "$(hostname)" "$input_stem" "$n_mols" "$label" "$root" "$rep" \
            "$stage" "$run" "$copy" "$total" "$rc" "$n_files" "$n_bytes" >> "$timings"
        flock -u 9
    done
done

echo
echo "Timings for this task:"
awk -F'\t' -v t="$task_id" -v j="$job_id" 'NR == 1 || ($1 == j && $2 == t)' "$timings" | column -t -s $'\t'
