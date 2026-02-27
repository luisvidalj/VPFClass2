#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------
# run_prodigal_parallel.sh
# Usage:
#   run_prodigal_parallel.sh INPUT_FNA OUT_PATH N_CPUS [SEQS_PER_SHARD|AUTO] [--keep-intermediate]
#
# OUT_PATH can be:
#   - a directory, e.g. .../1_prodigal
#   - or a .faa file path, e.g. .../1_prodigal/output.faa
#
# Main output: <OUT_DIR>/output.faa (or the provided .faa path)
# Requirements: prodigal-gv, awk, (optional) GNU parallel; if parallel is not available, uses xargs -P
# ---------------------------------------------------------------------

# ========== Parse args ==========
INPUT_FNA="${1:-}"
RAW_OUT="${2:-}"
N_CPUS_RAW="${3:-}"
SEQS_PER_SHARD_ARG="${4:-AUTO}"
KEEP_INTERMEDIATE=false
if [[ "${5:-}" == "--keep-intermediate" ]]; then
  KEEP_INTERMEDIATE=true
fi

# ========== Basic checks ==========
if [[ -z "$INPUT_FNA" || -z "$RAW_OUT" || -z "$N_CPUS_RAW" ]]; then
  echo "Usage: $0 INPUT_FNA OUT_PATH N_CPUS [SEQS_PER_SHARD|AUTO] [--keep-intermediate]" >&2
  exit 2
fi
if [[ ! -f "$INPUT_FNA" ]]; then
  echo "[ERROR] INPUT_FNA does not exist: $INPUT_FNA" >&2
  exit 1
fi
if ! command -v prodigal-gv >/dev/null 2>&1; then
  echo "[ERROR] 'prodigal-gv' was not found in PATH." >&2
  exit 1
fi
if ! command -v awk >/dev/null 2>&1; then
  echo "[ERROR] 'awk' is required." >&2
  exit 1
fi
# N_CPUS must be an integer >= 1
if ! [[ "$N_CPUS_RAW" =~ ^[0-9]+$ ]] || [[ "$N_CPUS_RAW" -lt 1 ]]; then
  echo "[ERROR] N_CPUS must be an integer >= 1 (received: $N_CPUS_RAW)" >&2
  exit 1
fi
N_CPUS="$N_CPUS_RAW"

# ========== Normalize OUT_PATH (accepts directory or .faa file) ==========
# Portable abspath helper (no GNU realpath dependency)
abspath_dir() {
  # $1: path to an existing or non-existing directory
  # Returns absolute path to that directory (creates parent temporarily if needed)
  local d="$1"
  if [[ -d "$d" ]]; then
    (cd "$d" && pwd)
  else
    local parent
    parent="$(dirname "$d")"
    local base
    base="$(basename "$d")"
    mkdir -p "$parent"
    (cd "$parent" && echo "$(pwd)/$base")
  fi
}

OUT_BASE_DIR=""
FINAL_FAA=""

if [[ "$RAW_OUT" == *.faa ]]; then
  # User provided a final output file path
  OUT_BASE_DIR="$(abspath_dir "$(dirname "$RAW_OUT")")"
  FINAL_FAA="$OUT_BASE_DIR/$(basename "$RAW_OUT")"
else
  # User provided a directory
  OUT_BASE_DIR="$(abspath_dir "$RAW_OUT")"
  FINAL_FAA="$OUT_BASE_DIR/output.faa"
fi

mkdir -p "$OUT_BASE_DIR"

# Subdirectories for intermediate files
SHARDS_DIR="$OUT_BASE_DIR/shards"
SHARD_OUTS="$OUT_BASE_DIR/shard_outputs"
LOGS_DIR="$OUT_BASE_DIR/logs"
mkdir -p "$SHARDS_DIR" "$SHARD_OUTS" "$LOGS_DIR"

# Idempotency: if the final output already exists, exit
if [[ -s "$FINAL_FAA" ]]; then
  echo "[INFO] $FINAL_FAA already exists. Exiting."
  exit 0
fi

# ========== Count sequences ==========
N_SEQS=$(grep -c '^>' "$INPUT_FNA" || true)
if [[ "$N_SEQS" -eq 0 ]]; then
  echo "[ERROR] No sequences were detected in $INPUT_FNA" >&2
  exit 1
fi

# ========== AUTO sharding heuristic (by number of sequences) ==========
MAX_SHARDS=256
FLOOR_SEQS_PER_SHARD=1000

if [[ "$SEQS_PER_SHARD_ARG" == "AUTO" || -z "$SEQS_PER_SHARD_ARG" ]]; then
  # TARGET_SHARDS: at least MIN_SHARDS and ~2*N_CPUS, capped
  MIN_SHARDS=$(( N_SEQS < N_CPUS ? N_SEQS : N_CPUS ))
  TARGET_SHARDS=$(( 2 * N_CPUS ))
  if [[ $TARGET_SHARDS -lt $MIN_SHARDS ]]; then
    TARGET_SHARDS=$MIN_SHARDS
  fi
  if [[ $TARGET_SHARDS -gt $MAX_SHARDS ]]; then
    TARGET_SHARDS=$MAX_SHARDS
  fi
  # ceil(N_SEQS / TARGET_SHARDS)
  SEQS_PER_SHARD=$(( (N_SEQS + TARGET_SHARDS - 1) / TARGET_SHARDS ))
  if [[ $SEQS_PER_SHARD -lt $FLOOR_SEQS_PER_SHARD ]]; then
    SEQS_PER_SHARD=$FLOOR_SEQS_PER_SHARD
  fi
  echo "[INFO] AUTO shards: N_SEQS=$N_SEQS, N_CPUS=$N_CPUS, TARGET_SHARDS=$TARGET_SHARDS, SEQS_PER_SHARD=$SEQS_PER_SHARD"
else
  # User forces SEQS_PER_SHARD explicitly
  if ! [[ "$SEQS_PER_SHARD_ARG" =~ ^[0-9]+$ ]] || [[ "$SEQS_PER_SHARD_ARG" -lt 1 ]]; then
    echo "[ERROR] SEQS_PER_SHARD must be an integer >= 1 or 'AUTO' (received: $SEQS_PER_SHARD_ARG)" >&2
    exit 1
  fi
  SEQS_PER_SHARD="$SEQS_PER_SHARD_ARG"
  echo "[INFO] USER shards: SEQS_PER_SHARD=$SEQS_PER_SHARD"
fi

# ========== Split by sequences (do not split contigs) ==========
echo "[INFO] Splitting $INPUT_FNA into shards of $SEQS_PER_SHARD sequences..."
awk -v n="$SEQS_PER_SHARD" -v outdir="$SHARDS_DIR" '
BEGIN { seq=0; fileidx=1 }
{
  if ($0 ~ /^>/) {
    seq++
    fileidx = int((seq-1)/n)+1
    fname = sprintf("%s/shard_%04d.fna", outdir, fileidx)
  }
  print $0 >> fname
}
' "$INPUT_FNA"

# Remove empty shards just in case
find "$SHARDS_DIR" -type f -name 'shard_*.fna' -size 0 -delete || true

# Shard list
mapfile -t SHARDS < <(LC_ALL=C ls -1 "$SHARDS_DIR"/shard_*.fna 2>/dev/null || true)
N_SHARDS="${#SHARDS[@]}"
if [[ "$N_SHARDS" -eq 0 ]]; then
  echo "[ERROR] No shards were generated." >&2
  exit 1
fi

echo "[INFO] Shards generated: $N_SHARDS"
if [[ "$N_SHARDS" -eq 1 ]]; then
  echo "[WARN] Only 1 shard; parallelization will not provide much speedup."
fi

# ========== Function to process a shard ==========
run_one_shard() {
  local shard_fna="$1"
  local base
  base="$(basename "$shard_fna" .fna)"         # e.g., shard_0001
  local faa_out="$SHARD_OUTS/${base}.faa"
  local faa_tagged="$SHARD_OUTS/${base}.tagged.faa"
  local log_err="$LOGS_DIR/${base}.stderr.log"
  local log_out="$LOGS_DIR/${base}.stdout.log"

  # Run prodigal-gv
  prodigal-gv -i "$shard_fna" -a "$faa_out" -p meta -q 1>"$log_out" 2>"$log_err"

  # Tag headers for uniqueness (append |shard_XXXX)
  awk -v tag="|${base}" 'BEGIN{OFS=""}
    /^>/ { sub(/\r$/, "", $0); print $0, tag; next }
         { print }' "$faa_out" > "$faa_tagged"
}

export -f run_one_shard
export SHARD_OUTS LOGS_DIR

# ========== Parallel backend ==========
echo "[INFO] Launching prodigal-gv in parallel with up to $N_CPUS jobs..."
if command -v parallel >/dev/null 2>&1; then
  # GNU parallel
  printf "%s\n" "${SHARDS[@]}" | parallel -j "$N_CPUS" --halt now,fail=1 run_one_shard {}
else
  # Fallback to xargs -P (POSIX-ish)
  printf "%s\0" "${SHARDS[@]}" | xargs -0 -n1 -P "$N_CPUS" bash -c 'run_one_shard "$0"'
fi

# ========== Verify per-shard outputs ==========
mapfile -t TAGGED < <(LC_ALL=C ls -1 "$SHARD_OUTS"/shard_*.tagged.faa 2>/dev/null || true)
if [[ "${#TAGGED[@]}" -ne "$N_SHARDS" ]]; then
  echo "[ERROR] Number of .tagged.faa files (${#TAGGED[@]}) != number of shards ($N_SHARDS). Check logs in $LOGS_DIR" >&2
  exit 1
fi

# ========== Deterministic concatenation ==========
echo "[INFO] Concatenating shard outputs into $FINAL_FAA ..."
LC_ALL=C ls -1 "$SHARD_OUTS"/shard_*.tagged.faa | sort -V | xargs cat -- > "$FINAL_FAA"

if [[ ! -s "$FINAL_FAA" ]]; then
  echo "[ERROR] $FINAL_FAA was not generated" >&2
  exit 1
fi

echo "[INFO] Prodigal-GV (parallel) finished. Proteins: $FINAL_FAA"

# ========== Optional cleanup ==========
if [[ "$KEEP_INTERMEDIATE" == false ]]; then
  rm -rf "$SHARDS_DIR" "$SHARD_OUTS"
  # Keep logs by default; remove them too if you do not want them
fi

exit 0
