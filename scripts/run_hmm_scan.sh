#!/usr/bin/env bash
# run_hmmscan_single.sh
#
# Usage:
#   ./run_hmmscan_single.sh PROT_FAA PROFILES_HMMS OUT_DIR NUM_CPUS
#
# Example:
#   ./run_hmmscan_single.sh proteins.faa profiles.hmms out 24
#
# Description:
#   Runs hmmscan in single-process, multi-threaded mode against a pressed
#   HMM database. This avoids splitting the FASTA file and launching
#   multiple hmmsearch jobs, improving I/O efficiency.

set -euo pipefail

# ---- 1) Input arguments ----
PROT_FAA=${1:?Missing PROT_FAA (FASTA of proteins)}
PROFILES_HMMS=${2:?Missing PROFILES_HMMS (.hmms profile database)}
OUT_DIR=${3:?Missing OUT_DIR}
NUM_CPUS=${4:?Missing NUM_CPUS (integer >0)}

mkdir -p "$OUT_DIR"

# ---- 2) Basic checks ----
command -v hmmscan >/dev/null || { echo "[ERROR] 'hmmscan' not found in PATH"; exit 1; }
command -v hmmpress >/dev/null || { echo "[ERROR] 'hmmpress' not found in PATH"; exit 1; }

[[ -s "$PROT_FAA" ]] || { echo "[ERROR] PROT_FAA not found: $PROT_FAA"; exit 1; }
[[ -s "$PROFILES_HMMS" ]] || { echo "[ERROR] PROFILES_HMMS not found: $PROFILES_HMMS"; exit 1; }

# ---- 3) Ensure hmmpress indices exist ----
need_press=0
for ext in h3f h3i h3m h3p; do
  if [[ ! -s "${PROFILES_HMMS}.${ext}" ]]; then
    need_press=1
    break
  fi
  # Re-press if .hmms file is newer than any index file
  if [[ "$PROFILES_HMMS" -nt "${PROFILES_HMMS}.${ext}" ]]; then
    need_press=1
    break
  fi
done

if [[ $need_press -eq 1 ]]; then
  echo "[INFO] Running hmmpress on ${PROFILES_HMMS}..."
  hmmpress "$PROFILES_HMMS"
fi

# ---- 4) Run hmmscan (single job, multi-threaded) ----
DOMOUT="${OUT_DIR}/all.domtbl"
TBLOUT="${OUT_DIR}/all.tbl"

echo "[INFO] Running hmmscan with --cpu ${NUM_CPUS}"
echo "[INFO] Output: ${DOMOUT} (domain table), ${TBLOUT} (target table)"

SECONDS=0

hmmscan --cpu "$NUM_CPUS" --noali \
  --domtblout "$DOMOUT" \
  --tblout "$TBLOUT" \
  "$PROFILES_HMMS" "$PROT_FAA"

DUR=$SECONDS
printf "[INFO] Total runtime: %02d:%02d:%02d (hh:mm:ss)\n" $((DUR/3600)) $(((DUR%3600)/60)) $((DUR%60))

echo "[OK] hmmscan finished successfully."
