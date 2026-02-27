#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# get_msl_models.sh  —  Download MSL-specific models for virus|complete
#
# Usage:
#   bash scripts/get_msl_models.sh virus  MSL40
#   bash scripts/get_msl_models.sh complete MSL39
#
# Installs:
#   tool_data/<set>/models/MSLxx/{Family,Genus}/{idx_to_label.json, model.pt}
#
# Requirements: curl, zstd, tar
# ==============================================================================

usage() {
  echo "Usage: $0 {virus|complete} MSLxx"
  echo "Example:  $0 virus MSL40"
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] '$1' not found in PATH."; exit 1; }
}

[[ $# -eq 2 ]] || usage
SET="$1"         # virus|complete
MSL="$2"         # e.g. MSL35, MSL36, ..., MSL40

[[ "${SET}" == "virus" || "${SET}" == "complete" ]] || usage
[[ "${MSL}" =~ ^MSL[0-9]+$ ]] || { echo "[ERROR] MSL must be of the form MSL35, MSL36, ..."; exit 1; }

need_cmd curl
need_cmd zstd
need_cmd tar

# Base directory (can be overridden with VPF_TOOL_DATA=/path)
TOOL_DATA_DIR="${VPF_TOOL_DATA:-tool_data}"

# Build URL depending on set and MSL
# Adjust the domain/path to match your hosting setup:
BASE_URL="${BASE_URL:-https://bioinfo.uib.es/~recerca/vpfclass2}"

if [[ "${SET}" == "virus" ]]; then
  PKG_NAME="models_virus_${MSL}.tar.zst"
  PKG_PATH="models_virus"
  DEST_DIR="${TOOL_DATA_DIR}/virus_markers/models"
else
  PKG_NAME="models_complete_${MSL}.tar.zst"
  PKG_PATH="models_complete"
  DEST_DIR="${TOOL_DATA_DIR}/complete_markers/models"
fi

URL="${BASE_URL}/${PKG_PATH}/${PKG_NAME}"

# If the destination folder for this MSL already exists, skip
if [[ -d "${DEST_DIR}/${MSL}" ]]; then
  echo "[SKIP] ${DEST_DIR}/${MSL} already exists. Nothing to do."
  exit 0
fi

echo "[INFO] Creating destination: ${DEST_DIR}"
mkdir -p "${DEST_DIR}"

echo "[INFO] Downloading: ${URL}"
TMP="$(mktemp -t vpf_msl_XXXX.tar.zst)"
curl -L --fail --retry 3 --continue-at - -o "${TMP}" "${URL}"

echo "[INFO] Extracting ${PKG_NAME} → ${DEST_DIR} (zstd -d | tar -x)…"
# The tar archive has MSLxx/ at its root, so we extract into DEST_DIR
zstd -d --long=31 -T0 < "${TMP}" | tar -C "${DEST_DIR}" -xf -

# Quick sanity check
if [[ ! -d "${DEST_DIR}/${MSL}/Genus" || ! -d "${DEST_DIR}/${MSL}/Family" ]]; then
  echo "[ERROR] Could not find ${DEST_DIR}/${MSL}/{Genus,Family} after extraction. Please verify the tar contents."
  rm -f "${TMP}"
  exit 1
fi

rm -f "${TMP}"
echo "[OK] Models ${SET} ${MSL} installed in ${DEST_DIR}/${MSL}"

