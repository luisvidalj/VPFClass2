#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# setup_tool_data.sh  —  Download and install tool_data in an idempotent way
#
# Usage:
#   bash scripts/setup_tool_data.sh virus
#   bash scripts/setup_tool_data.sh all
#
# What it does:
#   1) Creates tool_data/ if it does not exist and cd's into it
#   2) Downloads and extracts MSL_labelling/ if not already present
#   3) Depending on {virus|all}, downloads and extracts:
#        - virus_markers/ (vpf_data + models/MSL40/Family+Genus)
#        - complete_markers/ (vpf_data + models/MSL40/Family+Genus)
#   4) Respects existing files; does not re-extract if the expected root folder already exists
#   5) Returns to the repository root
#
# Requirements: curl, zstd, tar
#   - On macOS you can install zstd with Homebrew:  brew install zstd
# ==============================================================================

# --------------------------
# CONFIG: Download URLs
# (expected root folders inside the archives: MSL_labelling/, virus_markers/, complete_markers/)
# --------------------------
MSL_LABELLING_URL="${MSL_LABELLING_URL:-https://bioinfo.uib.es/~recerca/vpfclass2/MSL_labelling.tar.zst}"
VIRUS_PACKAGE_URL="${VIRUS_PACKAGE_URL:-https://bioinfo.uib.es/~recerca/vpfclass2/virus_markers_MSL40.tar.zst}"
COMPLETE_PACKAGE_URL="${COMPLETE_PACKAGE_URL:-https://bioinfo.uib.es/~recerca/vpfclass2/complete_markers_MSL40.tar.zst}"

# --------------------------
# Base destination
# --------------------------
TOOL_DATA_DIR="${VPF_TOOL_DATA:-tool_data}"

# Save starting directory to return at the end
REPO_ROOT="$(pwd)"

usage() {
  echo "Usage: $0 {virus|all}"
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] '$1' not found in PATH."; exit 1; }
}

# Download to a temporary file and extract via portable 'zstd -d | tar -x'.
# Does not re-extract if the expected root folder already exists.
download_and_extract_once() {
  local url="$1"             # URL of the .tar.zst
  local expected_root="$2"   # root folder expected after extraction (e.g., MSL_labelling, virus_markers, complete_markers)
  local dest="$3"            # destination directory (usually "." inside tool_data)

  if [[ -d "${dest}/${expected_root}" ]]; then
    echo "[SKIP] ${expected_root}/ already exists in ${dest}"
    return 0
  fi

  echo "[INFO] Downloading package: ${url}"
  local tmp_tar
  tmp_tar="$(mktemp -t vpf_dl_XXXX.tar.zst)"
  # Retries + resume (if supported by the server)
  curl -L --fail --retry 3 --continue-at - -o "${tmp_tar}" "${url}"

  echo "[INFO] Extracting into ${dest} with 'zstd -d | tar -x' (macOS/Linux compatible)…"
  mkdir -p "${dest}"
  # Note: --long=31 was used during compression; not strictly required for decompression,
  # but we include it for compatibility with older builds.
  zstd -d --long=31 -T0 < "${tmp_tar}" | tar -C "${dest}" -xf -

  # Verify the expected root folder exists after extraction
  if [[ ! -d "${dest}/${expected_root}" ]]; then
    echo "[ERROR] After extraction, ${dest}/${expected_root} was not found. Please verify the tar internal structure."
    rm -f "${tmp_tar}"
    exit 1
  fi

  rm -f "${tmp_tar}"
  echo "[OK] ${expected_root}/ installed."
}

verify_basic_layout() {
  echo "[CHECK] Quick checks…"

  # Taxonomy
  if ls -1 "${TOOL_DATA_DIR}/MSL_labelling"/MSL*/lineage.json >/dev/null 2>&1; then
    echo "[OK] MSL_labelling/*/lineage.json present"
  else
    echo "[WARN] lineage.json not found under MSL_labelling/*/"
  fi

  # virus_markers (if present)
  if [[ -d "${TOOL_DATA_DIR}/virus_markers" ]]; then
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/names_virus.txt" ]] || echo "[WARN] Missing names_virus.txt"
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/profiles_virus.hmms" ]] || echo "[WARN] Missing profiles_virus.hmms"
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/vpf_to_index_V.json" ]] || echo "[WARN] Missing vpf_to_index_V.json"
    if [[ -d "${TOOL_DATA_DIR}/virus_markers/models/MSL40/Genus" && -d "${TOOL_DATA_DIR}/virus_markers/models/MSL40/Family" ]]; then
      echo "[OK] virus_markers MSL40 models (Family+Genus) present"
    fi
  fi

  # complete_markers (if present)
  if [[ -d "${TOOL_DATA_DIR}/complete_markers" ]]; then
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/model_names.txt" ]] || echo "[WARN] Missing model_names.txt"
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/profiles.hmms" ]] || echo "[WARN] Missing profiles.hmms"
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/vpf_to_index.json" ]] || echo "[WARN] Missing vpf_to_index.json"
    if [[ -d "${TOOL_DATA_DIR}/complete_markers/models/MSL40/Genus" && -d "${TOOL_DATA_DIR}/complete_markers/models/MSL40/Family" ]]; then
      echo "[OK] complete_markers MSL40 models (Family+Genus) present"
    fi
  fi
}

# Install non-editable mode
install_pkg() {
  # Allow skipping installation (e.g., in CI) with SKIP_PIP=1
  if [[ "${SKIP_PIP:-0}" == "1" ]]; then
    echo "[SKIP] SKIP_PIP=1 → skipping 'pip install .'"
    return 0
  fi

  echo "[INFO] Installing VPF-Class 2 (non-editable)…"
  # Use the python from the active environment
  if python -m pip install -q .; then
    echo "[OK] Package installed."
  else
    echo "[ERROR] 'pip install .' failed. Make sure your conda env is ACTIVE and you have write permissions."
    exit 1
  fi
}


# --- Helper to install the package in editable mode -----------------------
install_editable_pkg() {
  # If SKIP_PIP=1, do not install
  if [[ "${SKIP_PIP:-0}" == "1" ]]; then
    echo "[SKIP] SKIP_PIP=1 → skipping 'pip install -e .'"
    return 0
  fi

  # Check if already importable
  if python - <<'PY' >/dev/null 2>&1
import importlib.util
import sys
sys.exit(0 if importlib.util.find_spec('vpf_classifier') else 1)
PY
  then
    echo "[OK] Package 'vpf_classifier' is already importable; installation is not required."
    return 0
  fi

  echo "[INFO] Installing VPF-Class 2 in editable mode…"
  # -q for quiet; remove -q for more logs
  if python -m pip install -q -e .; then
    echo "[OK] Editable installation completed."
  else
    echo "[ERROR] 'pip install -e .' failed. Make sure your conda env is ACTIVE and you have write permissions."
    exit 1
  fi
}


# --------------------------
# MAIN
# --------------------------
[[ $# -eq 1 ]] || usage
MODE="$1"
[[ "${MODE}" == "virus" || "${MODE}" == "all" ]] || usage

need_cmd curl
need_cmd zstd
need_cmd tar

echo "[INFO] Creating/using folder ${TOOL_DATA_DIR}"
mkdir -p "${TOOL_DATA_DIR}"
cd "${TOOL_DATA_DIR}"

# 1) Always install MSL_labelling
download_and_extract_once "${MSL_LABELLING_URL}" "MSL_labelling" "."

# 2) 'virus' or 'all' package
if [[ "${MODE}" == "virus" ]]; then
  download_and_extract_once "${VIRUS_PACKAGE_URL}" "virus_markers" "."
else
  download_and_extract_once "${COMPLETE_PACKAGE_URL}" "complete_markers" "."
fi

# 3) Basic checks
verify_basic_layout

# 4) Return to repo root
cd "${REPO_ROOT}"
echo "[DONE] tool_data is ready. Returning to ${REPO_ROOT}"

#install_editable_pkg
install_pkg
