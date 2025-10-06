#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# setup_tool_data.sh  —  Descarga e instala tool_data de forma idempotente
#
# Uso:
#   bash scripts/setup_tool_data.sh virus
#   bash scripts/setup_tool_data.sh all
#
# Hace:
#   1) Crea tool_data/ si no existe y entra ahí
#   2) Descarga y extrae MSL_labelling/ si no existe
#   3) Según {virus|all}, descarga y extrae:
#        - virus_markers/ (vpf_data + models/MSL40/Family+Genus)
#        - complete_markers/ (vpf_data + models/MSL40/Family+Genus)
#   4) Respeta lo existente; no re-extrae si la carpeta raíz ya está
#   5) Vuelve a la raíz del repo
#
# Requisitos: curl, zstd, tar
#   - En macOS puedes instalar zstd con Homebrew:  brew install zstd
# ==============================================================================

# --------------------------
# CONFIG: URLs de descarga
# (carpetas raíz internas de los tars: MSL_labelling/, virus_markers/, complete_markers/)
# --------------------------
MSL_LABELLING_URL="${MSL_LABELLING_URL:-https://bioinfo.uib.es/~recerca/vpfclass2/MSL_labelling.tar.zst}"
VIRUS_PACKAGE_URL="${VIRUS_PACKAGE_URL:-https://bioinfo.uib.es/~recerca/vpfclass2/virus_markers_MSL40.tar.zst}"
COMPLETE_PACKAGE_URL="${COMPLETE_PACKAGE_URL:-https://bioinfo.uib.es/~recerca/vpfclass2/complete_markers_MSL40.tar.zst}"

# --------------------------
# Destino base
# --------------------------
TOOL_DATA_DIR="${VPF_TOOL_DATA:-tool_data}"

# Guardamos el dir inicial para volver al final
REPO_ROOT="$(pwd)"

usage() {
  echo "Uso: $0 {virus|all}"
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Necesito '$1' en PATH"; exit 1; }
}

# Descarga a fichero temporal y extrae vía 'zstd -d | tar -x' de forma portable.
# No re-extrae si ya existe la carpeta raíz esperada.
download_and_extract_once() {
  local url="$1"             # URL del .tar.zst
  local expected_root="$2"   # carpeta raíz que debe aparecer tras extraer (p.ej. MSL_labelling, virus_markers, complete_markers)
  local dest="$3"            # destino (normalmente "." dentro de tool_data)

  if [[ -d "${dest}/${expected_root}" ]]; then
    echo "[SKIP] ${expected_root}/ ya existe en ${dest}"
    return 0
  fi

  echo "[INFO] Descargando paquete: ${url}"
  local tmp_tar
  tmp_tar="$(mktemp -t vpf_dl_XXXX.tar.zst)"
  # Reintentos + reanudación (si el servidor lo permite)
  curl -L --fail --retry 3 --continue-at - -o "${tmp_tar}" "${url}"

  echo "[INFO] Extrayendo en ${dest} con 'zstd -d | tar -x' (modo compatible macOS/Linux)…"
  mkdir -p "${dest}"
  # Nota: --long=31 fue usado en compresión; para descompresión no es estrictamente necesario,
  # pero lo incluimos por compatibilidad con builds antiguos.
  zstd -d --long=31 -T0 < "${tmp_tar}" | tar -C "${dest}" -xf -

  # Verificación de que la carpeta raíz esperada existe tras extraer
  if [[ ! -d "${dest}/${expected_root}" ]]; then
    echo "[ERROR] Tras extraer, no encuentro ${dest}/${expected_root}. Revisa la estructura interna del tar."
    rm -f "${tmp_tar}"
    exit 1
  fi

  rm -f "${tmp_tar}"
  echo "[OK] ${expected_root}/ instalado."
}

verify_basic_layout() {
  echo "[CHECK] Comprobaciones rápidas…"

  # Taxonomía
  if ls -1 "${TOOL_DATA_DIR}/MSL_labelling"/MSL*/lineage.json >/dev/null 2>&1; then
    echo "[OK] MSL_labelling/*/lineage.json presente"
  else
    echo "[WARN] No encuentro lineage.json en MSL_labelling/*/"
  fi

  # virus_markers (si existe)
  if [[ -d "${TOOL_DATA_DIR}/virus_markers" ]]; then
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/names_virus.txt" ]] || echo "[WARN] Falta names_virus.txt"
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/profiles_virus.hmms" ]] || echo "[WARN] Falta profiles_virus.hmms"
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/vpf_to_index_V.json" ]] || echo "[WARN] Falta vpf_to_index_V.json"
    if [[ -d "${TOOL_DATA_DIR}/virus_markers/models/MSL40/Genus" && -d "${TOOL_DATA_DIR}/virus_markers/models/MSL40/Family" ]]; then
      echo "[OK] virus_markers modelos MSL40 (Family+Genus) presentes"
    fi
  fi

  # complete_markers (si existe)
  if [[ -d "${TOOL_DATA_DIR}/complete_markers" ]]; then
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/model_names.txt" ]] || echo "[WARN] Falta model_names.txt"
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/profiles.hmms" ]] || echo "[WARN] Falta profiles.hmms"
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/vpf_to_index.json" ]] || echo "[WARN] Falta vpf_to_index.json"
    if [[ -d "${TOOL_DATA_DIR}/complete_markers/models/MSL40/Genus" && -d "${TOOL_DATA_DIR}/complete_markers/models/MSL40/Family" ]]; then
      echo "[OK] complete_markers modelos MSL40 (Family+Genus) presentes"
    fi
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

echo "[INFO] Creando/usar carpeta ${TOOL_DATA_DIR}"
mkdir -p "${TOOL_DATA_DIR}"
cd "${TOOL_DATA_DIR}"

# 1) Siempre MSL_labelling
download_and_extract_once "${MSL_LABELLING_URL}" "MSL_labelling" "."

# 2) Paquete 'virus' o 'all'
if [[ "${MODE}" == "virus" ]]; then
  download_and_extract_once "${VIRUS_PACKAGE_URL}" "virus_markers" "."
else
  download_and_extract_once "${COMPLETE_PACKAGE_URL}" "complete_markers" "."
fi

# 3) Comprobación básica
verify_basic_layout

# 4) Volver a la raíz del repo
cd "${REPO_ROOT}"
echo "[DONE] tool_data listo. Volviendo a ${REPO_ROOT}"
