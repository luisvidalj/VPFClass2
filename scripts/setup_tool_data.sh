#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# setup_tool_data.sh
#
# Uso:
#   bash scripts/setup_tool_data_skeleton.sh virus
#   bash scripts/setup_tool_data_skeleton.sh all
#
# Hace:
#   1) Crea tool_data/ si no existe
#   2) Descarga y extrae MSL_labelling/ (si no existe)
#   3) Según 'virus' o 'all', descarga y extrae:
#        - virus_markers/ (vpf_data + models/MSL40/{Family,Genus})
#        - complete_markers/ (vpf_data + models/MSL40/{Family,Genus})
#   4) Respeta lo ya existente (no sobreescribe si ya está)
#   5) Vuelve a la raíz del repo al terminar
#
# Requisitos: curl, tar (con --zstd), sha256sum (opcional)
# ==============================================================================

# --------------------------
# CONFIGURA AQUÍ TUS URLs
# (las rutas internas del .tar.zst deben empezar en tool_data/…)
# --------------------------
MSL_LABELLING_URL="${MSL_LABELLING_URL:-https://example.org/VPFClass2/MSL_labelling.tar.zst}"

# Paquete con: tool_data/virus_markers/{vpf_data,models/MSL40/Family,models/MSL40/Genus}
VIRUS_PACKAGE_URL="${VIRUS_PACKAGE_URL:-https://bioinfo.uib.es/~recerca/vpfclass2/virus_markers_MSL40.tar.zst}"

# Paquete con: tool_data/complete_markers/{vpf_data,models/MSL40/Family,models/MSL40/Genus}
COMPLETE_PACKAGE_URL="${COMPLETE_PACKAGE_URL:-https://example.org/VPFClass2/complete_markers_MSL40_Family_Genus.tar.zst}"

# (Opcional) archivo con checksums si lo publicas
SHA256SUMS_URL="${SHA256SUMS_URL:-}"

# --------------------------
# Destino base
# --------------------------
TOOL_DATA_DIR="${VPF_TOOL_DATA:-tool_data}"

# Guarda dir de inicio para volver luego
REPO_ROOT="$(pwd)"

usage() {
  echo "Uso: $0 {virus|all}"
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Necesito '$1' en PATH"; exit 1; }
}

download_and_extract_once() {
  local url="$1"
  local expected_root="$2"    # p.ej. "MSL_labelling" o "virus_markers" o "complete_markers"
  local dest="$3"             # p.ej. "$TOOL_DATA_DIR"

  # Si ya existe la carpeta objetivo, no hacemos nada
  if [[ -d "${dest}/${expected_root}" ]]; then
    echo "[SKIP] ${expected_root}/ ya existe en ${dest}"
    return 0
  fi

  echo "[INFO] Descargando paquete: ${url}"
  # Descarga a un tmp
  local tmp_tar
  tmp_tar="$(mktemp -t vpf_dl_XXXX.tar.zst)"
  curl -L --fail --retry 3 --continue-at - -o "${tmp_tar}" "${url}"

  echo "[INFO] Extrayendo en ${dest} (con zstd)…"
  mkdir -p "${dest}"
  tar --zstd -xf "${tmp_tar}" -C "${dest}"

  # Verificamos que lo esperado apareció
  if [[ ! -d "${dest}/${expected_root}" ]]; then
    echo "[ERROR] Tras extraer, no encuentro ${dest}/${expected_root}. ¿El tar tiene rutas correctas (tool_data/${expected_root})?"
    rm -f "${tmp_tar}"
    exit 1
  fi

  rm -f "${tmp_tar}"
  echo "[OK] ${expected_root}/ instalado."
}

verify_basic_layout() {
  echo "[CHECK] Comprobaciones rápidas…"

  # Taxonomía (requerida)
  if ! ls -1 "${TOOL_DATA_DIR}/MSL_labelling"/MSL*/lineage.json >/dev/null 2>&1; then
    echo "[WARN] No encuentro lineage.json en MSL_labelling/*/"
  else
    echo "[OK] MSL_labelling/*/lineage.json presente"
  fi

  # virus_markers (si existe)
  if [[ -d "${TOOL_DATA_DIR}/virus_markers" ]]; then
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/names_virus.txt" ]] || echo "[WARN] Falta names_virus.txt"
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/profiles_virus.hmms" ]] || echo "[WARN] Falta profiles_virus.hmms"
    [[ -f "${TOOL_DATA_DIR}/virus_markers/vpf_data/vpf_to_index_V.json" ]] || echo "[WARN] Falta vpf_to_index_V.json"
    [[ -d "${TOOL_DATA_DIR}/virus_markers/models/MSL40/Genus" ]] && [[ -d "${TOOL_DATA_DIR}/virus_markers/models/MSL40/Family" ]] \
      && echo "[OK] virus_markers modelos MSL40 (Family+Genus) presentes"
  fi

  # complete_markers (si existe)
  if [[ -d "${TOOL_DATA_DIR}/complete_markers" ]]; then
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/model_names.txt" ]] || echo "[WARN] Falta model_names.txt"
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/profiles.hmms" ]] || echo "[WARN] Falta profiles.hmms"
    [[ -f "${TOOL_DATA_DIR}/complete_markers/vpf_data/vpf_to_index.json" ]] || echo "[WARN] Falta vpf_to_index.json"
    [[ -d "${TOOL_DATA_DIR}/complete_markers/models/MSL40/Genus" ]] && [[ -d "${TOOL_DATA_DIR}/complete_markers/models/MSL40/Family" ]] \
      && echo "[OK] complete_markers modelos MSL40 (Family+Genus) presentes"
  fi
}

# --------------------------
# MAIN
# --------------------------
[[ $# -eq 1 ]] || usage
MODE="$1"
[[ "${MODE}" == "virus" || "${MODE}" == "all" ]] || usage

need_cmd curl
need_cmd tar

echo "[INFO] Creando/usar carpeta ${TOOL_DATA_DIR}"
mkdir -p "${TOOL_DATA_DIR}"

# Movernos a tool_data/ para extraer paquetes con rutas tool_data/… coherentes
cd "${TOOL_DATA_DIR}"

# 1) Siempre MSL_labelling
download_and_extract_once "${MSL_LABELLING_URL}" "MSL_labelling" "."

# 2) Paquete 'virus' o 'all'
if [[ "${MODE}" == "virus" ]]; then
  download_and_extract_once "${VIRUS_PACKAGE_URL}" "virus_markers" "."
else
  download_and_extract_once "${COMPLETE_PACKAGE_URL}" "complete_markers" "."
fi

# (Opcional) Verificación de checksums globales si se publica SHA256SUMS
if [[ -n "${SHA256SUMS_URL}" ]]; then
  echo "[INFO] Verificando SHA256SUMS…"
  tmp_sum="$(mktemp -t vpf_sha_XXXX.txt)"
  curl -L --fail -o "${tmp_sum}" "${SHA256SUMS_URL}" || { echo "[WARN] No pude descargar SHA256SUMS"; true; }
  if [[ -s "${tmp_sum}" ]]; then
    # Nota: esto verifica los archivos existentes; omite los que no estén
    sha256sum -c "${tmp_sum}" || echo "[WARN] Algún checksum no coincide. Revisa tus descargas."
  fi
  rm -f "${tmp_sum}"
fi

# 3) Comprobación básica de layout
verify_basic_layout

# 4) Volver a la raíz del repo
cd "${REPO_ROOT}"
echo "[DONE] tool_data listo. Volviendo a ${REPO_ROOT}"
