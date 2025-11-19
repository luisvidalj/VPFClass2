#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# get_msl_models.sh  —  Descarga modelos por MSL para virus|complete
#
# Uso:
#   bash scripts/get_msl_models.sh virus  MSL40
#   bash scripts/get_msl_models.sh complete MSL39
#
# Deja:
#   tool_data/<set>/models/MSLxx/{Family,Genus}/{idx_to_label.json, model.pt}
#
# Requisitos: curl, zstd, tar
# ==============================================================================

usage() {
  echo "Uso: $0 {virus|complete} MSLxx"
  echo "Ej:  $0 virus MSL40"
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Necesito '$1' en PATH"; exit 1; }
}

[[ $# -eq 2 ]] || usage
SET="$1"         # virus|complete
MSL="$2"         # p.ej. MSL35, MSL36, ..., MSL40

[[ "${SET}" == "virus" || "${SET}" == "complete" ]] || usage
[[ "${MSL}" =~ ^MSL[0-9]+$ ]] || { echo "[ERROR] MSL debe ser tipo MSL35, MSL36, ..."; exit 1; }

need_cmd curl
need_cmd zstd
need_cmd tar

# Directorio base (puedes sobreescribir con VPF_TOOL_DATA=/ruta)
TOOL_DATA_DIR="${VPF_TOOL_DATA:-tool_data}"

# Construcción de URL en función del set y MSL
# Ajusta el dominio/ruta a tu hosting real:
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

# Si ya existe la carpeta de destino para ese MSL, salta
if [[ -d "${DEST_DIR}/${MSL}" ]]; then
  echo "[SKIP] ${DEST_DIR}/${MSL} ya existe. Nada que hacer."
  exit 0
fi

echo "[INFO] Creando destino: ${DEST_DIR}"
mkdir -p "${DEST_DIR}"

echo "[INFO] Descargando: ${URL}"
TMP="$(mktemp -t vpf_msl_XXXX.tar.zst)"
curl -L --fail --retry 3 --continue-at - -o "${TMP}" "${URL}"

echo "[INFO] Extrayendo ${PKG_NAME} → ${DEST_DIR} (zstd -d | tar -x)…"
# El tar tiene raíz MSLxx/, así que extraemos con -C DEST_DIR
zstd -d --long=31 -T0 < "${TMP}" | tar -C "${DEST_DIR}" -xf -

# Verificación rápida
if [[ ! -d "${DEST_DIR}/${MSL}/Genus" || ! -d "${DEST_DIR}/${MSL}/Family" ]]; then
  echo "[ERROR] No encuentro ${DEST_DIR}/${MSL}/{Genus,Family} tras extraer. Revisa el tar."
  rm -f "${TMP}"
  exit 1
fi

rm -f "${TMP}"
echo "[OK] Modelos ${SET} ${MSL} instalados en ${DEST_DIR}/${MSL}"
