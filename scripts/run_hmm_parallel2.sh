#!/usr/bin/env bash
set -euo pipefail

# Uso: ./run_hmm_parallel.sh PROT_FAA HMMS_FILE OUT_DIR NUM_CPUS

# 1) Parámetros
PROT_FAA=${1:-}     # archivo .faa con proteínas predichas
HMMS_FILE=${2:-}    # archivo .hmms con todos los modelos
OUT_DIR=${3:-}      # carpeta para guardar splits y outputs
NUM_CPUS=${4:-}     # número total de CPUs disponibles (ej. 20)

# 2) Verificaciones básicas
if [[ -z "$PROT_FAA" || -z "$HMMS_FILE" || -z "$OUT_DIR" || -z "$NUM_CPUS" ]]; then
  echo "[ERROR] Uso: $0 PROT_FAA HMMS_FILE OUT_DIR NUM_CPUS" >&2
  exit 1
fi

if ! command -v hmmsearch >/dev/null 2>&1; then
  echo "[ERROR] No se encontró 'hmmsearch' en el PATH." >&2
  exit 1
fi

if [[ ! -s "$PROT_FAA" ]]; then
  echo "[ERROR] Archivo de proteínas no existe o está vacío: $PROT_FAA" >&2
  exit 1
fi
if [[ ! -s "$HMMS_FILE" ]]; then
  echo "[ERROR] Archivo HMM no existe o está vacío: $HMMS_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
LOG_ERR="$OUT_DIR/hmmsearch.stderr.log"
: > "$LOG_ERR"  # truncar log

# 3) Calcular número de splits (2 CPUs por split)
SPLITS=$(( NUM_CPUS / 2 ))
if [[ $SPLITS -lt 1 ]]; then SPLITS=1; fi

# 4) Contar secuencias y calcular tamaño por split
NUM_SEQS=$(grep -c "^>" "$PROT_FAA" || true)
if [[ $NUM_SEQS -eq 0 ]]; then
  echo "[WARN] No se encontraron secuencias ('>') en $PROT_FAA. No se ejecutará hmmsearch." >&2
  exit 0
fi
SEQ_PER_SPLIT=$(( (NUM_SEQS + SPLITS - 1) / SPLITS ))

echo "[INFO] Ejecutando HMMER en paralelo… (seqs=$NUM_SEQS, splits=$SPLITS, ~${SEQ_PER_SPLIT}/split)"
echo "[INFO] Salida silenciosa; errores (si los hay) en: $LOG_ERR"

# 5) Dividir el .faa en OUT_DIR/split_*.faa
#    (cada archivo comienza en un encabezado '>'; AWK rota de archivo cada n secuencias)
awk -v n="$SEQ_PER_SPLIT" -v dir="$OUT_DIR" '
  /^>/ { if (count++ % n == 0) { file = dir "/split_" ++i ".faa" } }
  { print > file }
' "$PROT_FAA"

# 6) Ejecutar hmmsearch (silencioso) en background para cada split
SECONDS=0
for file in "$OUT_DIR"/split_*.faa; do
  base="${file%.faa}"
  # --tblout a ${base}.tbl
  # Redirecciones: stdout -> /dev/null, stderr -> LOG_ERR (append)
  hmmsearch --cpu 2 --tblout "${base}.tbl" "$HMMS_FILE" "$file" >/dev/null 2>>"$LOG_ERR" &
done
wait

# 7) Reporte de duración
DUR=$SECONDS
printf "[INFO] HMMER terminado en %02d:%02d:%02d\n" $((DUR/3600)) $(((DUR%3600)/60)) $((DUR%60))

