#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------
# run_prodigal_parallel.sh
# Uso:
#   run_prodigal_parallel.sh INPUT_FNA OUT_PATH N_CPUS [SEQS_PER_SHARD|AUTO] [--keep-intermediate]
#
# OUT_PATH puede ser:
#   - una carpeta, p.ej. .../1_prodigal
#   - o una ruta a archivo .faa, p.ej. .../1_prodigal/output.faa
#
# Salida principal: <OUT_DIR>/output.faa (o la ruta .faa indicada)
# Requisitos: prodigal-gv, awk, (opcional) GNU parallel; si no hay parallel, usa xargs -P
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

# ========== Checks básicos ==========
if [[ -z "$INPUT_FNA" || -z "$RAW_OUT" || -z "$N_CPUS_RAW" ]]; then
  echo "Uso: $0 INPUT_FNA OUT_PATH N_CPUS [SEQS_PER_SHARD|AUTO] [--keep-intermediate]" >&2
  exit 2
fi
if [[ ! -f "$INPUT_FNA" ]]; then
  echo "[ERROR] No existe INPUT_FNA: $INPUT_FNA" >&2
  exit 1
fi
if ! command -v prodigal-gv >/dev/null 2>&1; then
  echo "[ERROR] No se encontró 'prodigal-gv' en PATH." >&2
  exit 1
fi
if ! command -v awk >/dev/null 2>&1; then
  echo "[ERROR] 'awk' es requerido." >&2
  exit 1
fi
# N_CPUS debe ser entero >=1
if ! [[ "$N_CPUS_RAW" =~ ^[0-9]+$ ]] || [[ "$N_CPUS_RAW" -lt 1 ]]; then
  echo "[ERROR] N_CPUS debe ser un entero >= 1 (recibido: $N_CPUS_RAW)" >&2
  exit 1
fi
N_CPUS="$N_CPUS_RAW"

# ========== Normalización OUT_PATH (acepta carpeta o archivo .faa) ==========
# Función abspath portable (sin depender de realpath GNU)
abspath_dir() {
  # $1: path a directorio existente o no
  # Devuelve ruta absoluta del directorio (crea temporalmente si hace falta)
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
  # Usuario pasó una ruta a archivo final
  OUT_BASE_DIR="$(abspath_dir "$(dirname "$RAW_OUT")")"
  FINAL_FAA="$OUT_BASE_DIR/$(basename "$RAW_OUT")"
else
  # Usuario pasó carpeta
  OUT_BASE_DIR="$(abspath_dir "$RAW_OUT")"
  FINAL_FAA="$OUT_BASE_DIR/output.faa"
fi

mkdir -p "$OUT_BASE_DIR"

# Subcarpetas para temporales
SHARDS_DIR="$OUT_BASE_DIR/shards"
SHARD_OUTS="$OUT_BASE_DIR/shard_outputs"
LOGS_DIR="$OUT_BASE_DIR/logs"
mkdir -p "$SHARDS_DIR" "$SHARD_OUTS" "$LOGS_DIR"

# Idempotencia: si ya existe el resultado final, salimos
if [[ -s "$FINAL_FAA" ]]; then
  echo "[INFO] Ya existe $FINAL_FAA. Saliendo."
  exit 0
fi

# ========== Conteo de secuencias ==========
N_SEQS=$(grep -c '^>' "$INPUT_FNA" || true)
if [[ "$N_SEQS" -eq 0 ]]; then
  echo "[ERROR] No se han detectado secuencias en $INPUT_FNA" >&2
  exit 1
fi

# ========== Heurística AUTO para sharding (por nº de secuencias) ==========
MAX_SHARDS=256
FLOOR_SEQS_PER_SHARD=1000

if [[ "$SEQS_PER_SHARD_ARG" == "AUTO" || -z "$SEQS_PER_SHARD_ARG" ]]; then
  # TARGET_SHARDS: al menos MIN_SHARDS y aprox 2*N_CPUS, con tope
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
  # Usuario fuerza SEQS_PER_SHARD explícito
  if ! [[ "$SEQS_PER_SHARD_ARG" =~ ^[0-9]+$ ]] || [[ "$SEQS_PER_SHARD_ARG" -lt 1 ]]; then
    echo "[ERROR] SEQS_PER_SHARD debe ser un entero >=1 o 'AUTO' (recibido: $SEQS_PER_SHARD_ARG)" >&2
    exit 1
  fi
  SEQS_PER_SHARD="$SEQS_PER_SHARD_ARG"
  echo "[INFO] USER shards: SEQS_PER_SHARD=$SEQS_PER_SHARD"
fi

# ========== Split por secuencias (no corta contigs) ==========
echo "[INFO] Dividiendo $INPUT_FNA en shards de $SEQS_PER_SHARD secuencias..."
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

# Elimina shards vacíos por si acaso
find "$SHARDS_DIR" -type f -name 'shard_*.fna' -size 0 -delete || true

# Lista de shards
mapfile -t SHARDS < <(LC_ALL=C ls -1 "$SHARDS_DIR"/shard_*.fna 2>/dev/null || true)
N_SHARDS="${#SHARDS[@]}"
if [[ "$N_SHARDS" -eq 0 ]]; then
  echo "[ERROR] No se generaron shards." >&2
  exit 1
fi

echo "[INFO] Shards generados: $N_SHARDS"
if [[ "$N_SHARDS" -eq 1 ]]; then
  echo "[WARN] Solo 1 shard; la paralelización no aportará gran aceleración."
fi

# ========== Función para procesar un shard ==========
run_one_shard() {
  local shard_fna="$1"
  local base
  base="$(basename "$shard_fna" .fna)"         # p.ej., shard_0001
  local faa_out="$SHARD_OUTS/${base}.faa"
  local faa_tagged="$SHARD_OUTS/${base}.tagged.faa"
  local log_err="$LOGS_DIR/${base}.stderr.log"
  local log_out="$LOGS_DIR/${base}.stdout.log"

  # Ejecutar prodigal-gv
  prodigal-gv -i "$shard_fna" -a "$faa_out" -p meta -q 1>"$log_out" 2>"$log_err"

  # Etiquetar cabeceras para unicidad (añade |shard_XXXX)
  awk -v tag="|${base}" 'BEGIN{OFS=""}
    /^>/ { sub(/\r$/, "", $0); print $0, tag; next }
         { print }' "$faa_out" > "$faa_tagged"
}

export -f run_one_shard
export SHARD_OUTS LOGS_DIR

# ========== Backend de paralelización ==========
echo "[INFO] Lanzando prodigal-gv en paralelo con hasta $N_CPUS trabajos..."
if command -v parallel >/dev/null 2>&1; then
  # GNU parallel
  printf "%s\n" "${SHARDS[@]}" | parallel -j "$N_CPUS" --halt now,fail=1 run_one_shard {}
else
  # Fallback a xargs -P (POSIX-ish)
  printf "%s\0" "${SHARDS[@]}" | xargs -0 -n1 -P "$N_CPUS" bash -c 'run_one_shard "$0"'
fi

# ========== Verificación de salidas por shard ==========
mapfile -t TAGGED < <(LC_ALL=C ls -1 "$SHARD_OUTS"/shard_*.tagged.faa 2>/dev/null || true)
if [[ "${#TAGGED[@]}" -ne "$N_SHARDS" ]]; then
  echo "[ERROR] Nº de .tagged.faa (${#TAGGED[@]}) != Nº de shards ($N_SHARDS). Revisa logs en $LOGS_DIR" >&2
  exit 1
fi

# ========== Concatenación determinista ==========
echo "[INFO] Concatenando salidas a $FINAL_FAA ..."
LC_ALL=C ls -1 "$SHARD_OUTS"/shard_*.tagged.faa | sort -V | xargs cat -- > "$FINAL_FAA"

if [[ ! -s "$FINAL_FAA" ]]; then
  echo "[ERROR] No se generó $FINAL_FAA" >&2
  exit 1
fi

echo "[INFO] Prodigal-GV (paralelo) finalizado. Proteínas: $FINAL_FAA"

# ========== Limpieza opcional ==========
if [[ "$KEEP_INTERMEDIATE" == false ]]; then
  rm -rf "$SHARDS_DIR" "$SHARD_OUTS"
  # Conservamos logs por defecto; borra si no los quieres
fi

exit 0
