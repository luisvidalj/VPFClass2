#!/bin/bash

# Uso: ./run_hmm_parallel.sh PROT_DIR HMMS_DIR OUT_DIR NUM_CPUS

# 1. Parámetros
PROT_DIR=$1       # archivo .faa con proteínas predichas
HMMS_DIR=$2       # archivo .hmms con todos los modelos juntos
OUT_DIR=$3        # carpeta para guardar splits y outputs
NUM_CPUS=$4       # número total de CPUs disponibles (ej. 20)

# 2. Verificaciones básicas
if [[ -z "$PROT_DIR" || -z "$HMMS_DIR" || -z "$OUT_DIR" || -z "$NUM_CPUS" ]]; then
    echo "[ERROR] Uso: $0 PROT_DIR HMMS_DIR OUT_DIR NUM_CPUS"
    exit 1
fi

if ! command -v hmmsearch &> /dev/null; then
    echo "[ERROR] No se encontró 'hmmsearch' en el PATH."
    exit 1
fi

# 3. Calcular número de splits
SPLITS=$((NUM_CPUS / 2))
if [[ $SPLITS -lt 1 ]]; then
    SPLITS=1
fi

# 4. Contar secuencias
NUM_SEQS=$(grep -c ">" "$PROT_DIR")
SEQ_PER_SPLIT=$(echo "($NUM_SEQS + $SPLITS - 1) / $SPLITS" | bc)

echo "[INFO] Número total de secuencias: $NUM_SEQS"
echo "[INFO] Número de splits: $SPLITS (usando 2 CPUs por split)"
echo "[INFO] Secuencias por archivo: $SEQ_PER_SPLIT"
echo

# 5. Crear carpeta de salida
mkdir -p "$OUT_DIR"

# 6. Dividir el archivo .faa
awk -v n="$SEQ_PER_SPLIT" -v dir="$OUT_DIR" '/^>/ {
    if (count++ % n == 0) {
        file = dir "/split_" ++i ".faa"
    }
} {
    print > file
}' "$PROT_DIR"

# 7. Ejecutar hmmsearch en paralelo
echo "[INFO] Ejecutando hmmsearch sobre cada split..."

SECONDS=0

for file in "$OUT_DIR"/split_*.faa; do
    BASENAME=$(basename "${file%.faa}")
    hmmsearch --cpu 2 --tblout "$OUT_DIR/${BASENAME}.tbl" "$HMMS_DIR" "$file" &
done

wait

# 8. Reporte de duración
DURATION=$SECONDS
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

printf "[INFO] Tiempo total: %02d:%02d:%02d (hh:mm:ss)\n" $HOURS $MINUTES $SECONDS

