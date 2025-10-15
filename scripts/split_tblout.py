#!/usr/bin/env python3
from pathlib import Path

def split_tblout(tblout_path, n_parts=2, out_dir=None):
    """
    Divide un archivo .tblout de HMMER3 en n_parts archivos válidos.
    Mantiene la cabecera (#...) en cada subarchivo y corta sólo entre targets distintos.
    """
    tblout_path = Path(tblout_path)
    out_dir = Path(out_dir or tblout_path.parent)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 1) Leer todo el archivo
    with tblout_path.open() as f:
        lines = f.readlines()

    # --- 2) Separar cabecera y datos
    header = [l for l in lines if l.startswith("#")]
    data   = [l for l in lines if not l.startswith("#")]

    if not data:
        print("[WARN] No data lines found.")
        return

    # --- 3) Agrupar por 'target name' (primera columna)
    blocks = []
    current_target = None
    current_block = []
    for line in data:
        target = line.split()[0]
        if target != current_target:
            # nuevo bloque
            if current_block:
                blocks.append(current_block)
            current_block = [line]
            current_target = target
        else:
            current_block.append(line)
    if current_block:
        blocks.append(current_block)

    # --- 4) Calcular tamaño aproximado por archivo
    total_blocks = len(blocks)
    per_file = (total_blocks // n_parts) or 1

    # --- 5) Escribir los archivos resultantes
    for i in range(n_parts):
        start = i * per_file
        end = (i + 1) * per_file if i < n_parts - 1 else total_blocks
        subset = list(chain.from_iterable(blocks[start:end]))

        out_file = out_dir / f"{tblout_path.stem}_part{i+1}.tbl"
        with out_file.open("w") as f:
            f.writelines(header)
            f.writelines(subset)

        print(f"[OK] Written {out_file} ({len(subset)} lines, {len(blocks[start:end])} targets)")

if __name__ == "__main__":
    from itertools import chain
    import sys

    if len(sys.argv) < 3:
        print("Usage: split_tblout.py <file.tblout> <n_parts> [out_dir]")
        sys.exit(1)

    split_tblout(sys.argv[1], int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
