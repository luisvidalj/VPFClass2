#!/usr/bin/env python3
from pathlib import Path
from itertools import islice

def split_tblout_streaming(tblout_path, n_parts=20, out_dir=None):
    tblout_path = Path(tblout_path)
    out_dir = Path(out_dir or tblout_path.parent)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Leer solo la cabecera primero
    with tblout_path.open() as f:
        header = []
        while True:
            pos = f.tell()
            line = f.readline()
            if not line.startswith("#"):
                f.seek(pos)  # retrocede a la primera línea de datos
                break
            header.append(line)

        # 2) Primer pase: contar cuántos targets distintos hay
        total_targets = 0
        last_target = None
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            target = line.split()[0]
            if target != last_target:
                total_targets += 1
                last_target = target

    print(f"[INFO] Found {total_targets:,} targets in total")

    # 3) Calcular cuántos targets por split
    per_file = (total_targets // n_parts) or 1

    # 4) Segundo pase: escribir los splits al vuelo
    with tblout_path.open() as f:
        # saltar la cabecera
        while True:
            pos = f.tell()
            line = f.readline()
            if not line.startswith("#"):
                f.seek(pos)
                break

        part_idx = 1
        targets_in_current = 0
        current_target = None
        out_file = out_dir / f"{tblout_path.stem}_part{part_idx}.tbl"
        out_handle = out_file.open("w")
        out_handle.writelines(header)

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            target = line.split()[0]
            # Nuevo target -> contar
            if target != current_target:
                targets_in_current += 1
                current_target = target
                # ¿Cambiamos de archivo?
                if targets_in_current > per_file and part_idx < n_parts:
                    out_handle.close()
                    print(f"[OK] Written {out_file}")
                    part_idx += 1
                    out_file = out_dir / f"{tblout_path.stem}_part{part_idx}.tbl"
                    out_handle = out_file.open("w")
                    out_handle.writelines(header)
                    targets_in_current = 1  # ya contamos este target

            out_handle.write(line)

        out_handle.close()
        print(f"[OK] Written {out_file}")
        print("[DONE] Splitting complete!")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: split_tblout_streaming.py <file.tblout> <n_parts> [out_dir]")
        sys.exit(1)

    split_tblout_streaming(sys.argv[1], int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
