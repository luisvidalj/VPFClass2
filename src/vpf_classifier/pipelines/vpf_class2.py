from __future__ import annotations

from pathlib import Path
import time
import shutil
import json
import torch
import torch.nn.functional as F
import numpy as np
import sys
import argparse

from scipy.sparse import load_npz, csr_matrix, issparse
from dataclasses import dataclass
from scipy.sparse import save_npz, csr_matrix, issparse

from vpf_classifier.parsers.fasta_parser import FastaParser
from vpf_classifier.parsers.prodigal_parser import Prodigal
from vpf_classifier.parsers.vpf_parser import VPF_parser

import pandas as pd



##### --------------------------- Helpers --------------------------------------------------

def _save_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _load_idx_to_label(path: Path) -> list[str]:
    data = json.loads((path / "idx_to_label.json").read_text(encoding="utf-8"))
    if isinstance(data,list):
        return data
    # dict {idx:lab} -> llista de tuples ordenades segons el primer element de la tupla (l'idx)
    return [lab for _,lab in sorted(data.items(), key=lambda kv: int(kv[0]))]



def _scipy_csr_to_torch_coo_batch(X_csr_batch: csr_matrix, device:str = "cpu") -> torch.Tensor:
    """
    Convert (BxD) Scipy CSR to a torch.sparse_coo_tensor coalesced (as this is what SparseNN is expecting)
    (Lo que esperam de sa collate_fn que hem dissenyat a sa xarxa)
    """

    assert issparse(X_csr_batch)
    coo = X_csr_batch.tocoo(copy=False)
    indices = np.vstack((coo.row, coo.col)).astype(np.int64, copy=False)
    values = coo.data.astype(np.float32, copy=False)
    shape = (X_csr_batch.shape[0], X_csr_batch.shape[1])

    i = torch.from_numpy(indices)
    v = torch.from_numpy(values)
    sp = torch.sparse_coo_tensor(i,v, size=shape, device="cpu").coalesce()
    if device != "cpu":
        sp = sp.to(device=device)

    return sp


def _load_lineage_by_genus(model_dir_p: Path) -> dict:
    """
    Load tool_data/MSL_labelling/MSLXX/lineage.json
    """

    # model_dir_p = tool_data/<base>/models/MSLXX
    try:
        msl_tag = model_dir_p.name 
        tool_data_root = model_dir_p.parents[3] #.../tool_data
    except Exception:
        return {}
    
    lineage_path = tool_data_root / "MSL_labelling" / msl_tag / "lineage.json"
    print(f"Just in case {tool_data_root}")
    if lineage_path.exists():
        return json.loads(lineage_path.read_text(encoding="utf-8"))
    return {}



# def _predict_in_batches(
#         model: torch.nn.Module,
#         X_csr: csr_matrix,
#         idx_to_label: list[str],
#         device: str = "cpu",
#         batch_size: int = 1024,
#         topk: int = 3,
# ) -> list[dict]:
#     """
#     Function to infer predictions from the pre-loaded model.
#     """
#     # per com esta construit el codi no tocaria mai trobar aquest error
#     assert issparse(X_csr), "Error when calling _predict_in_batches due to not sparse X_Csr" # aixo amollara un error si no es compleix la condició
#     n = X_csr.shape[0]
#     records: list[dict] = []
#     k = min(topk, len(idx_to_label))

#     model.eval()
#     model.to(device)

#     with torch.no_grad():
#         for start in range(0,n, batch_size):
#             end = min(start + batch_size, n)



# ## Meter a partir de predict_ in batches




##### --------------------------- Directories parser --------------------------------------------------

@dataclass
class RunDirs:
    root: Path
    logs: Path
    fasta: Path
    prodigal: Path
    hmmer: Path
    features: Path
    preds: Path

def _make_run_dirs(outdir_p: Path) -> RunDirs:
    rd = RunDirs(
        root = outdir_p,
        logs = outdir_p / "logs",
        fasta = outdir_p / "fasta",
        prodigal = outdir_p / "1_prodigal",
        hmmer = outdir_p / "2_hmmer",
        features = outdir_p / "features",
        preds=outdir_p / "0_preds"
    )
    for d in (rd.prodigal, rd.hmmer, rd.features, rd.preds):
        d.mkdir(parents=True, exist_ok=True)

    return rd


import sys
from pathlib import Path
from typing import Optional

def _ensure_sparse_genus_on_path(model_dir_p: Path) -> None:
    """
    Inserta '.../models' en sys.path si existe '.../models/sparse_genus',
    para que 'import sparse_genus' funcione al hacer torch.load(model).
    """
    try:
        # Caso típico: model_dir_p = repo/models/complete_markers/MSL40
        models_root = model_dir_p.parents[3]
        models_root = Path(models_root) / "models"
        if (models_root / "sparse_genus").exists():
            path_str = str(models_root.resolve())
            if path_str not in sys.path:
                sys.path.insert(0, path_str)
                print(f"[DEBUG] Added to sys.path: {path_str}")
    except Exception:
        # Fallback: buscar subiendo desde CWD
        cwd = Path.cwd()
        for r in [cwd] + list(cwd.parents):
            cand = r / "models" / "sparse_genus"
            if cand.exists():
                path_str = str((r / "models").resolve())
                if path_str not in sys.path:
                    sys.path.insert(0, path_str)
                    # print(f"[DEBUG] Added to sys.path: {path_str}")
                break


from vpf_classifier.utils.config import ROOT_DIR
import sys

_MODELS_CODE_ADDED = False

def _ensure_models_code_on_path_once() -> None:
    """
    Insert '<repo>/models' into sys.path to see sparse_genus and sparse_family modules
    """

    global _MODELS_CODE_ADDED

    if _MODELS_CODE_ADDED:
        return
    
    models_dir = (ROOT_DIR / "models").resolve()
    path_str = str(models_dir)

    if models_dir.exists() and path_str not in sys.path:
        sys.path.insert(0,path_str) # insert es como append pero donde quieras. Igual funcionaria con append
    _MODELS_CODE_ADDED = True


def _ensure_sparse_task_on_path(model_dir_p: Path, module_dir_names: str) -> None:


    def _add_models_path(models_dir: Path) -> None:
        path_str = str(models_dir.resolve())
        if path_str not in sys.path:
            sys.path.insert(0, path_str)

    def _find_models_root_upwards(start:Path) -> Optional[Path]:
        
        pass



###### ------------------------ Main function ----------------------------------------------

def run_user_pipeline(
        fasta: str | None = None,
        outdir: str | None = None,
        model_dir: str | None = None,
        vpf_dict: str | None = None,
        hmm_models: str | None = None,
        e_value_threshold: float = 1e-3,
        num_cpus: int = 8,
        device: str = "cuda",  # ojo luego si no hay GPU
        topk: int = 3,
        #batch_size: int = 100 (no se si afegir-ho)
):
    # ========= SETUP ========================================

    # 1) Input Validation
    if fasta is None:
        raise ValueError("You must provide --fasta")
    if model_dir is None:
        raise ValueError("Models directory not found")
    if vpf_dict is None:
        raise ValueError("Tool_data information not found")

    fasta_p = Path(fasta).expanduser().resolve()
    model_dir_p = Path(model_dir).expanduser().resolve()
    vpf_dict_p = Path(vpf_dict).expanduser().resolve()

    if outdir is None:
        outdir_p = Path("runs") / f"inference_{time.strftime('%Y%m%d-%H%M%S')}"
    else:
        outdir_p = Path(outdir).expanduser().resolve()
    outdir_p.mkdir(parents=True, exist_ok=True)

    # subcarpetas por run
    run_dirs = _make_run_dirs(outdir_p)
    print(f"MODEL DIR: {model_dir_p}")

    # 2) Existence checks
    if not fasta_p.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_p}")
    if not model_dir_p.exists():
        raise FileNotFoundError(f"Model bundle dir not found: {model_dir_p}")
    if not (model_dir_p).exists():
        raise FileNotFoundError(f"Missing model directories in {model_dir_p}")
    if not vpf_dict_p.exists():
        raise FileNotFoundError(f"vpf_to_index JSON not found: {vpf_dict_p}")

    # 3) Required Software
    prodigal_path = shutil.which("prodigal-gv")   
    hmmsearch_path = shutil.which("hmmsearch")

    print("=========================== RUNTIME SETUP ===========================")
    print(f"[SETUP] FASTA:        {fasta_p}")
    print(f"[SETUP] OUTDIR:       {outdir_p}")
    print(f"[SETUP] MODEL DIR:    {model_dir_p}")
    print(f"[SETUP] VPF DICT:     {vpf_dict_p}")
    print(f"[SETUP] e-value thr.: {e_value_threshold}")
    print(f"[SETUP] num_cpus:     {num_cpus}")
    print(f"[SETUP] device:       {device}")
    print(f"[SETUP] prodigal-gv:     {prodigal_path or 'NOT FOUND'}")
    print(f"[SETUP] hmmsearch:    {hmmsearch_path or 'NOT FOUND'}")
    print(f"[SETUP] Generating subdirs:")
    #print(f"    - fasta:     {run_dirs.fasta}")
    print(f"    - prodigal:  {run_dirs.prodigal}")
    print(f"    - hmmer:     {run_dirs.hmmer}")
    #print(f"    - features:  {run_dirs.features}")
    print(f"    - preds:     {run_dirs.preds}")
    print("=====================================================================")

    if prodigal_path is None:
        raise EnvironmentError("Prodigal not found. Make sure the conda environment is activated.")
    if hmmsearch_path is None:
        raise EnvironmentError("hmmsearch not found. Make sure the conda environment is activated.")

    # 4) Run metadata (setup)
    meta = {
        "fasta": str(fasta_p),
        "outdir": str(outdir_p),
        "model_dir": str(model_dir_p),
        "vpf_dict": str(vpf_dict_p),
        "e_value_threshold": float(e_value_threshold),
        "num_cpus": int(num_cpus),
        "device": device,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "stage": "setup_done",
        "run_dirs": {
            #"fasta": str(run_dirs.fasta),
            "prodigal": str(run_dirs.prodigal),
            "hmmer": str(run_dirs.hmmer),
            #"features": str(run_dirs.features),
            "preds": str(run_dirs.preds),
            # "logs": str(run_dirs.logs),
        },
    }
    (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print("--------------------------------------------------------------------")
    # ========= Fasta -> DataFrame ===========================================
    print("[PIPELINE] 1/5 - Parsing FASTA...")
    fasta_parser = FastaParser(fna_path=str(fasta_p))
    fasta_parser.parse_fasta_to_dataframe()

    full_accessions = fasta_parser.ncbi_df['Accession'].tolist()

    # Guardar cabeceras en subcarpeta FASTA
    # headers_csv = run_dirs.fasta / "fasta_headers.csv"
    # _save_csv(fasta_parser.ncbi_df, headers_csv)
    # print(f"[PIPELINE] FASTA parsed. Headers saved in: {headers_csv}")

    # Update metadata
    meta.update({
        "stage": "fasta_parsed",
        "n_sequences": int(len(getattr(fasta_parser, "ncbi_df", []))),
        #"fasta_headers_csv": str(headers_csv),
    })
    (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    # ======== Prodigal =========================================================
    print("--------------------------------------------------------------------")
    # (si en algún momento ofrecemos opción de cargar .faa directamente, aquí haríamos el bypass)
    print("[PIPELINE] 2/5 - Prodigal...")
    #fasta_parser.run_prodigal(output_dir=str(run_dirs.prodigal))
    prodigal = Prodigal(parser=fasta_parser, output_dir=str(run_dirs.prodigal))  # <- ahora en subcarpeta
    df_proteins = prodigal.parse_prodigal()

    # Guardar resumen de proteínas en subcarpeta PROTEINS
    # proteins_csv = run_dirs.prodigal / "prodigal_proteins.csv"
    # _save_csv(df_proteins, proteins_csv)
    # print(f"[PIPELINE] Prodigal completed. Summary of viral proteins: {proteins_csv}")

    # Registrar salidas presentes (si no existen, listas vacías)
    prodigal_files = {
        "faa": [str(p.resolve()) for p in run_dirs.prodigal.glob("*.faa")],
        "gff": [str(p.resolve()) for p in run_dirs.prodigal.glob("*.gff")],
        "fna": [str(p.resolve()) for p in run_dirs.prodigal.glob("*.fna")],
    }

    meta.update({
        "stage": "prodigal_done",
        "n_proteins": int(len(df_proteins)),
        "prodigal_outputs": prodigal_files,
        #"prodigal_proteins_csv": str(proteins_csv),
    })
    (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


    # ========= HMMER (Bloque 4) ==========================================
    print("--------------------------------------------------------------------")
    print(f"[PIPELINE] 3/5 - HMMER + filtering (e_value ≤ {e_value_threshold})…")

    if hmm_models is None:
        raise ValueError("You must either provide hmm_models or ensure that it is in the default path.")
    hmm_models_p = Path(hmm_models).expanduser().resolve()
    if not hmm_models_p.exists():
        raise FileNotFoundError(f"HMM profiles not found: {hmm_models_p}")
    
    """
    Aqui chatgpt sugiere que ejecute

    fasta_parser.run_hmmer(
        hmm_models=hmm_models_p,
        output_dir=run_dirs.hmmer,
        num_cpus=int(num_cpus)
    )

    Pero en teoria esto deberia parser con simplemente llamar a vpf.parse_multiple_hmm()
    """

    # VPF_parser debería reutilizar .tbl si ya existen en run_dirs.hmmer
    vpf = VPF_parser(
        parser=fasta_parser,
        hmm_file=str(hmm_models_p),                    # <-- importante
        e_value_threshold=float(e_value_threshold),
        num_cpus=int(num_cpus),
        user = True,
        vpf_dict_path=Path(vpf_dict_p),                # <-- importante (evita el "hack" de _vpf_dict_path)
        hmm_output_dir=run_dirs.hmmer,                 # <-- importante: nuestros .tbl van a outdir/hmmer/
        vector_norm="l2"
    )

    vpf.parse_multiple_hmm_parallel()
    

    # Guardamos los hits por virus 
    # hits_csv = run_dirs.hmmer / "hmm_hits.csv"
    # vpf.df_virus_hmm.to_csv(hits_csv, index=False)
    # print(f"[PIPELINE] HMMER completed. Hits per virus saved in: {hits_csv}")

    # Resumen seguro (si faltan columnas, no reventamos)
    try:
        n_pairs = int(len(vpf.df_virus_hmm))
    except Exception:
        n_pairs = None
    try:
        n_viruses_with_hits = int(vpf.df_virus_hmm["Accession"].nunique())
    except Exception:
        n_viruses_with_hits = None

    meta.update({
        "stage": "hmmsearch_done",
        "e_value_threshold": float(e_value_threshold),
        "num_cpus": int(num_cpus),
        "n_hit_pairs": n_pairs,
        "n_viruses_with_hits": n_viruses_with_hits,
        #"hmm_hits_csv": str(hits_csv),
        "hmm_models": str(hmm_models_p),
    })
    (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


    
# ========= 5: Features Building (VPF counts aggregated by Virus) =======
    print("--------------------------------------------------------------------")
    print("[PIPELINE] 4/5 - Building features (counts per VPF)…")

    # 1) Matriz dispersa CSR fija construida por VPF_parser
    sparse_mat = vpf.get_sparse_matrix_from_dataframe()  # csr_matrix
    if not issparse(sparse_mat):
        sparse_mat = csr_matrix(sparse_mat)

    accessions = vpf.df_virus_hmm["Accession"].tolist()
    print(accessions[1:5])
    print(f"[INFO] Detected Accessions at .tbl files: {len(accessions)}")
    missing_acccessions = [acc for acc in full_accessions if acc not in accessions]
    vpf_to_index = getattr(vpf, "vpf_to_index", None)
    if vpf_to_index is None:
        print("Here 1")
        # fallback: leer del json pasado por el usuario
        with open(vpf_dict_p, "r", encoding="utf-8") as f:
            vpf_to_index = json.load(f)

    # 2) Guardados
    features_npz_path = run_dirs.features / "features_counts_sparse.npz"
    save_npz(features_npz_path, sparse_mat)
    print(f"[INFO] Features saved")

    accessions_path = run_dirs.features / "accessions.txt"
    with open(accessions_path, "w", encoding="utf-8") as f:
        for acc in accessions:
            f.write(f"{acc}\n")

    # Copiamos/guardamos el mapeo VPF->índice que se usó realmente
    vpf_map_path = run_dirs.features / "vpf_to_index.used.json"
    with open(vpf_map_path, "w", encoding="utf-8") as f:
        json.dump(vpf_to_index, f, ensure_ascii=False, indent=2)

    # 3) Estadísticas
    #n_rows, n_cols = sparse_mat.shape
    #nnz = int(sparse_mat.nnz)
    #nonzero_fraction = float(nnz / (n_rows * n_cols)) if (n_rows and n_cols) else 0.0
    # features_stats = {
    #     "shape": [int(n_rows), int(n_cols)],
    #     "nnz": nnz,
    #     "nonzero_fraction": nonzero_fraction,
    #     "n_sequences_for_features": int(len(accessions)),
    #     "dict_size": int(len(vpf_to_index)),
    #     "note": "counts per VPF (fixed length), built by VPF_parser._add_vpf_counts_sparse_fixed",
    # }
    #(run_dirs.features / "features_stats.json").write_text(json.dumps(features_stats, indent=2), encoding="utf-8")

    # print(f"[INFO] Features (sparse CSR) saved: {features_npz_path}")
    # print(f"[PIPELINE] Virus Accessions: {accessions_path}")
    # print(f"[PIPELINE] VPF map: {vpf_map_path}")
    # print(f"[PIPELINE] Stats: {features_stats}")

    # 4) Meta
    meta.update({
        "stage": "features_built",
        "features_counts_sparse_npz": str(features_npz_path),
        "accessions_txt": str(accessions_path),
        "vpf_to_index_used_json": str(vpf_map_path),
        #"features_stats": features_stats,
    })
    (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


    # =============== 6: Model + Inference ========================================================
    print("--------------------------------------------------------------------")
    print("[PIPELINE] 5/5 - Loading model(s) and generating predictions...")

    # 1) features (CSR), accessions and labels
    features_npz_path = run_dirs.features /"features_counts_sparse.npz" 
    accessions_path = run_dirs.features / "accessions.txt"

    # idx_to_label = _load_idx_to_label(path=model_dir_p)

    if not features_npz_path.exists():
        raise FileNotFoundError(f"Expected features at {features_npz_path}")
    X_csr = load_npz(features_npz_path)

    if not issparse(X_csr):
        X_csr = csr_matrix(X_csr)

    if not accessions_path.exists():
        raise FileNotFoundError(f"Expected accessions list at {accessions_path}")
    # la idea del codi d'abaix es recollir en una llista totes els accessions (un per linia al .txt)
    accessions2 = [line.strip() for line in accessions_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    print(f"Accessions abans de guardar: {len(accessions)}")
    print(f"Accessions abans de guardar: {len(accessions2)}")

    # model_path = model_dir_p / "model.pt"
    in_features = int(X_csr.shape[1])
    # num_genus   = int(len(idx_to_label))

    # 1b) Resolve task-specific model directories (genus/family) with backward-compat
    # Should we care about majus/minus????????
    genus_dir = (model_dir_p / "Genus")
    family_dir = (model_dir_p / "Family")

    available_tasks: list[tuple[str, Path]] = []
    # Prefer new layout (subfolders). If not present, fall back to former (model.pt directly)
    if (genus_dir / "model.pt").exists():
        available_tasks.append(("genus",genus_dir))
    elif (model_dir_p / "model.pt").exists():
        available_tasks.append(("genus",model_dir_p)) # former: genus-only layout

    if (family_dir / "model.pt").exists():
        available_tasks.append(("family", family_dir))

    
    if not available_tasks:
        raise FileNotFoundError(
        f"No models found under {model_dir_p}. "
        f"Expected either {model_dir_p/'model.pt'} or {genus_dir/'model.pt'} / {family_dir/'model.pt'}."
    )

    print(f"[INFO] Detected tasks: {', '.join(t for t,_ in available_tasks)}")

        
    # 2) Device check
    if device == "cuda" and not torch.cuda.is_available():
        print("[WARN] CUDA not available. Falling back to CPU.")
        device = "cpu" 
    
    batch_size = 100

    # 3) Inference poer task (same features)
    task_dfs: dict[str, pd.DataFrame] = {}

    # 3.0) Assure import of sparse modules
    _ensure_models_code_on_path_once()

    for task, task_dir in available_tasks:
        # 3.1) Load labels for this task
        idx_to_label = _load_idx_to_label(path=task_dir)
        num_labels = int(len(idx_to_label))
        k = min(topk, num_labels) # ESTO PUEDE SER PROBLEMATICO AUNQUE NUNCA PASARA

        model_path = task_dir / "model.pt"

        try:
            # Cargamos SIEMPRE el objeto completo; evita el FutureWarning al fijar weights_only=False
            model = torch.load(model_path, map_location=device, weights_only=False)
        except TypeError:
            # Para compatibilidad con PyTorch<2.5 que no tiene weights_only
            model = torch.load(model_path, map_location=device)

        if not isinstance(model, torch.nn.Module):
            raise RuntimeError(
            "This pipeline expects a full model object saved via torch.save(model, 'model.pt'). "
            "It looks like the file contains a state_dict; please re-export saving the full object."
            )

        model.eval()

        # 3.4) Friendly dimensional check
        fc1 = getattr(model, "fc1", None)
        if fc1 is not None and hasattr(fc1, "in_features"):
            if int(fc1.in_features) != in_features:
                print(f"[WARN] Model in_features={int(fc1.in_features)} but features have {in_features}. "
                  "If it does not match, check that vpf_to_index and the model correspond to the same version..")


        # 3.5) Batches --> Predictions
        print(f"[INFO] Running task={task} (C={num_labels}) …")
        records = []
        with torch.no_grad():
            for start in range(0, X_csr.shape[0], batch_size):
                end = min(start + batch_size, X_csr.shape[0])
                X_batch_csr = X_csr[start:end]
                x_sp = _scipy_csr_to_torch_coo_batch(X_csr_batch=X_batch_csr, device=device)

                logits = model(x_sp)
                probs = F.softmax(logits, dim=1)

                top1_scores, top1_idx = probs.max(dim=1)
                tk_scores, tk_idx = probs.topk(k=k, dim=1) # torch.topk()

                if task == "genus":
                    pred_col = "pred_genus"
                    score_col = "score_genus"
                    topk_col = "topk_genus"
                
                else: # family
                    pred_col = "pred_family"
                    score_col = "score_family"
                    topk_col = "topk_family"

                for i in range(end - start):
                    pred = idx_to_label[int(top1_idx[i].item())] 
                    score = float(top1_scores[i].item())
                    tkl = [idx_to_label[int(j.item())] for j in tk_idx[i]]
                    tks = [f"{float(s.item()):.2f}" for s in tk_scores[i]]
                    records.append({
                        pred_col: pred,
                        score_col: score,
                        topk_col: "; ".join(f"{g} ({s})" for g, s in zip(tkl, tks))
                    })

                del x_sp, logits, probs, top1_scores, top1_idx, tk_scores, tk_idx

            task_dfs[task] = pd.DataFrame(records)

    # 4) Store preds in a CSV

    # 4.0) Accessions length check
    n_rows = int(X_csr.shape[0])
    if len(accessions) != n_rows:
        raise ValueError(
            "[ERROR] Row mismatch: features rows != accessions .\n"
            f" - X_csr rows: {n_rows}\n"
            f" - accessions: {len(accessions)}\n"
            "Asegurate de que el orden es el que toca"
        )
    
    preds_df = pd.DataFrame({"Accession": accessions})


    for _, df_task in task_dfs.items():
        preds_df = pd.concat([preds_df, df_task], axis=1)


    # Orden de rangos con Mayúsculas (como en tu lineage.json)
    RANKS = ["Realm","Subrealm","Kingdom","Subkingdom","Phylum","Subphylum",
            "Class","Subclass","Order","Suborder","Family","Subfamily","Genus","Subgenus","Species"]


    def apply_lineage_and_tidy(
        preds_df: pd.DataFrame,
        model_dir_p: Path,  # .../tool_data/<base>/models/MSLxx
        lineage_filename: str = "lineage.json",  # en tool_data/MSL_labelling/MSLxx/lineage.json
    ) -> pd.DataFrame:
        """
        - Elige fuente (Genus vs Family) comparando score_genus vs score_family.
        - Reconstruye Lineage (hasta Genus o hasta Family según la fuente elegida).
        - Redondea scores a 3 decimales y ordena columnas.
        """

        df = preds_df.copy()

        # Asegura columnas esperadas (si faltan, vacías para no romper)
        for c in ["pred_genus","score_genus","topk_genus","pred_family","score_family","topk_family"]:
            if c not in df.columns:
                df[c] = np.nan

        # 2) Carga lineage por GENUS (con claves de rangos Mayúsculas)
        msl_tag = model_dir_p.name
        try:
            msl_tag = model_dir_p.name
            tool_data_root = model_dir_p.parents[2]
            lineage_path = tool_data_root / "MSL_labelling" / msl_tag / lineage_filename
            # print(lineage_path)
            lineage_by_genus = json.loads(lineage_path.read_text(encoding="utf-8"))
        except Exception:
            lineage_by_genus = {}

        # Preconstruye un mapa Family -> “nodo” (tomado del primer género que pertenezca a esa familia)
        family_node = {}
        for g, node in lineage_by_genus.items():
            fam = node.get("Family")
            if fam and fam not in family_node:
                # Copia superficial con rangos Mayúsculas
                family_node[fam] = {k: v for k, v in node.items() if k in RANKS}

        def _join_lineage(node: dict | None, up_to: str) -> str | None:
            """
            Build lineage from Realm to 'up__to' (keeping missing ranks)
            """
            if not node:
                return None
            out = []
            for r in RANKS:
                val = node.get(r,"")
                out.append(val)
                if r == up_to:
                    break
            return "; ".join(out)
            
        # 3) Fila a fila: escoger fuente y reconstruir
        def _build_lineage(row):
            g = row.get("pred_genus")
            sg = row.get("score_genus")
            f = row.get("pred_family")
            sf = row.get("score_family")

            # normaliza
            g = None if g in (None, "", "Unknown") else g
            f = None if f in (None, "", "Unknown") else f
            sg = float(sg) if pd.notna(sg) else float("-inf")
            sf = float(sf) if pd.notna(sf) else float("-inf")


            # --- 1: If genus_score > family_score ---
            if sg >= sf and g in lineage_by_genus:
                node = lineage_by_genus.get(g)
                return _join_lineage(node, up_to="Genus")
            # --- 2: If genus_score > family_score ---
            elif f in family_node:
                node = family_node.get(f)
                lineage_to_family = _join_lineage(node, up_to="Family")

                add_genus = False
                if g:
                    if g not in lineage_by_genus:
                        add_genus = True
                    else:
                        fam_of_g = lineage_by_genus[g].get("Family")
                        add_genus = (fam_of_g == f)

                if add_genus:
                    return lineage_to_family + "; " + g
                else:
                    return lineage_to_family
            
            else:
                # últimos recursos: si no hay familia válida pero hay genus en el mapa
                if g in lineage_by_genus:
                    return _join_lineage(lineage_by_genus.get(g), up_to="Genus")
            return None

        df["Lineage"] = df.apply(_build_lineage, axis=1)

        # 4) Redondeo de scores (3 decimales) si existen
        for sc in ("score_genus","score_family"):
            if sc in df.columns:
                df[sc] = pd.to_numeric(df[sc], errors="coerce").round(3)

        # 5) Reordenar columnas para una salida limpia
        order = [c for c in ["Accession",
                            "pred_genus","score_genus","topk_genus",
                            "pred_family","score_family","topk_family",
                            "Lineage"]
                if c in df.columns]
        # Mantén también cualquier otra columna que hubiera (por si añadiste extras)
        rest = [c for c in df.columns if c not in order]
        df = df[order + rest]

        # 6) Cambiar la prediccion de familia si el score lo requiere
        genus_to_family = {genus: data.get("Family") for genus, data in lineage_by_genus.items()}
        mask = df["score_genus"] > df["score_family"]
        df.loc[mask, "pred_family"] = df.loc[mask, "pred_genus"].map(genus_to_family)
        df.loc[mask, "score_family"] = df.loc[mask, "score_genus"]

        return df


    preds_df = apply_lineage_and_tidy(preds_df=preds_df, model_dir_p=model_dir_p)
    preds_csv = run_dirs.preds / "prediction.csv"
    print(preds_df.head())

    if missing_acccessions:
        missing_df = pd.DataFrame({
            "Accession": missing_acccessions
        })
        print(f"[INFO] {len(missing_acccessions)} virus can not be annotated")

        for col in preds_df.columns:
            if col not in missing_df.columns:
                missing_df[col] = np.nan

        preds_df = pd.concat([preds_df, missing_df], ignore_index=True)
        _save_csv(df=preds_df, path=preds_csv)

    print(f"[PIPELINE] Predictions stored in: {preds_csv}")


    # 5) Meta final
    meta.update({
        "stage": "done",
        "n_predicted": int(len(preds_df)),
        "predictions_csv": str(preds_csv),
        "device_used": device,
    })
    (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    if run_dirs.features.exists():
        shutil.rmtree(run_dirs.features)


    
    return fasta_p, outdir_p, model_dir_p, vpf_dict_p, fasta_parser





    





# FALTA AÑADIR LA CLI Y ACLARAR SI LAS RUTAS FUNCIONAN. UNA VEZ FUNCIONEN, HAY QUE VER COMO SUBIR LAS COSAS Y SIMULAR UN EJERCICIO
# DE USUARIO REAL










































# def run_user_pipeline(
#         fasta: str | None = None,
#         outdir: str | None = None,
#         model_dir: str | None = None,
#         vpf_dict: str | None = None,
#         e_value_threshold: float = 1e-3,
#         num_cpus: int = 8,
#         device: str = "cuda", # A esto hay que echarle un ojo
#         save_features: bool = True,
#         topk: int = 3,
# ):
    
#     # ========= SETUP ========================================

#     # 1) Input Validation
    
#     if fasta is None:
#         raise ValueError("You must provide --fasta")
#     if model_dir is None:
#         raise ValueError("You must provide --model-dir")
#     if vpf_dict is None:
#         raise ValueError("You must provide --vpf-dict")
    
#     fasta_p = Path(fasta).expanduser().resolve()
#     model_dir_p = Path(model_dir).expanduser().resolve()
#     vpf_dict_p = Path(vpf_dict).expanduser().resolve()

#     if outdir is None:
#         outdir_p = Path("runs") / f"inference_{time.strftime('%Y%m%d-%H%M%S')}"
#     else:
#         outdir_p = Path(outdir_p).expanduser().resolve()
#     outdir_p.mkdir(parents=True, exist_ok=True)

#     run_dirs = _make_run_dirs(outdir_p)

#     # 2) Existence checks:
#     if not fasta_p.exists():
#         raise FileNotFoundError(f"FASTA not found: {fasta_p}")
#     if not model_dir_p.exists():
#         raise FileNotFoundError(f"Model bundle dir not found: {model_dir_p}")
#     if not (model_dir_p / "model.pt").exists():
#         raise FileNotFoundError(f"Missing model.pt in {model_dir_p}")
#     if not (model_dir_p / "idx_to_label.json").exists():
#         raise FileNotFoundError(f"Missing idx_to_label.json in {model_dir_p}")
#     if not vpf_dict_p.exists():
#         raise FileNotFoundError(f"vpf_to_index JSON not found: {vpf_dict_p}")
    
#     # 3) Required Software
#     prodigal_path = shutil.which("prodigal-g")
#     hmmsearch_path = shutil.which("hmmsearch")

#     print("[PIPELINE] ========= RUNTIME SETUP =========")
#     print(f"[PIPELINE] FASTA:        {fasta_p}")
#     print(f"[PIPELINE] OUTDIR:       {outdir_p}")
#     print(f"[PIPELINE] MODEL DIR:    {model_dir_p}")
#     print(f"[PIPELINE] VPF DICT:     {vpf_dict_p}")
#     print(f"[PIPELINE] e-value thr.: {e_value_threshold}")
#     print(f"[PIPELINE] num_cpus:     {num_cpus}")
#     print(f"[PIPELINE] device:       {device}")
#     print(f"[PIPELINE] prodigal:     {prodigal_path or 'NOT FOUND'}")
#     print(f"[PIPELINE] hmmsearch:    {hmmsearch_path or 'NOT FOUND'}")
#     print("[PIPELINE] =================================")

#     if prodigal_path is None:
#         raise EnvironmentError("Prodigal not found. Make sure the conda environment is activated")
#     if hmmsearch_path is None:
#         raise EnvironmentError("hmmsearch not found. Make sure the conda environment is activated")
    
#     # 4) Run metadata
#     meta = {
#         "fasta": str(fasta_p),
#         "outdir": str(outdir_p),
#         "model_dir": str(model_dir_p),
#         "vpf_dict": str(vpf_dict_p),
#         "e_value_threshold": float(e_value_threshold),
#         "num_cpus": int(num_cpus),
#         "device": device,
#         "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
#         "stage": "setup_done"
#     }

#     (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

#     # ========= Fasta -> DataFrame ===========================================
#     print("[PIPELINE] 1/5 - Parsing FASTA...")
#     fasta_parser = FastaParser(fna_path=str(fasta_p))
#     fasta_parser.parse_fasta_to_dataframe()

#     # Save header's map
#     headers_csv = outdir / "fasta_headers.csv"
#     _save_csv(fasta_parser.ncbi_df, headers_csv)
#     print(f"[PIPELINE] FASTA parsed. Headers saved in: {headers_csv}")

#     # Update de the metadata
#     meta.update({
#         "stage": "fasta_parsed",
#         "n_sequences": int(len(getattr(fasta_parser, "ncbi_df", [])))
#     })
#     (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


#     # ======== Prodigal =========================================================
#     #### Hay que pensar una manera en la que el usuario pueda elegir introducir proteinas
#     print("[PIPELINE] 2/5 - Prodigal (reusing if .faa is found)...")
#     prodigal = Prodigal(parser=fasta_parser, output_dir=str(outdir_p))
#     df_proteins = prodigal.parse_prodigal()

#     # Guardam el df de proteines a un csv? Per si el usuari ho vol consultar?
#     proteins_csv = outdir_p
#     _save_csv(df_proteins,proteins_csv)
#     print(f"[PIPELINE] Prodigal completed. Summary of viral proteins: {proteins_csv}")

#     try:
#         out_faa = Path(getattr(prodigal, "out_faa"))
#         out_gff = Path(getattr(prodigal, "out_faa"))
#         out_fna = Path(getattr(prodigal, "out_faa"))
#         prodigal_files = {
#             "faa": str(out_faa.resolve()) if out_faa.exists() else "",
#             "gff": str(out_gff.resolve()) if out_gff.exists() else "",
#             "fna": str(out_fna.resolve()) if out_fna.exists() else ""
#         }
#     except Exception:
#         prodigal_files = {
#             "faa": [str(p) for p in outdir_p.glob("*.faa")],
#             "gff": [str(p) for p in outdir_p.glob("*.gff")],
#             "fna": [str(p) for p in outdir_p.glob("*.fna")]
#         }

#     meta.update({
#         "stage": "prodigal_done",
#         "n_proteins": int(len(df_proteins)),
#         "prodigal_outputs": prodigal_files,
#     })
#     (outdir_p / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

#     # ========= Hmmsearch ============================================================
#     print(f"[PIPELINE] 3/5 - HMMER + parse ({e_value_threshold}) ")





#     return fasta_p, outdir_p, model_dir_p, vpf_dict_p, fasta_p