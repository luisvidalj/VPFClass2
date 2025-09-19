# vpf_classifier/pipelines/infer_pipelines.py

from __future__ import annotations
import os
import json
import time
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz

from vpf_classifier.parsers.fasta_parser import FastaParser
from vpf_classifier.parsers.prodigal_parser import Prodigal
from vpf_classifier.parsers.vpf_parser import VPF_parser
from vpf_classifier.utils.runtime_paths import RuntimePaths


#### ------------------- Utility functions ------------------------------------

# keys:str
# values: Any
def _write_json(path: Path, obj: Dict[str, Any]) -> None:
    """ Guardar un dict como JSON con identacion y UTF-8"""
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False))

def _tool_versions() -> Dict[str, Optional[str]]:
    """
    Intenta capturar versiones de prodigal/hmmsearch para el manifest.
    Si algo falla (no están instalados, etc.), devuelve None.
    """
    import subprocess
    def ver(cmd: str, args: List[str]) -> Optional[str]:
        try:
            out = subprocess.check_output([cmd, *args], stderr=subprocess.STDOUT, text=True, timeout=10)
            return out.strip().splitlines()[0][:200]
        except Exception:
            return None

    return {
        "prodigal": ver("prodigal", ["-v"]),
        "hmmsearch": ver("hmmsearch", ["-h"]),
    }

def _save_features(outdir: Path, X: csr_matrix, accessions: List[str]) -> None:
    """
    Guarda la matriz dispersa de features y el orden de accessions.
    """
    save_npz(outdir / "features.npz", X.tocsr())
    (outdir / "features.index.csv").write_text("\n".join(accessions))


#### ------------------- Predictions placeholder ------------------------------------ (c/p)

def _empty_predictions(df_virus_hmm: pd.DataFrame, reason: str = "NO_MODEL_BUNDLE") -> pd.DataFrame:
    """
    Crea un DataFrame de predicciones 'placeholder' (hasta integrar modelo).
    Incluye contadores útiles: nº de HMM hits y nº de proteínas con hits.
    """
    n_hits = df_virus_hmm.get("hmms_hits", [[]]).apply(
        lambda x: len(x) if isinstance(x, (list, tuple)) else 0
    )
    n_prot = df_virus_hmm.get("protein_accessions", [[]]).apply(
        lambda x: len(set(x)) if isinstance(x, (list, tuple)) else 0
    )

    return pd.DataFrame({
        "Accession": df_virus_hmm["Accession"],
        "pred_genus": "Unclassified",
        "pred_score": np.nan,
        "topk_genus": "",
        "status": reason,
        "n_vpf_hits": n_hits,
        "n_proteins_with_hits": n_prot,
    })

#### --- Model bundle loading utilities -----------------------------------------
from dataclasses import dataclass

# --- Model bundle loading (shared vpf_to_index.json lives outside the bundle) --

from dataclasses import dataclass

@dataclass
class ModelBundle:
    bundle_dir: Path
    target_level: str                 # "genus" | "family"
    config: Dict[str, Any]            # config.json
    idx_to_label: Dict[int, str]      # índice -> nombre clase
    feature_norm: str                 # "none" | "tfidf" | "standardize"
    feature_params: Dict[str, Any]    # parámetros del normalizador
    taxonomy_table: Optional[pd.DataFrame]
    model: Any                        # torch.nn.Module cargado
    device: str                       # "cpu"/"cuda"
    vpf_dict_sha256_expected: Optional[str]  # hash opcional que guarda el modelo

def _resolve_model_dir(model_id: Optional[str], bundle_dir: Optional[str]) -> Path:
    if bundle_dir:
        p = Path(bundle_dir).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"[ERROR] bundle_dir not found: {p}")
        return p

    if model_id:
        candidates = [
            Path("models") / model_id,
            Path("models") / model_id.lower(),
        ]
        if "_" in model_id:
            target, msl = model_id.split("_", 1)
            candidates += [
                Path("models") / target / msl,
                Path("models") / target.lower() / msl.upper(),
            ]
        for c in candidates:
            if c.exists():
                return c.resolve()
        raise FileNotFoundError(f"[ERROR] Could not resolve model_id '{model_id}' under ./models/")
    raise ValueError("[ERROR] Provide either bundle_dir or model_id.")

def _safe_read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] Missing JSON file: {path}")
    return json.loads(path.read_text())

def load_model_bundle(
    model_id: Optional[str] = None,
    bundle_dir: Optional[str] = None,
    force_device: str = "auto"
) -> ModelBundle:
    """
    Carga el bundle de modelo (PERO el vpf_to_index.json NO va aquí; es global y se pasa por CLI/RuntimePaths).
    config.json puede incluir 'vpf_dict_sha256' para comprobar compatibilidad con el archivo global.
    """
    bdir = _resolve_model_dir(model_id, bundle_dir)

    cfg = _safe_read_json(bdir / "config.json")
    idx_to_label = _safe_read_json(bdir / "idx_to_label.json")

    feature_norm = cfg.get("feature_norm", "none").lower()
    feature_params = cfg.get("feature_params", {})
    target_level  = cfg.get("target_level", "genus").lower()
    vpf_sha_expected = cfg.get("vpf_dict_sha256", None)

    tax_path = bdir / "taxonomy_table.csv"
    taxonomy_table = pd.read_csv(tax_path) if tax_path.exists() else None

    import torch
    device = "cuda" if (force_device in ("auto", "cuda") and torch.cuda.is_available()) else "cpu"
    model_path = bdir / "model.pt"
    if not model_path.exists():
        raise FileNotFoundError(f"[ERROR] Missing model.pt in bundle: {model_path}")

    model = torch.load(model_path, map_location=device)
    if hasattr(model, "eval"):
        model.eval()

    return ModelBundle(
        bundle_dir=bdir,
        target_level=target_level,
        config=cfg,
        idx_to_label={int(k): v for k, v in idx_to_label.items()},
        feature_norm=feature_norm,
        feature_params=feature_params,
        taxonomy_table=taxonomy_table,
        model=model,
        device=device,
        vpf_dict_sha256_expected=vpf_sha_expected
    )
# --- VPF dict compatibility check -------------------------------------------

import hashlib

def _sha256_of_file(path: Path) -> str:
    """Devuelve el SHA-256 hex de un archivo sin cargarlo entero en memoria."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):  # 1MB
            h.update(chunk)
    return h.hexdigest()

def _check_vpf_dict_compat(vpf_dict_path: Path, bundle: ModelBundle) -> None:
    """
    Si el bundle declara un 'vpf_dict_sha256' en su config, comprobamos que el archivo global
    vpf_to_index.json que está usando el usuario coincide. Es solo un WARNING si no coincide.
    """
    if bundle.vpf_dict_sha256_expected:
        actual = _sha256_of_file(vpf_dict_path)
        if actual != bundle.vpf_dict_sha256_expected:
            print(f"[WARNING] vpf_to_index.json hash mismatch.\n"
                  f"  expected: {bundle.vpf_dict_sha256_expected}\n"
                  f"  actual:   {actual}\n"
                  f"  (Proceeding anyway; ensure you are using the right mapping.)")



#### ------------------- Funcion principal ------------------------------------ (c/p)
# Los argumentos de despues del * -> Solo se pueden pasasr con el nombre (no de manera posicional)

def classify_fasta(
    fna_path: str,
    out_dir: Optional[str] = None,
    *,
    e_value_threshold: float = 1e-3,
    num_cpus: int = 8,
    hmm_models: Optional[str] = None,
    vpf_dict: Optional[str] = None,
    save_features: bool = True,
) -> Tuple[pd.DataFrame, Path]:
    """
    Pipeline de inferencia (esqueleto, sin cargar modelo aún).

    Flujo:
      1) FASTA -> DataFrame de contigs (FastaParser)
      2) Prodigal -> proteínas .faa/.gff (se reutiliza si existen) (Prodigal)
      3) HMMER -> parseo de .tbl y construcción de features VPF (VPF_parser)
      4) Predicciones placeholder (sin modelo)
      5) Guardado de outputs: predictions.csv + run_manifest.json (+ features opcionales)

    Returns:
        (predictions_df, outdir_path)
    """
    run_start = time.time()

    # --- 0) Resolver rutas para este run (no dependemos de Files aquí) ---
    #     - Cada run crea su carpeta con subcarpetas para prodigal y hmmer_tbl
    paths = RuntimePaths.resolve(
        fasta=fna_path,
        outdir=out_dir,
        hmm_models=hmm_models,
        vpf_dict=vpf_dict,
        files_fallback=None  # inferencia: no usamos Files a menos que tú lo decidas más adelante
    )
    outdir = paths.outdir
    logs_path = outdir / "logs.txt"

    # Pequeño log con parámetros y versiones de herramientas
    with logs_path.open("w") as lg:
        lg.write(f"[PIPELINE] Inference run started: {time.ctime(run_start)}\n")
        lg.write(f"FASTA: {paths.fasta}\n")
        lg.write(f"OUTDIR: {paths.outdir}\n")
        lg.write(f"Prodigal dir: {paths.prodigal_dir}\n")
        lg.write(f"HMM out dir: {paths.hmm_out_dir}\n")
        lg.write(f"hmm_models: {paths.hmm_models}\n")
        lg.write(f"vpf_dict: {paths.vpf_dict}\n")
        lg.write(f"e_value_threshold={e_value_threshold}, num_cpus={num_cpus}\n")
        lg.write(f"tool_versions={_tool_versions()}\n")

    # --- 1) FASTA parsing ---
    fasta_parser = FastaParser(fna_path=paths.fasta)
    ncbi_df = fasta_parser.parse_fasta_to_dataframe(return_df=True)

    # --- 2) Prodigal (idempotente) ---
    #     Usamos la carpeta específica del run para .faa/.gff
    prod = Prodigal(parser=fasta_parser, output_dir=paths.prodigal_dir)
    df_virus_prot = prod.parse_prodigal()

    # 3) Cargar bundle (elige por bundle_dir o model_id)
    bundle = load_model_bundle(model_id=model_id, bundle_dir=bundle_dir, force_device=device)

    # 3.b) (nuevo) Chequear compatibilidad del vpf_to_index.json GLOBAL con el bundle
    _check_vpf_dict_compat(paths.vpf_dict, bundle)

    # --- 4) HMMER -> VPF features ---
    #     Enviamos:
    #       - hmm_file=paths.hmm_models (para que use el .hmms que tú elijas)
    #       - vpf_dict_path=paths.vpf_dict (diccionario VPF alineado con el modelo futuro)
    #       - hmm_output_dir=paths.hmm_out_dir (los .tbl del run)
    vpf = VPF_parser(
        parser=fasta_parser,
        hmm_file=paths.hmm_models,
        e_value_threshold=e_value_threshold,
        num_cpus=num_cpus,
        vpf_dict_path=paths.vpf_dict,
        hmm_output_dir=paths.hmm_out_dir,
    )
    vpf.parse_multiple_hmm()
    df_virus_hmm = vpf.df_virus_hmm.copy()

    # Guardar features dispersas e índice de accessions (útil para depurar y para la inferencia real después)
    if save_features and hasattr(vpf, "vpf_sparse_matrix"):
        accessions = df_virus_hmm["Accession"].tolist()
        _save_features(outdir, vpf.vpf_sparse_matrix, accessions)

    # --- 4) Predicciones placeholder ---
    preds_df = _empty_predictions(df_virus_hmm, reason="NO_MODEL_BUNDLE")
    preds_path = outdir / "predictions.csv"
    preds_df.to_csv(preds_path, index=False)

    # --- 5) Manifest del run ---
    manifest = {
        "started_at": time.ctime(run_start),
        "finished_at": time.ctime(),
        "duration_sec": round(time.time() - run_start, 2),
        "n_contigs": int(ncbi_df.shape[0]),
        "n_with_features": int(df_virus_hmm.shape[0]),
        "params": {
            "e_value_threshold": e_value_threshold,
            "num_cpus": num_cpus,
            "hmm_models": str(paths.hmm_models),
            "vpf_dict": str(paths.vpf_dict),
        },
        "tools": _tool_versions(),
        "outputs": {
            "predictions_csv": str(preds_path),
            "features_npz": str(outdir / "features.npz") if save_features else None,
            "features_index": str(outdir / "features.index.csv") if save_features else None,
            "logs": str(logs_path),
            "prodigal_dir": str(paths.prodigal_dir),
            "hmmer_tbl_dir": str(paths.hmm_out_dir),
        },
        "notes": "Model bundle not integrated yet; predictions are placeholders.",
    }
    _write_json(outdir / "run_manifest.json", manifest)

    print(f"[DONE] Inference skeleton finished. Outputs at: {outdir}")
    return preds_df, outdir
