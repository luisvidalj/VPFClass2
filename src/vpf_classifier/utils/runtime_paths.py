from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import time

# @dataclass
# class RuntimePaths:
#     """
#     Minimal runtime path resolver for user-facing runs.
#     For now we only need:
#         - fasta: input FASTA path (resolved to absolute)
#         - outdir: output directory (auto-created)
#     """
#     fasta: Path
#     outdir: Path

#     @classmethod
#     def resolve(cls, fasta: str, outdir: str | None = None) -> "RuntimePaths":
#         # 1) Normalize Fasta path to absolute

         
#         fasta_p = Path(fasta).expanduser().resolve()

#         # 2) Create an outdir if not passed
#         if outdir:
#             outdir_p = Path(outdir).expanduser().resolve()
#         else:
#             stamp = time.strftime("%Y%m%d_%H%M%S")
#             outdir_p = Path("outputs") / f"{stamp}_run"

#         # 3) Create the output folder if it does not exist
#         outdir_p.mkdir(parents=True, exist_ok=True)

#         return cls(fasta=fasta_p, outdir=outdir_p)
           


from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import os, time

@dataclass
class RuntimePaths:
    fasta: Path
    outdir: Path
    prodigal_dir: Path
    hmm_out_dir: Path
    hmm_models: Path
    vpf_dict: Path

    @classmethod
    def resolve(
        cls,
        fasta: str,
        outdir: str | None = None,
        hmm_models: str | None = None,
        vpf_dict: str | None = None,
        *,
        files_fallback=None  # opcional: pasa vpf_classifier.utils.config.Files como respaldo
    ) -> "RuntimePaths":
        fasta_p = Path(fasta).expanduser().resolve()
        run_dir = Path(outdir).expanduser().resolve() if outdir else Path("outputs") / (time.strftime("%Y%m%d_%H%M%S") + "_run")
        prodigal_dir = run_dir / "prodigal"
        hmm_out_dir  = run_dir / "hmmer_tbl"

        # Prioridad: argumentos > variables de entorno > fallback Files > defaults del paquete
        hmm_models_p = (
            Path(hmm_models).expanduser().resolve() if hmm_models else
            (Path(os.environ["VPFCLASS_HMM_MODELS"]).expanduser().resolve() if "VPFCLASS_HMM_MODELS" in os.environ else
             (files_fallback.HMM_MODELS if files_fallback else Path("models") / "markers.hmms"))
        )
        vpf_dict_p = (
            Path(vpf_dict).expanduser().resolve() if vpf_dict else
            (Path(os.environ["VPFCLASS_VPF_DICT"]).expanduser().resolve() if "VPFCLASS_VPF_DICT" in os.environ else
             (files_fallback.HMM_DICT if files_fallback else Path("models") / "vpf_to_index.json"))
        )

        for d in (run_dir, prodigal_dir, hmm_out_dir):
            d.mkdir(parents=True, exist_ok=True)

        return cls(
            fasta=fasta_p,
            outdir=run_dir,
            prodigal_dir=prodigal_dir,
            hmm_out_dir=hmm_out_dir,
            hmm_models=hmm_models_p,
            vpf_dict=vpf_dict_p,
        )
