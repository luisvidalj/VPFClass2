"""
CLI entrypoint for VPF-Class 2.

Resource resolution priority:
1) CLI flags (e.g., --tool-data, --markers, --msl)
2) Environment variables (e.g., VPFCLASS_TOOL_DATA)
3) Repository autodetection (search for a local 'tool_data/' directory)
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

from vpf_classifier.pipelines.vpf_class2 import run_user_pipeline


def _find_tool_data_root(given: str | None) -> Path:
    """
    1) CLI --tool-data
    2) ENV VPFCLASS_TOOL_DATA
    3) Autodetect: 'tool_data'
    """
    # 0) given by line command
    if given:
        return Path(given).expanduser().resolve()
    
    # 1) If the user has specified the path as env variable (Readme explanation)
    env_td = os.getenv("VPFCLASS_TOOL_DATA")
    if env_td:
        return Path(env_td).expanduser().resolve()
    cwd = Path.cwd()

    # 2) Trying to autodetect the 'tool_data' folder
    for r in [cwd] + list(cwd.parents):
        cand = r / "tool_data"
        if cand.is_dir():
            return cand.resolve()
        
    raise FileNotFoundError("tool_dir/ path not found. Make sure you have downloaded the required data")

def _vpf_dir_for(markers_root: Path) -> Path:
    """Return the directory containing VPF resources."""
    cand = markers_root / "vpf_data"
    return cand if cand.is_dir() else cand

def resolve_from_markers(tool_data_root: Path, markers: str, msl: str):
    """
    tool_data/{complete_markers|virus_markers}/models/MSL{msl}
    tool_data/{...}/{vpf_data}/{vpf_to_index.json, profiles.hmms}
    """
    base_map = {"all": "complete_markers", "virus": "virus_markers"}
    if markers not in base_map:
        raise ValueError(f"Unknown markers='{markers}'. Expected one of {list(base_map.keys())}.")
    root = tool_data_root / base_map[markers]
    model_dir  = root / "models" / f"MSL{msl}"
    vpf_dir    = _vpf_dir_for(root)

    if markers == "virus":
        expected_json = vpf_dir / "vpf_to_index_V.json"
        expected_hmms = vpf_dir / "profiles_virus.hmms"
    else: # markers == "all"
        expected_json = vpf_dir / "vpf_to_index.json"
        expected_hmms = vpf_dir / "profiles.hmms"

    

    vpf_dict = expected_json if expected_json.exists() else None
    hmm_models = expected_hmms if expected_hmms.exists() else None


    return model_dir, vpf_dict, hmm_models

def _resolve_resources(
    markers: str,
    msl: str,
    tool_data: str | None,
):
    """
    Prioridad:
      1) --markers/--msl dentro de --tool-data (o env/autodetect) [ESTRICTO]
         -> si faltan archivos, devolvemos igualmente esas rutas (el CHECK dirá MISSING)
      3) (Compat opcional) ENV antiguos y defaults del repo sólo si NO se usó markers/msl
    """
    # Finding tool_data directory and solving paths based on: markers + msl
    td_root = _find_tool_data_root(tool_data)
    md_mark, vd_mark, hm_mark = resolve_from_markers(td_root, markers, msl)
    return md_mark, vd_mark, hm_mark, td_root


def main():
    p = argparse.ArgumentParser(prog="vpfclass2", description="VPF-Class 2 command line")
    sub = p.add_subparsers(dest="cmd", required=True) 

    # vpf-class predictions
    pred = sub.add_parser("predict", help="Taxonomic prediction of viral contigs")
    pred.add_argument("--fasta", required=True, help="Input FASTA/contigs (.fna/.fa/.fasta)")
    pred.add_argument("--outdir", default=None, help="Output directory")
    pred.add_argument("--num-cpus", type=int, default=8)
    pred.add_argument("--device", choices=["cpu","cuda"], default="cuda")
    pred.add_argument("--topk", type=int, default=3)
    pred.add_argument("--batch-size", type=int, default=1024, help="Batch for inference (optional, if exposed)")
    
    pred.add_argument("--markers", choices=["all","virus"], default="virus",
                      help="Choose 'all' to use the 227k profile markers or 'virus' to use the viral specific ones. Default 'virus'")
    pred.add_argument("--msl", choices=["40", "39", "38", "37", "36", "35", "34", "33", "32", "31"], default="40",
                  help="MSL Release the model has been trained with. Availables: 40,39,38,37,36,35,34,33,32,31")
    pred.add_argument("--tool-data", default=None,
                  help="Root of the data bundle (if not trying to allocate at 'tool_data/').")
    


    chk = sub.add_parser("check", help="Validate environment and resources")
    chk.add_argument("--markers", choices=["all", "virus"], default="virus")
    chk.add_argument("--msl", choices=["40", "39", "38", "37", "36", "35", "34", "33", "32", "31"], default="40")
    chk.add_argument("--tool-data", default=None)

    args = p.parse_args()

    if args.cmd == "check":
        ok = True
        print("[CHECK] prodigal-gv:", shutil.which("prodigal-gv") or "NOT FOUND"); ok &= bool(shutil.which("prodigal-gv"))
        print("[CHECK] hmmsearch:", shutil.which("hmmsearch") or "NOT FOUND"); ok &= bool(shutil.which("hmmsearch"))

        md, vd, hm, td = _resolve_resources(
            markers=args.markers, msl=args.msl, tool_data=args.tool_data
        )


        if md:
            ok_genus = (md.exists() and (md / "Genus"/"model.pt").exists() and (md / "Genus" / "idx_to_label.json").exists())
            ok_fam  = (md.exists() and (md / "Family"/"model.pt").exists() and (md / "Family" / "idx_to_label.json").exists())
            ok_md = ok_genus and ok_fam
            print("[CHECK] model_dir:", md, "--> OK" if ok_md else "MISSING: model.pt or and idx_to_label.json")
        else:
            ok_md = False
            print("[CHECK] model_dir: not found (use --model-dir or defaults)")

        if vd:
            ok_vd = vd.exists()
            print("[CHECK] vpf_to_index.json:", vd, "--> OK" if ok_vd else "MISSING: vpf_to_index.json")
        else:
            ok_vd = False
            print("[CHECK] vpf_to_index.json: not found (use --vpf-dict or defaults)")

        if hm:
            ok_hm = hm.exists()
            print("[CHECK] profiles.hmms:", hm, "--> OK" if ok_hm else "MISSING")
        else:
            ok_hm = False
            print("[CHECK] profiles.hmms: not found (use --hmm-models or defaults)")

        sys.exit(0 if (ok_md and ok_vd and ok_hm) else 1)

    elif args.cmd == "predict":
        md, vd, hm, td = _resolve_resources(
            markers=args.markers, msl=args.msl, tool_data=args.tool_data
        )
        if not md or not vd or not hm:
            print("[ERROR] Could not resolve required resources. Try 'vpfclass check' or pass flags explicitly.")
            sys.exit(1)

        run_user_pipeline(
            fasta=args.fasta,
            outdir=args.outdir,
            model_dir=str(md),
            vpf_dict=str(vd),
            hmm_models=str(hm),
            e_value_threshold=1e-3,
            num_cpus=args.num_cpus,
            device=args.device,
            topk=args.topk,
            # extra (si ya expusiste estos parámetros en la pipeline; p.ej. batch_size/hmm_models)
            # batch_size=args.batch_size,
        )

if __name__ == "__main__":
    main()