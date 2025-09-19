"""
Idea del cli.py: Intenta trobar rutes amb els arxius corresponents seguint la següent prioritat:
1. CLI -> Si l'usuari empra les flags corresponents i especifica les rutes demanades, s'empraran
2. Variables d'entorn -> Si l'usuari defineix variables d'entorn (explicarem com), no fara falta indiqui flags
3. Autodeteccio al repositori -> pedr defecte hauria de cercar unes rutes concretes. Si l'usuari ha seguit les instruccions del readme
i no ha indicat cap de les altres opcions, tocaria funcionar (aquesta es la que m'agrada mes)
"""

import argparse
import os, shutil, sys
from pathlib import Path

from vpf_classifier.pipelines.vpf_class2 import run_user_pipeline


############################ Helpers 2 ######################################################

from pathlib import Path
import os

def _first_existing(*cands: Path | None) -> Path | None:
    for c in cands:
        if c and c.exists():
            return c
    return None

def _find_tool_data_root(given: str | None) -> Path:
    """
    1) CLI --tool-data
    2) ENV VPFCLASS_TOOL_DATA
    3) Autodetect: 'tool_data'
    """
    # 0) given by line command
    if given:
        print("GIVEN",given)
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
    """Acepta 'vpf_data' (preferida) o 'vpf_models' (legacy)."""
    for name in ("vpf_data"):
        cand = markers_root / name
        if cand.is_dir():
            return cand
    return markers_root / "vpf_data"

def resolve_from_markers(tool_data_root: Path, markers: str, msl: str):
    """
    tool_data/{complete_markers|virus_markers}/models/MSL{msl}
    tool_data/{...}/{vpf_data}/{vpf_to_index.json, profiles.hmms}
    """
    base = {"all": "complete_markers", "virus": "virus_markers"}[markers]
    root = tool_data_root / base
    print("BUSCO EN: ", root)
    model_dir  = root / "models" / f"MSL{msl}"
    print("MODELOS EN: ", model_dir)
    vpf_dir    = _vpf_dir_for(root)
    print("VPF DIR: ", vpf_dir)

    json_files = list(vpf_dir.glob("*.json"))
    hmms_files = list(vpf_dir.glob("*.hmms"))
    vpf_dict = json_files[0] if json_files else None
    hmm_models = hmms_files[0] if hmms_files else None

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
##### 1) CODIGO CUANDO PERMITIAMOS RESOLVER POR RUTA ESPECIFICADA ###################################
    # md_cli = Path(model_dir).expanduser().resolve() if model_dir else None
    # vd_cli = Path(vpf_dict).expanduser().resolve() if vpf_dict else None
    # hm_cli = Path(hmm_models).expanduser().resolve() if hmm_models else None
    # if md_cli or vd_cli or hm_cli:
    #     # Si el usuario pasó rutas explícitas, las usamos tal cual (las otras pueden resolverse con None)
    #     return md_cli, vd_cli, hm_cli, _find_tool_data_root(tool_data)
##################################################################################################

    # Finding tool_data directory and solving paths based on: markers + msl
    td_root = _find_tool_data_root(tool_data)
    print(td_root)
    md_mark, vd_mark, hm_mark = resolve_from_markers(td_root, markers, msl)
    return md_mark, vd_mark, hm_mark, td_root



#############################################################################################


# AJUSTAR/QUITAR
# def _default_repo_paths() -> dict:
#     """
#     Crea rutes "default". Es a dir, la comanda sense especificar
#     les rutes, les cercara aqui.
#     """
#     # Aixo crea una llista de rutes [cwd (pwd), pare, abuelo, ...]
#     cwd = Path.cwd()
#     roots = [cwd] + list(cwd.parents) 

#     # i lo que volem despres es trobar sa ruta mes propera que conte data i models (tocaria si han seguit ses instruccions)
#     # de manera que s'arrel des repo sera d'on pengin ses carpetes data/ i models/ (pot ser enlloc de .exists() tocaria emprar .is_dir()
#     # ja que volem que siguin carpetes)
#     repo_root = next((r for r in roots if (r / "data").exists() and (r / "models").exists()), cwd)

#     return {
#         "model_dir": repo_root / "models" / "complete_markers" / "Genus_MSL40",
#         "vpf_dict": repo_root / "data" / "vpf_models" / "vpf_to_index.json",
#         "hmm_models": repo_root / "data" / "vpf_models" / "profiles.hmms"
#     }

# def _resolver_or_default(model_dir: str | None, vpf_dict: str | None, hmm_models: str | None):
    
#     d = _default_repo_paths()

#     # Si s'han introdut rutes -> Seran la primera opcio
#     md_cli = Path(model_dir).expanduser().resolve() if model_dir else None
#     vd_cli = Path(vpf_dict).expanduser().resolve() if vpf_dict else None
#     hm_cli = Path(hmm_models).expanduser().resolve() if hmm_models else None 

#     # Si s'han definit variables d'entorn (segona opcio)
#     md_env = Path(os.getenv("VPFCLASS_MODELS_DIR")).expanduser().resolve() if os.getenv("VPFCLASS_MODELS_DIR") else None
#     vd_env = Path(os.getenv("VPFCLASS_VPF_DICT")).expanduser().resolve() if os.getenv("VPFCLASS_VPF_DICT") else None
#     hm_env = Path(os.getenv("VPFCLASS_HMM_MODELS")).expanduser().resolve() if os.getenv("VPFCLASS_HMM_MODELS") else None

#     # Els default
#     md_def, vd_def, hm_def = d["model_dir"], d["vpf_dict"], d["hmm_models"]

#     def _first_existing(*cands):
#         for c in cands:
#             if c and c.exists():
#                 return c
#         return None
    
#     # se podria deixar mes estricte dividint en file i dir 
#     return (
#         _first_existing(md_cli, md_env, md_def),
#         _first_existing(vd_cli, vd_env, vd_def),
#         _first_existing(hm_cli, hm_env, hm_def)
#     )




# def _resolve_path(given: str | None, env_var: str | None, default: str| None) -> Path | None:
#     cand = given or (os.getenv(env_var) if env_var else None) or default
#     return Path(cand).expanduser().resolve() if cand else None

def main():
    p = argparse.ArgumentParser(prog="vpfclass2", description="VPF-Class 2 command line")
    sub = p.add_subparsers(dest="cmd", required=True) # crea subcomandos (args.cmd)

    # vpf-class predictions
    pred = sub.add_parser("predict", help="Taxonomic prediction of viral contigs")
    pred.add_argument("--fasta", required=True, help="Input FASTA/contigs (.fna/.fa/.fasta)")
    pred.add_argument("--outdir", default=None, help="Output directory")
    # pred.add_argument("--model-dir", default=None, help="Model bundle directory (model.pt, idx_to_label.json)")
    # pred.add_argument("--vpf-dict", default=None, help="vpf_to_index.json")
    # pred.add_argument("--hmm-models", default=None, help="profiles.hmms")
    pred.add_argument("--num-cpus", type=int, default=8)
    pred.add_argument("--device", choices=["cpu","cuda"], default="cuda")
    pred.add_argument("--topk", type=int, default=3)
    pred.add_argument("--batch-size", type=int, default=1024, help="Batch for inference (optional, if exposed)")
    
    pred.add_argument("--markers", choices=["all","virus"], default="all",
                      help="Choose 'all' to use the 227k profile markers or 'virus' to use the viral specific ones. Default 'all'")
    pred.add_argument("--msl", choices=["40", "39", "37", "35"], default="40",
                  help="MSL Release the model has been trained with. Availables: 40,39,37,35")
    pred.add_argument("--tool-data", default=None,
                  help="Root of the data bundle (if not trying to allocate at 'tool_data/').")
    


    chk = sub.add_parser("check", help="Validate environment and resources")
    # chk.add_argument("--model-dir", default=None)
    # chk.add_argument("--vpf-dict", default=None)
    # chk.add_argument("--hmm-models", default=None)
    chk.add_argument("--markers", choices=["all", "virus"], default="all")
    chk.add_argument("--msl", choices=["40", "39", "37", "35"], default="40")
    chk.add_argument("--tool-data", default=None)

    args = p.parse_args()

    if args.cmd == "check":
        ok = True
        print("[CHECK] Prodgial:", shutil.which("prodigal-gv") or "NOT FOUND"); ok &= bool(shutil.which("prodigal-gv"))
        print("[CHECK] hmmsearch:", shutil.which("hmmsearch") or "NOT FOUND"); ok &= bool(shutil.which("hmmsearch"))

        # md, vd, hm = _resolver_or_default(args.model_dir, args.vpf_dict, args.hmm_models)
        md, vd, hm, td = _resolve_resources(
            markers=args.markers, msl=args.msl, tool_data=args.tool_data
        )


        if md:
            ok_md = (md.exists() and (md/"model.pt").exists() and (md/"idx_to_label.json").exists())
            print("[CHECK] model_dir:", md, "OK" if ok_md else "MISSING: model.pt or and idx_to_label.json")
        else:
            ok_md = False
            print("[CHECK] model_dir: not found (use --model-dir or defaults)")

        if vd:
            ok_vd = vd.exists()
            print("[CHECK] vpf_to_index.json:", vd, "OK" if ok_vd else "MISSING: vpf_to_index.json")
        else:
            ok_vd = False
            print("[CHECK] vpf_to_index.json: not found (use --vpf-dict or defaults)")

        if hm:
            ok_hm = hm.exists()
            print("[CHECK] profiles.hmms:", hm, "OK" if ok_hm else "MISSING")
        else:
            ok_hm = False
            print("[CHECK] profiles.hmms: not found (use --hmm-models or defaults)")

        sys.exit(0 if (ok_md and ok_vd and ok_hm) else 1)

    elif args.cmd == "predict":
        # md, vd, hm = _resolver_or_default(args.model_dir, args.vpf_dict, args.hmm_models)
        md, vd, hm, td = _resolve_resources(
            markers=args.markers, msl=args.msl, tool_data=args.tool_data
        )
        if not md or not vd or not hm:
            print("[ERROR] Could not resolve required resources. Try 'vpfclass check' or pass flags explicitly.")

        # Pasamos al pipeline lo que ya resolvimos
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