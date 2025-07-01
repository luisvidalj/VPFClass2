# config.py
from pathlib import Path
import os

# --------------------------------------------
# Current MSL release (ej: MSL35, MSL38, MSL40...)
# --------------------------------------------
ICTV_RELEASE = os.environ.get("ICTV_RELEASE", "mini_test")

# --------------------------------------------
# Base directory
# --------------------------------------------
ROOT_DIR = Path(__file__).resolve().parents[3]         # VPFTaxoClass/
DATA_DIR = ROOT_DIR / "data"             # data/
SCRIPTS_DIR = ROOT_DIR / "scripts"

class Files:
    DATA_DIR = DATA_DIR

    # (GenBank) FASTA file
    FASTA = DATA_DIR / f"raw/ICTV_databases/{ICTV_RELEASE}.fasta"

    # Prodigal output directory
    PRODIGAL = DATA_DIR / f"raw/prodigal_output/{ICTV_RELEASE}"

    # HMM Models (genomad hmm models)
    HMM_MODELS = DATA_DIR / "vpf_models/profiles.hmms"


    # HMM outputs
    HMM_OUTPUT = DATA_DIR / f"raw/output_hmm_join/{ICTV_RELEASE}"
    HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output/{ICTV_RELEASE}"
    HMM_TBL = HMM_OUTPUT_MULTIPLE  # Ruta donde buscará archivos .tbl

    # ICTV taxonomic CSV file (opcional si tienes uno por release)
    ML = DATA_DIR / f'ICTV_databases/MSLs/{ICTV_RELEASE}.csv'



class Constants:
    prot_embedding_esm = 'esm2_t6_8M'
    graph_vpf_threshold = 0.2
    e_value_threshold = 1e-3
