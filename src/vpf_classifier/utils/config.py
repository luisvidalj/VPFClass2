# config.py
from pathlib import Path
import os

# --------------------------------------------
# Current MSL release (ej: MSL35, MSL38, MSL40...)
# --------------------------------------------
ICTV_RELEASE = os.environ.get("ICTV_RELEASE", "MSL33")

# --------------------------------------------
# Base directory
# --------------------------------------------
ROOT_DIR = Path(__file__).resolve().parents[3]         # VPFTaxoClass/
DATA_DIR = ROOT_DIR / "data"             # data/
SCRIPTS_DIR = ROOT_DIR / "scripts"


# ### SI QUEREMOS USAR LAS 227K -> USAR ESTE ####
# class Files:
#     ICTV_RELEASE = ICTV_RELEASE
#     DATA_DIR = DATA_DIR

#     # (GenBank) FASTA file
#     FASTA = DATA_DIR / f"raw/ICTV_databases/Fastas/{ICTV_RELEASE}.fasta"

#     # Prodigal output directory
#     PRODIGAL = DATA_DIR / f"raw/prodigal_gv_output/{ICTV_RELEASE}"
#     FAA = PRODIGAL / f"{ICTV_RELEASE}.faa"
#     GFF = PRODIGAL / f"{ICTV_RELEASE}.gff"

#     # HMM Models (genomad hmm models)

#     #HMM_MODELS = DATA_DIR / "vpf_models/profiles.hmms"
#     HMM_MODELS = Path("/home/uib/proyectos/VPFTaxoClass/tool_data/complete_markers/vpf_data/profiles.hmms")
#     HMM_DICT = Path("/home/uib/proyectos/VPFTaxoClass/tool_data/complete_markers/vpf_data/vpf_to_index.json")
#     #HMM_DICT = DATA_DIR / "vpf_models/vpf_to_index.json"



#     # HMM outputs
#     #HMM_OUTPUT = DATA_DIR / f"raw/output_hmm_join/{ICTV_RELEASE}"
#     HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output_complete/{ICTV_RELEASE}"
#     #HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output_V/{ICTV_RELEASE}"

#     # Si queremos usar las vpfs antiguas
#     # HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output_old/{ICTV_RELEASE}"

#     HMM_TBL = HMM_OUTPUT_MULTIPLE  # Ruta donde buscará archivos .tbl
    

#     # ICTV taxonomic CSV file (opcional si tienes uno por release)
#     ML = DATA_DIR / f'raw/ICTV_databases/MSLs/{ICTV_RELEASE}.csv'


### SI QUEREMOS USAR LAS 160_VV -> USAR ESTA: ####
class Files:
    ICTV_RELEASE = ICTV_RELEASE
    DATA_DIR = DATA_DIR

    # (GenBank) FASTA file
    FASTA = DATA_DIR / f"raw/ICTV_databases/Fastas/{ICTV_RELEASE}.fasta"

    # Prodigal output directory
    PRODIGAL = DATA_DIR / f"raw/prodigal_gv_output/{ICTV_RELEASE}"
    FAA = PRODIGAL / f"{ICTV_RELEASE}.faa"
    GFF = PRODIGAL / f"{ICTV_RELEASE}.gff"

    # HMM Models (genomad hmm models)

    HMM_MODELS = Path("/home/uib/proyectos/VPFTaxoClass/tool_data/virus_markers/vpf_data/profiles_virus.hmms")
    HMM_DICT = Path("/home/uib/proyectos/VPFTaxoClass/tool_data/virus_markers/vpf_data/vpf_to_index_V.json")
    #HMM_MODELS = DATA_DIR / "vpf_models/profiles_virus.hmms"
    #HMM_DICT = DATA_DIR / "vpf_models/vpf_to_index_V.json"


    # HMM outputs
    HMM_OUTPUT = DATA_DIR / f"raw/output_hmm_join/{ICTV_RELEASE}"
    #HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output/{ICTV_RELEASE}"
    HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmmsearch3_ouput/{ICTV_RELEASE}"

    # Si queremos usar las vpfs antiguas
    # HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output_old/{ICTV_RELEASE}"

    HMM_TBL = HMM_OUTPUT_MULTIPLE  # Ruta donde buscará archivos .tbl
    

    # ICTV taxonomic CSV file (opcional si tienes uno por release)
    ML = DATA_DIR / f'raw/ICTV_databases/MSLs/{ICTV_RELEASE}.csv'

# ### SOLO PARA COMPARACION VIEJA 25K: ####
# class Files:
#     ICTV_RELEASE = ICTV_RELEASE
#     DATA_DIR = DATA_DIR

#     # (GenBank) FASTA file
#     FASTA = DATA_DIR / f"raw/ICTV_databases/Fastas/{ICTV_RELEASE}.fasta"

#     # Prodigal output directory
#     PRODIGAL = DATA_DIR / f"raw/prodigal_gv_output/{ICTV_RELEASE}"
#     FAA = PRODIGAL / f"{ICTV_RELEASE}.faa"
#     GFF = PRODIGAL / f"{ICTV_RELEASE}.gff"

#     # HMM Models (genomad hmm models)

#     HMM_MODELS = Path("/home/uib/proyectos/VPFTaxoClass/data/raw/hmm_25k/vpf_data/final_list.hmms")
#     HMM_DICT = Path("/home/uib/proyectos/VPFTaxoClass/data/raw/hmm_25k/vpf_data/vpf_to_index_25k.json")
#     #HMM_MODELS = DATA_DIR / "vpf_models/profiles_virus.hmms"
#     #HMM_DICT = DATA_DIR / "vpf_models/vpf_to_index_V.json"


#     # HMM outputs
#     #HMM_OUTPUT = DATA_DIR / f"raw/output_hmm_join/{ICTV_RELEASE}"
#     #HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output/{ICTV_RELEASE}"
#     HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_25k/output/"

#     # Si queremos usar las vpfs antiguas
#     # HMM_OUTPUT_MULTIPLE = DATA_DIR / f"raw/hmm_output_old/{ICTV_RELEASE}"

#     HMM_TBL = HMM_OUTPUT_MULTIPLE  # Ruta donde buscará archivos .tbl
    

#     # ICTV taxonomic CSV file (opcional si tienes uno por release)
#     ML = DATA_DIR / f'raw/ICTV_databases/MSLs/{ICTV_RELEASE}.csv'


class Constants:
    e_value_threshold = 1e-3
