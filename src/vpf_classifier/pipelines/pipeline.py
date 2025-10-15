import pandas as pd
from vpf_classifier.parsers.fasta_parser import FastaParser
from vpf_classifier.parsers.prodigal_parser import Prodigal
from vpf_classifier.parsers.vpf_parser import VPF_parser
from vpf_classifier.parsers.ictv_parser import clean_ictv_csv,merge_vpf_with_ictv
from vpf_classifier.utils.config import Files



def run_pipeline(e_value_threshold=1e-3, num_cpus=12, vector_norm=None):

    import shutil

    print("[DEBUG] Prodigal path:", shutil.which("prodigal"))
    print("[DEBUG] HMMER path:", shutil.which("hmmsearch"))

    print(f"[PIPELINE] Starting run for ICTV release :{Files.ICTV_RELEASE}")
    print("[STEP 1] Pasing FASTA file...")
    fasta_parser = FastaParser(fna_path=Files.FASTA)
    fasta_parser.parse_fasta_to_dataframe()

    print("[STEP 2] Running Prodigal (or reusing previous output)...")
    prodigal = Prodigal(parser=fasta_parser,output_dir=Files.PRODIGAL)
    df_proteins = prodigal.parse_prodigal()

    print("[STEP 3] Running HMMER (if needed) and parsing .tbl hits...")
    vpf = VPF_parser(parser=fasta_parser,hmm_file=Files.HMM_MODELS ,e_value_threshold=e_value_threshold, num_cpus=num_cpus,vector_norm=vector_norm)
    vpf.parse_multiple_hmm_parallel()

    print(f"[STEP 4] Merging VPF hits with taxonomic metadata from ICTV ({Files.ICTV_RELEASE})")
    ictv_df = pd.read_csv(Files.ML, sep=";")
    expanded_ictv = clean_ictv_csv(ictv_df=ictv_df, msl_tag=Files.ICTV_RELEASE)
    merged_df = merge_vpf_with_ictv(vpf_database=vpf.df_virus_hmm, ictv_df=expanded_ictv)

    print('[DONE] 4 datasets are outputed. The last one contains the metadata from ICTV')

    return fasta_parser.ncbi_df, df_proteins, vpf.df_virus_hmm, merged_df


