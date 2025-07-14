from vpf_classifier.parsers.fasta_parser import FastaParser
from vpf_classifier.parsers.prodigal_parser import Prodigal
from vpf_classifier.parsers.vpf_parser import VPF_parser
from vpf_classifier.utils.config import Files


def run_pipeline(e_value_threshold=1e-3, num_cpus=12):
    print(f"[PIPELINE] Starting run for ICTV release :{Files.ICTV_RELEASE}")
    print("[STEP 1] Pasing FASTA file...")
    fasta_parser = FastaParser(fna_path=Files.FASTA)
    fasta_parser.parse_fasta_to_dataframe()

    print("[STEP 2] Running Prodigal (or reusing previous output)...")
    prodigal = Prodigal(parser=fasta_parser)
    df_proteins = prodigal.parse_prodigal()

    print("[STEP 3] Running HMMER (if needed) and parsing .tbl hits...")
    vpf = VPF_parser(parser=fasta_parser, e_value_threshold=e_value_threshold, num_cpus=num_cpus)
    vpf.parse_multiple_hmm()

    print("[DONE] VPF data available in .df_virus_hmm")
    return fasta_parser.ncbi_df, df_proteins, vpf.df_virus_hmm


