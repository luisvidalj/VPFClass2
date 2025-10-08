# src/ictv_classifier/parsers/refseq_parser.py

from Bio import SeqIO
import pandas as pd
from pathlib import Path
import os
import glob
from vpf_classifier.utils.config import Files
from vpf_classifier.utils.config import SCRIPTS_DIR
from vpf_classifier.utils.config import DATA_DIR
from typing import Optional
import subprocess


class FastaParser:
    def __init__(self, fna_path: Optional[Path] = None):
        if fna_path is not None:
            self.fna_path = fna_path
        else:
            print("ME CAGO EN MI SANTA VIDA")
            self.fna_path = Files.FASTA
        self.faa = None
        self.ncbi_df = None

    def parse_fasta_to_dataframe(self, return_df=True) -> pd.DataFrame:
        records = list(SeqIO.parse(self.fna_path, "fasta"))
        data = []

        for record in records:
            accession = record.id
            description = record.description.replace(accession, "").strip()
            sequence = str(record.seq)
            length = len(sequence)

            data.append({
                'Accession': accession,
                'Description': description,
                'Sequence': sequence,
                'Length': length
            })
        
        self.ncbi_df = pd.DataFrame(data)
        return self.ncbi_df if return_df else None
    

    def run_prodigal(self, output_dir: Path = Files.PRODIGAL, num_cpus: int = 10,
                     seqs_x_split: int=1000) -> Path:
        """
        Runs Prodigal on the input FASTA file if not already executed.

        Parameters:
        - output_dir (Path): directory where Prodigal outputs will be stored.

        Returns:
        - Path to the .faa file with predicted proteins.
        """
        os.makedirs(output_dir, exist_ok=True)
        
        outdir = Path(output_dir)

        # Busca archivos que terminen en .faa
        faa_files = list(outdir.glob("*.faa"))

        # Si hay al menos uno, coge el primero
        faa_out = faa_files[0] if faa_files else None
        
        # gff_out = output_dir / Files.GFF

        if faa_out: #and gff_out.exists():
            print(f"[INFO] Prodigal output already found at {output_dir}")
            self.faa = faa_out
        else:
            
            faa_out = Path(output_dir) / "output.faa"

            
            script_path = SCRIPTS_DIR / "run_prodigalgv_parallel.sh"
            cmd2 = [
                    str(script_path),
                    str(self.fna_path),
                    str(Path(output_dir)),
                    str(num_cpus),
                    str(seqs_x_split),
                    "--keep-intermediate"
                ]
            
            cmd = [
                "prodigal-gv",
                "-i", str(self.fna_path),
                "-a", str(faa_out),
                "-p", "meta",  # Metagenomic mode
                "-q"           # Quiet mode (suppress output)
            ]

            subprocess.run(cmd2, check=True, stdout=subprocess.DEVNULL)

            try:
                try:
                    print("[INFO] Running (parallel) Prodigal in metagenomic mode...")
                    subprocess.run(cmd2, check=True, stdout=subprocess.DEVNULL)
                    print(f"[INFO] Parallel prodigal-gv finished. File(s) saved to: {output_dir}")
                except:
                    print("[INFO] Running Prodigal in metagenomic mode...")
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
                    print(f"[INFO] Prodigal finished. File(s) saved to: {output_dir}")
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"[ERROR] Prodigal execution failed: {e}")

        self.faa = faa_out
        return self.faa




    def run_hmmer(self, hmm_models: Path = Files.HMM_MODELS, 
                  output_dir: Path = Files.HMM_OUTPUT_MULTIPLE, 
                  num_cpus: int=10):
        """
        Calls the external shell script to run hmmsearch in parallel.
        """
        if self.faa is None:
            print(self.faa)
            raise RuntimeError("Protein file (.faa) is missing.")
         
        script_path = SCRIPTS_DIR / "run_hmm_parallel2.sh"
        # print(f'Busco en {script_path}')
        # print(f'ROOT: {DATA_DIR}')
        if not script_path.exists():
            raise FileNotFoundError(f"HMMER execution script not found at {script_path}")
         
        cmd = [str(script_path), str(self.faa), str(hmm_models), str(output_dir), str(num_cpus)]
         
        print(f"[INFO] Launching parallel hmmsearch using {num_cpus} CPUs...")
        subprocess.run(cmd, check=True)
        print("[INFO] HMMER execution completed.")