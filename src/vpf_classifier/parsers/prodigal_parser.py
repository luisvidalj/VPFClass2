from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from vpf_classifier.parsers.fasta_parser import FastaParser


from pathlib import Path
from vpf_classifier.utils.config import Files
from pathlib import Path
from typing import Optional
import glob

class Prodigal:
    # def __init__(self, parser: FastaParser):
    #     """
    #     Wraps a FastaParser instance and processes predicted protein sequences from Prodigal output.
    #     Will run Prodigal only if no .faa file is found in the expected directory.
    #     """
    #     self.parser = parser
    #     faa_candidates = sorted(glob.glob(str(Files.PRODIGAL / "*.faa")))

    #     if faa_candidates:
    #         self.parser.run_prodigal()
    #         self.output_faa = Path(faa_candidates[0])
    #         if len(faa_candidates) > 1:
    #             print(f"[WARNING] Multiple .faa files found in {Files.PRODIGAL}. Using: {self.output_faa.name}")
    #         else:
    #             print(f"[INFO] Using existing Prodigal .faa: {self.output_faa.name}")
            
    #     else:
    #         print("[INFO] No .faa file found. Calling Prodigal...")
    #         self.output_faa = self.parser.run_prodigal()

    #     self.ncbi = self.parser.parse_fasta_to_dataframe(return_df=True)
    #     self.df_prots = None
    #     self.df_virus_prot = None

    def __init__(self, parser: FastaParser, output_dir: Optional[Path] = None):
        """
        Wraps a FastaParser instance and processes predicted protein sequences from Prodigal output.
        Delegates to FastaParser.run_prodigal(), which is idempotent:
        - If .faa/.gff already exist in output_dir, reuses them.
        - Otherwise, runs Prodigal and returns the new paths.
        """
        self.parser = parser
        # Si no te pasan carpeta, usa la de Files (modo entrenamiento).
        target_dir = output_dir if output_dir is not None else Files.PRODIGAL
        print(target_dir)
        self.output_faa = self.parser.run_prodigal(output_dir=target_dir)

        # FASTA metadata (contigs) para poder hacer merges más tarde
        self.ncbi = self.parser.parse_fasta_to_dataframe(return_df=True)

        # Placeholders
        self.df_prots = None
        self.df_virus_prot = None
    
    def parse_prodigal(self):
        """
        Parses the .faa file produced by Prodigal and returns a DataFrame aggregated by virus (accession).
        """
        data = {}
        for record in SeqIO.parse(self.output_faa, "fasta"):
            header_parts = record.description.split('#')
            seq_id = header_parts[0].strip()
            protein_seq = str(record.seq)

            if seq_id not in data:
                data[seq_id] = []
            data[seq_id].append(protein_seq)

        self.df_prots = pd.DataFrame([
            {
                "Protein_accession": seq_id,
                "Genes": len(proteins),
                "Proteins": str(*proteins),  # one protein per protein_accession assumption
                "Length": len(str(proteins))
            }
            for seq_id, proteins in data.items()
        ])

        self._aggregate_by_virus()
        self._merge_taxonomy()

        return self.df_virus_prot
    
    def _aggregate_by_virus(self):
        """
        Aggregates the protein-level DataFrame into a virus-level DataFrame based on accession.
        """
        aux_df = self.df_prots.copy()
        aux_df["Accession"] = aux_df["Protein_accession"].str.rsplit("_", n=1).str[0]

        self.df_virus_prot = aux_df.groupby("Accession").agg(
            protein_sequences=("Proteins", list),
            protein_accessions=("Protein_accession", list),
            protein_lengths=("Proteins", lambda x: [len(seq) for seq in x])
        ).reset_index()

    def _merge_taxonomy(self):
        """
        Merges taxonomic metadata from the input FASTA file into the virus-level protein DataFrame.
        """
        self.df_virus_prot = self.ncbi.merge(self.df_virus_prot, on="Accession", how="left")


    # # METHODS FOR DOWNSTREAM ANALYSIS (not tested)

    # def distr_number_proteins(self, taxonomy: str = None):
    #     """
    #     Plots the distribution of the number of proteins predicted per virus.

    #     Args:
    #         taxonomy (str): Optional. Taxonomic level to group by (e.g., 'Family', 'Genus').
    #     """
    #     if self.df_virus_prot is None:
    #         raise ValueError("parse_prodigal() must be executed first.")

    #     plot_df = self.df_virus_prot.copy()

    #     # Count number of proteins per virus
    #     plot_df["num_proteins"] = plot_df["protein_accessions"].apply(lambda x: len(x) if isinstance(x, list) else 0)

    #     if taxonomy:
    #         if taxonomy not in plot_df.columns:
    #             raise ValueError(f"Taxonomic level '{taxonomy}' not found in dataframe.")
    #         plot_df = plot_df.groupby(taxonomy)["num_proteins"].sum().reset_index()

    #     plt.figure(figsize=(10, 5))
    #     sns.histplot(plot_df["num_proteins"], bins=50, kde=True)
    #     plt.xlabel("Number of Proteins per Virus")
    #     plt.ylabel("Frequency")
    #     plt.title("Distribution of Predicted Proteins per Virus (Prodigal)")
    #     plt.grid(True)
    #     plt.tight_layout()
    #     plt.show()


    # def number_of_proteins(self, taxonomy: str = None, cmap: str = "Blues", order: str = "proteins"):
    #     """
    #     Visualizes the frequency distribution of number of proteins predicted per virus using a heatmap.

    #     Args:
    #         taxonomy (str): Optional. Grouping level (e.g., 'Family').
    #         cmap (str): Matplotlib colormap.
    #         order (str): 'proteins' to sort by number of proteins, 'counts' to sort by frequency.
    #     """
    #     if self.df_virus_prot is None:
    #         raise ValueError("parse_prodigal() must be executed first.")

    #     plot_df = self.df_virus_prot.copy()
    #     plot_df["num_proteins"] = plot_df["protein_accessions"].apply(lambda x: len(x) if isinstance(x, list) else 0)

    #     if taxonomy:
    #         if taxonomy not in plot_df.columns:
    #             raise ValueError(f"Taxonomic level '{taxonomy}' not found in dataframe.")
    #         plot_df = plot_df.groupby(taxonomy)["num_proteins"].sum().reset_index()

    #     if order == "proteins":
    #         counts = plot_df["num_proteins"].value_counts().sort_index(ascending=False)
    #     elif order == "counts":
    #         counts = plot_df["num_proteins"].value_counts().sort_values(ascending=False)
    #     else:
    #         raise ValueError("Parameter 'order' must be 'proteins' or 'counts'.")

    #     original_len = len(counts)
    #     counts = counts.head(150) if original_len > 150 else counts

    #     num_categories = len(counts)
    #     num_plots = int(np.ceil(num_categories / 50))
    #     fig, axes = plt.subplots(num_plots, 1, figsize=(10, 5 * num_plots), constrained_layout=True)

    #     if num_plots == 1:
    #         axes = [axes]

    #     vmin, vmax = counts.min(), counts.max()

    #     for i, ax in enumerate(axes):
    #         start, end = i * 50, min((i + 1) * 50, num_categories)
    #         subset = counts.iloc[start:end].to_frame()

    #         sns.heatmap(
    #             subset.T, ax=ax, cmap=cmap,
    #             cbar=(i == num_plots - 1),
    #             annot=True, fmt="d", linewidths=0.5,
    #             vmin=vmin, vmax=vmax,
    #             annot_kws={"rotation": 90}
    #         )
    #         ax.set_xlabel("Number of Proteins per Virus" if not taxonomy else f"per {taxonomy}")
    #         ax.set_title(f"Distribution of Predicted Proteins\nLower Frequencies Hidden: {original_len - len(counts)}")

        



