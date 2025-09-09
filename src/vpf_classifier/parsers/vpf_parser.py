from Bio import SearchIO
import pandas as pd
import numpy as np
from collections import Counter
from scipy.sparse import csr_matrix, vstack, lil_matrix
from sklearn.metrics.pairwise import cosine_similarity
from typing import Optional
import glob
import os
import json

from vpf_classifier.utils.config import Files
from vpf_classifier.utils.config import Constants
from vpf_classifier.parsers.fasta_parser import FastaParser

class VPF_parser:
    def __init__(self, parser: FastaParser, hmm_file=Files.HMM_MODELS, e_value_threshold: Optional[float] = Constants.e_value_threshold, num_cpus = 20):
        self.parser = parser
        self.evalue = e_value_threshold
        self.df_hmm = None
        self.df_virus_hmm = None

        # Check for existing HMMER .tbl files
        output_tbls = glob.glob(str(Files.HMM_OUTPUT_MULTIPLE / "*.tbl"))
        if output_tbls:
            print(f"[INFO] Existing HMMER output found: {len(output_tbls)} .tbl files")
        else:
            print("[INFO] No HMMER output found. Running hmmsearch assuming 12 CPUs for parallelization")
            self.parser.run_hmmer(hmm_models=hmm_file, output_dir=Files.HMM_OUTPUT_MULTIPLE,num_cpus=num_cpus)




    def parse_multiple_hmm(self, unique_hit=False, hmm_output_folder: Optional[str] = Files.HMM_OUTPUT_MULTIPLE):
        attribs = ['id', 'bitscore', 'cluster_num', 'evalue']
        hits = {key: [] for key in ['hmm_name'] + attribs} 

        print(f"[INFO] Parsing HMMER .tbl files from {Files.HMM_OUTPUT_MULTIPLE}...")
        output_files = glob.glob(os.path.join(hmm_output_folder, "*.tbl"))

        for file in output_files:
            with open(file) as handle:
                for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
                    for hit in queryresult.hits:
                        hits['hmm_name'].append(queryresult.id)
                        for attrib in attribs:
                            hits[attrib].append(getattr(hit,attrib))

        self.df_hmm = pd.DataFrame(hits)

        if self.evalue:
            original = self.df_hmm.shape[0]
            self.df_hmm = self.df_hmm[self.df_hmm['evalue'] < self.evalue]
            print(f"[INFO] Filtered HMM hits by e-value < {self.evalue}: {self.df_hmm.shape[0]} / {original} kept")


        self._aggregate_by_virus()
        self._merge_taxonomy()
        # self._add_vpf_counts_sparse_optimized()
        self._add_vpf_counts_sparse_fixed()
        print("[INFO] Aggregated VPF hit matrix available at .df_virus_hmm")


    def _aggregate_by_virus(self):
        aux_df = self.df_hmm.copy()
        aux_df["Accession"] = aux_df["id"].str.rsplit("_", n=1).str[0]
        self.df_virus_hmm = aux_df.groupby("Accession").agg(
            hmms_hits=('hmm_name', list),
            scores=('bitscore', list),
            evalues=('evalue', list),
            protein_accessions=('id', list),
            cluster_nums=('cluster_num', list)
        ).reset_index()


    def _merge_taxonomy(self):
        ncbi_df = self.parser.parse_fasta_to_dataframe(return_df=True)
        self.df_virus_hmm = self.df_virus_hmm.merge(ncbi_df, on="Accession", how="left")

    def _add_vpf_counts_sparse_optimized(self):
        if self.df_virus_hmm is None:
            raise ValueError("HMM parsing must be done first")

        vpf_unicos = sorted({vpf for row in self.df_virus_hmm['hmms_hits'] for vpf in row})
        vpf_to_index = {vpf: i for i, vpf in enumerate(vpf_unicos)}
        num_virus = len(self.df_virus_hmm)
        num_vpfs = len(vpf_unicos)

        sparse_matrix = lil_matrix((num_virus, num_vpfs), dtype=np.int32)

        for i, row in enumerate(self.df_virus_hmm['hmms_hits']):
            conteos = Counter(row)
            for vpf, count in conteos.items():
                j = vpf_to_index[vpf]
                sparse_matrix[i, j] = count

        self.vpf_sparse_matrix = sparse_matrix.tocsr()
        self.vpf_unicos = vpf_unicos
        self.vpf_to_index = vpf_to_index
        self.df_virus_hmm['hmms_conteos'] = [self.vpf_sparse_matrix.getrow(i) for i in range(num_virus)]


    
    def _add_vpf_counts_sparse_fixed(self):
        """
        Uses a precomputed vpf_to_index dictionary to ensure all vectors have the same length.
        Builds a fixed-length sparse matrix of VPF counts.
        """
        dict_path = Files.HMM_DICT
        if not dict_path.exists():
            raise FileNotFoundError(f"[ERROR] Expected dictionary not found at {dict_path}")

        with open(dict_path) as f:
            vpf_to_index = json.load(f)

        self.vpf_to_index = vpf_to_index
        vpf_dim = len(vpf_to_index)
        num_virus = len(self.df_virus_hmm)

        sparse_matrix = lil_matrix((num_virus, vpf_dim), dtype=np.int32)

        missing_vpfs = set()

        for i, row in enumerate(self.df_virus_hmm['hmms_hits']):
            counts = Counter(row)
            for vpf, count in counts.items():
                if vpf in vpf_to_index:
                    j = vpf_to_index[vpf]
                    sparse_matrix[i, j] = count
                else:
                    missing_vpfs.add(vpf)

        if missing_vpfs:
            print(f"[WARNING] {len(missing_vpfs)} VPFs found in hits not present in vpf_to_index")

        self.vpf_sparse_matrix = sparse_matrix.tocsr()
        self.df_virus_hmm['hmms_conteos'] = [self.vpf_sparse_matrix.getrow(i) for i in range(num_virus)]


    def get_sparse_matrix_from_dataframe(self):
        if 'hmms_conteos' not in self.df_virus_hmm.columns:
            raise ValueError('parse_multiple_hmm must be run first')
        sparse_matrix = vstack(self.df_virus_hmm['hmms_conteos'].values)
        if not isinstance(sparse_matrix, csr_matrix):
            sparse_matrix = sparse_matrix.tocsr()
        return sparse_matrix

    def calculate_cosine_similarity(self, normalize=True):
        if not hasattr(self, 'vpf_sparse_matrix'):
            raise ValueError('Sparse matrix not initialized')

        matrix = self.vpf_sparse_matrix
        if normalize:
            norms = np.sqrt(matrix.power(2).sum(axis=1))
            norms[norms == 0] = 1
            matrix = matrix.multiply(1 / norms)

        return cosine_similarity(matrix)