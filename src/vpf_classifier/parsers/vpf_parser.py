from Bio import SearchIO
import pandas as pd
import numpy as np
from collections import Counter
from scipy.sparse import csr_matrix, vstack, lil_matrix
from sklearn.metrics.pairwise import cosine_similarity
from typing import Optional
import glob, os
from itertools import chain
import json
from pathlib import Path
from typing import Optional
from multiprocessing import get_context

from vpf_classifier.utils.config import Files
from vpf_classifier.utils.config import Constants
from vpf_classifier.parsers.fasta_parser import FastaParser

class VPF_parser:
    # def __init__(self, parser: FastaParser, hmm_file=Files.HMM_MODELS, e_value_threshold: Optional[float] = Constants.e_value_threshold, num_cpus = 20):
    #     self.parser = parser
    #     self.evalue = e_value_threshold
    #     self.df_hmm = None
    #     self.df_virus_hmm = None

    #     # Check for existing HMMER .tbl files
    #     output_tbls = glob.glob(str(Files.HMM_OUTPUT_MULTIPLE / "*.tbl"))
    #     if output_tbls:
    #         print(f"[INFO] Existing HMMER output found: {len(output_tbls)} .tbl files")
    #     else:
    #         print("[INFO] No HMMER output found. Running hmmsearch assuming 12 CPUs for parallelization")
    #         self.parser.run_hmmer(hmm_models=hmm_file, output_dir=Files.HMM_OUTPUT_MULTIPLE,num_cpus=num_cpus)

    def __init__(
        self,
        parser: FastaParser,
        hmm_file = Files.HMM_MODELS,
        e_value_threshold: Optional[float] = Constants.e_value_threshold,
        num_cpus: int = 20,
        # --- NUEVO ---
        vpf_dict_path: Optional[Path] = None,
        hmm_output_dir: Optional[Path] = None,
        vector_norm: Optional[str] = None, # l2, l1 or None
    ):
        self.parser = parser
        self.evalue = e_value_threshold
        self.df_hmm = None
        self.df_virus_hmm = None

        # --- Almacenar overrides ---
        self._vpf_dict_path = vpf_dict_path
        self._hmm_output_dir = Path(hmm_output_dir) if hmm_output_dir is not None else Files.HMM_OUTPUT_MULTIPLE

        # Donam opcio a normalitzar
        self.vector_norm = vector_norm

        # Miramos si ya hay .tbl en la carpeta seleccionada; si no, ejecutamos hmmsearch ahí
        output_tbls = glob.glob(str(self._hmm_output_dir / "*.tbl"))
        if not output_tbls:
            os.makedirs(self._hmm_output_dir, exist_ok=True)
            print(f"[INFO] No HMMER output found. Running hmmsearch -> {self._hmm_output_dir} (cpus={num_cpus})")
            self.parser.run_hmmer(
                hmm_models=hmm_file,
                output_dir=self._hmm_output_dir,
                num_cpus=num_cpus
            )
        else:
            print(f"[INFO] HMMER ouput already found.")


    # def parse_multiple_hmm(self, unique_hit=False, hmm_output_folder: Optional[str] = Files.HMM_OUTPUT_MULTIPLE):
    #     attribs = ['id', 'bitscore', 'cluster_num', 'evalue']
    #     hits = {key: [] for key in ['hmm_name'] + attribs} 

    #     print(f"[INFO] Parsing HMMER .tbl files from {Files.HMM_OUTPUT_MULTIPLE}...")
    #     output_files = glob.glob(os.path.join(hmm_output_folder, "*.tbl"))
    def parse_multiple_hmm(self, unique_hit: bool = False, hmm_output_folder: Optional[str] = None):
        attribs = ['id', 'bitscore', 'cluster_num', 'evalue']
        hits = {key: [] for key in ['hmm_name'] + attribs}

        folder = Path(hmm_output_folder) if hmm_output_folder else self._hmm_output_dir
        print(f"[INFO] Parsing HMMER .tbl files from {folder}...")
        output_files = glob.glob(str(folder / "*.tbl"))

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
        # print("[INFO] Aggregated VPF hit matrix available at .df_virus_hmm")


    def _parse_tbl_file_worker(args):
        file, evalue_threshold = args
        # Salida: tuplas con el mismo esquema que usas en parse_multiple_hmm
        # columnas --> ['hmm_name','id','bitscore','cluster_num','evalue']
        rows = []
        total_rows = 0
        with open(file) as handle:
            for qres in SearchIO.parse(handle, 'hmmer3-tab'):
                qid = qres.id  # == hmm_name en tu implementación actual
                for hit in qres.hits:
                    ev = getattr(hit, 'evalue', None)
                    if (evalue_threshold is None) or (ev is not None and ev < evalue_threshold):
                        total_rows +=1 
                        rows.append((
                            qid,
                            getattr(hit, 'id', None),
                            float(getattr(hit, 'bitscore', 0.0)),
                            getattr(hit, 'cluster_num', None),
                            float(ev) if ev is not None else None
                        ))
                    elif (evalue_threshold is not None) and (ev >= evalue_threshold):
                        total_rows += 1
        return rows, total_rows
    

    def parse_multiple_hmm_parallel(
    self,
    unique_hit: bool = False,                   # reservado por si más adelante
    hmm_output_folder: Optional[str] = None,
    num_cpus: Optional[int] = None,
    ):
        """
        Versión paralela de parse_multiple_hmm que genera EXACTAMENTE el mismo
        DataFrame (columnas y tipos), salvo por el orden de las filas.
        """

        # Misma definición de columnas que en tu versión secuencial
        columns = ['hmm_name', 'id', 'bitscore', 'cluster_num', 'evalue']

        folder = Path(hmm_output_folder) if hmm_output_folder else self._hmm_output_dir
        print(f"[INFO] (Parallel) Parsing HMMER .tbl files from {folder}...")
        output_files = sorted(glob.glob(str(folder / "*.tbl")))  # ordenado por reproducibilidad

        if not output_files:
            print("[WARN] No .tbl files found; resulting df_hmm will be empty.")
            self.df_hmm = pd.DataFrame(columns=columns)
            # Mantiene el mismo flujo aguas abajo
            self._aggregate_by_virus()
            self._merge_taxonomy()
            self._add_vpf_counts_sparse_fixed()
            return

        # CPUs
        if num_cpus is None or num_cpus < 1:
            num_cpus = max(1, (os.cpu_count() or 1) - 1)

        # Usamos 'fork' si existe (Linux/macOS) para evitar problemas con __main__;
        # si no existe (p.ej. Windows), get_context lanzará y caemos a 'spawn'.
        try:
            ctx = get_context("fork")
        except ValueError:
            ctx = get_context()  # fallback (spawn/forkserver según plataforma)

        # Map paralelo por archivo
        with ctx.Pool(processes=num_cpus) as pool:
            results = pool.map(
                VPF_parser._parse_tbl_file_worker,
                [(f, self.evalue) for f in output_files]
            )

        # Flatten de todas las filas
        row_list, total_counts = zip(*results)
        all_rows = list(chain.from_iterable(row_list))
        total_before = sum(total_counts)

        # Construimos el DataFrame con el mismo esquema que el método original
        self.df_hmm = pd.DataFrame.from_records(all_rows, columns=columns)

        # Nota: ya hemos filtrado por evalue en el worker. Si quieres mantener el
        # log informativo y asegurar igualdad semántica, re-aplicamos el filtro:
        if self.evalue is not None and len(self.df_hmm) > 0:
            self.df_hmm = self.df_hmm[self.df_hmm['evalue'] < self.evalue]
            print(f"[INFO] Filtered HMM hits by e-value < {self.evalue}: {self.df_hmm.shape[0]} / {total_before} kept")

        # Resto del pipeline idéntico
        self._aggregate_by_virus()
        self._merge_taxonomy()
        self._add_vpf_counts_sparse_fixed()




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
        # dict_path = Files.HMM_DICT
        dict_path = Path(self._vpf_dict_path) if self._vpf_dict_path else Files.HMM_DICT
        if not dict_path.exists():
            raise FileNotFoundError(f"[ERROR] Expected dictionary not found at {dict_path}")

        with open(dict_path) as f:
            vpf_to_index = json.load(f)

        self.vpf_to_index = vpf_to_index
        vpf_dim = len(vpf_to_index)
        num_virus = len(self.df_virus_hmm)

        sparse_matrix = lil_matrix((num_virus, vpf_dim), dtype=np.float32)

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
        self.vpf_sparse_matrix_counts = self.vpf_sparse_matrix.copy()

        if self.vector_norm in ("l1","l2"):
            print("NORMALITZAM")
            self.vpf_sparse_matrix = self._normalize_sparse_rows(self.vpf_sparse_matrix,
                                                                        mode = self.vector_norm)
        
        self.df_virus_hmm['hmms_conteos'] = [self.vpf_sparse_matrix.getrow(i) for i in range(num_virus)]


    def get_sparse_matrix_from_dataframe(self):
        if 'hmms_conteos' not in self.df_virus_hmm.columns:
            raise ValueError('parse_multiple_hmm must be run first')
        sparse_matrix = vstack(self.df_virus_hmm['hmms_conteos'].values)
        if not isinstance(sparse_matrix, csr_matrix):
            sparse_matrix = sparse_matrix.tocsr()
        return sparse_matrix
    

    def _normalize_sparse_rows(self, M: csr_matrix, mode: str = "l2") -> csr_matrix:
        """Normaliza cada fila de M por su norma (L2 o L1). Las filas todo-cero quedan iguales."""
        if not isinstance(M, csr_matrix):
            M = M.tocsr()
        n = M.shape[0]
        if mode == "l2":
            # norma L2 por fila
            row_norms = np.sqrt(M.power(2).sum(axis=1)).A.ravel()
        elif mode == "l1":
            # norma L1 por fila
            row_norms = np.asarray(M.sum(axis=1)).ravel()
        else:
            raise ValueError(f"Unknown norm mode: {mode}")
        
        print(row_norms)

        # Evita división por 0: filas todo-cero se dejan tal cual
        row_norms[row_norms == 0] = 1.0
        inv = 1.0 / row_norms
        return M.multiply(inv[:, None]).tocsr()


    def calculate_cosine_similarity(self, normalize=True):
        if not hasattr(self, 'vpf_sparse_matrix'):
            raise ValueError('Sparse matrix not initialized')

        matrix = self.vpf_sparse_matrix
        if normalize:
            norms = np.sqrt(matrix.power(2).sum(axis=1))
            norms[norms == 0] = 1
            matrix = matrix.multiply(1 / norms)

        return cosine_similarity(matrix)