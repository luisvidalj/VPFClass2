# sparse_family/utils.py
import os
from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split
import torch
import json
import numpy as np
from scipy.sparse import csr_matrix
from typing import List, Tuple
from vpf_classifier.utils.config import ROOT_DIR


def csr_to_sparse_tensor(csr_mat: csr_matrix) -> torch.sparse.FloatTensor:
    """
    Converts a SciPy CSR sparse matrix to a PyTorch sparse COO tensor.

    Args:
        csr_mat (csr_matrix): Input sparse matrix

    Returns:
        torch.sparse.FloatTensor: COO-format sparse tensor
    """
    coo = csr_mat.tocoo()
    indices = torch.tensor(np.vstack((coo.row, coo.col)), dtype=torch.long)
    values = torch.tensor(coo.data, dtype=torch.float32)
    shape = coo.shape
    return torch.sparse_coo_tensor(indices, values, shape)


def collate_fn_family(batch: List[Tuple[torch.Tensor, int, str]]) -> Tuple[torch.Tensor, torch.Tensor, List[str]]:
    """
    Custom collate function to batch sparse tensors for FAMILY classification.

    Args:
        batch: List of (x_sparse, y_family, accession) tuples

    Returns:
        - x_batch: Batched sparse input tensor (torch.sparse_coo_tensor)
        - y_family_batch: Tensor of family indices (LongTensor)
        - accessions: List of virus accessions
    """
    x_sparse_list, y_family_list, accessions = zip(*batch)
    y_family_batch = torch.LongTensor(y_family_list)

    indices_list = []
    values_list = []

    input_size = x_sparse_list[0].size(1)

    for batch_idx, x_sparse in enumerate(x_sparse_list):
        x_sparse = x_sparse.coalesce()
        indices = x_sparse.indices()  # shape (2, nnz); [row, col]
        values = x_sparse.values()

        # Prepend batch index as the new "row"
        new_indices = torch.vstack([
            torch.full((1, indices.shape[1]), batch_idx, dtype=torch.long),
            indices[1:]  # keep column indices
        ])

        indices_list.append(new_indices)
        values_list.append(values)

    indices_batch = torch.cat(indices_list, dim=1)
    values_batch = torch.cat(values_list)

    x_batch = torch.sparse_coo_tensor(
        indices=indices_batch,
        values=values_batch,
        size=(len(batch), input_size),
        dtype=torch.float32
    ).coalesce()

    return x_batch, y_family_batch, list(accessions)


def get_or_create_split(df: pd.DataFrame,
                        msl_tag: str,
                        test_size: float = 0.2,
                        stratify_col: str = "Family",
                        one_sample: bool = True):
    """
    Creates or loads a train/test split based on viral accessions for a given MSL tag.

    Args:
        df (pd.DataFrame): DataFrame with an 'Accession' column and a stratification column.
        msl_tag (str): Release tag (e.g., "MSL35")
        test_size (float): Proportion of the dataset to assign to the test split.
        stratify_col (str): Column name used for stratified splitting (default: Family).
        one_sample (bool): If True, añade las familias con un único ejemplar al train.

    Returns:
        train_df (pd.DataFrame): Subset of df for training
        test_df (pd.DataFrame): Subset of df for testing
    """
    base_dir = ROOT_DIR / Path("tests/sparse_family") / msl_tag
    base_dir.mkdir(parents=True, exist_ok=True)

    train_file = base_dir / "train_accessions.txt"
    test_file = base_dir / "test_accessions.txt"

    config_file = base_dir / "split_config.json"
    config_file.write_text(json.dumps({"one_sample": one_sample, "test_size": test_size,
                                       "stratify_col": stratify_col}, indent=2))

    if train_file.exists() and test_file.exists():
        print(f"[INFO] Using an existing train/test split for {msl_tag} (Family)")
        train_acc = train_file.read_text().splitlines()
        test_acc = test_file.read_text().splitlines()
    else:
        print(f"[INFO] Creating a new train/test split for {msl_tag} (Family)")

        n_missing = df[stratify_col].isna().sum()
        if n_missing:
            print(f"[WARNING] {n_missing} sequences have {stratify_col} missing values. Removing them...")
        df_valid = df[df[stratify_col].notna()].copy()

        fam_counts = df_valid[stratify_col].value_counts()
        valid_taxa = fam_counts[fam_counts >= 2].index
        print(f"[INFO] Subsetting {stratify_col} represented by ≥2 samples...")
        print(f"{stratify_col} in {msl_tag}: {df_valid[stratify_col].nunique()}")
        print(f"{stratify_col} with ≥2 samples in {msl_tag}: {len(valid_taxa)}")

        filtered_data = df_valid[df_valid[stratify_col].isin(valid_taxa)].copy()
        complementary = df_valid.loc[~df_valid[stratify_col].isin(valid_taxa)].copy()

        train_acc, test_acc = train_test_split(
            filtered_data["Accession"],
            test_size=test_size,
            stratify=filtered_data[stratify_col],
            random_state=2000
        )

        if one_sample and not complementary.empty:
            train_acc = pd.concat([train_acc, complementary["Accession"]], axis=0, ignore_index=True)

        train_file.write_text("\n".join(map(str, train_acc)))
        test_file.write_text("\n".join(map(str, test_acc)))

    train_df = df[df["Accession"].isin(train_acc)].copy()
    test_df = df[df["Accession"].isin(test_acc)].copy()

    print(f"[INFO]: \n {len(train_df)} training sequences \n {len(test_df)} test sequences")
    return train_df, test_df

#################################################################################################################
#################################################################################################################
                                     # Fuera del proyecto
#################################################################################################################
#################################################################################################################

from Bio import SeqIO
from pathlib import Path

def extract_fasta_subset(all_fasta_path: str, accession_file: str, output_fasta: str):
    """
    Extrae un subconjunto de secuencias FASTA dado un archivo con accessions completos (como AB000403.1).
    """
    # Cargar accessions
    accessions = set(Path(accession_file).read_text().splitlines())

    # Iterar por el FASTA y seleccionar los que coincidan (con el id completo)
    selected_records = [
        record for record in SeqIO.parse(all_fasta_path, "fasta")
        if record.id in accessions
    ]

    print(f"[INFO] Selected {len(selected_records)} records from {all_fasta_path}")
    SeqIO.write(selected_records, output_fasta, "fasta")
    print(f"[INFO] Saved to {output_fasta}")

