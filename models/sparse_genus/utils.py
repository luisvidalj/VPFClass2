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



def collate_fn_genus(batch: List[Tuple[torch.Tensor, int, str]]) -> Tuple[torch.Tensor, torch.Tensor, List[str]]:
    """
    Custom collate function to batch sparse tensors for genus classification.

    Args:
        batch: List of (x_sparse, y_genus, accession) tuples

    Returns:
        - x_batch: Batched sparse input tensor (torch.sparse_coo_tensor)
        - y_genus_batch: Tensor of genus indices (LongTensor)
        - accessions: List of virus accessions
    """
    x_sparse_list, y_genus_list, accessions = zip(*batch)
    y_genus_batch = torch.LongTensor(y_genus_list)

    indices_list = []
    values_list = []

    input_size = x_sparse_list[0].size(1)

    for batch_idx, x_sparse in enumerate(x_sparse_list):
        x_sparse = x_sparse.coalesce()
        indices = x_sparse.indices()
        values = x_sparse.values()

        # Prepend batch index to row coordinates
        new_indices = torch.vstack([
            torch.full((1, indices.shape[1]), batch_idx, dtype=torch.long),
            indices[1:]
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

    return x_batch, y_genus_batch, list(accessions)





def get_or_create_split(df: pd.DataFrame, msl_tag: str, test_size: float = 0.2, stratify_col: str = "Genus", one_sample: bool = True):
    """
    Creates or loads a train/test split based on viral accessions for a given MSL tag.

    Args:
        df (pd.DataFrame): DataFrame with an 'Accession' column and a stratification column.
        msl_tag (str): Release tag (e.g., "MSL35")
        test_size (float): Proportion of the dataset to assign to the test split.
        stratify_col (str): Column name used for stratified splitting.

    Returns:
        train_df (pd.DataFrame): Subset of df for training
        test_df (pd.DataFrame): Subset of df for testing
    """

    base_dir = ROOT_DIR / Path("tests/sparse_genus") / msl_tag
    base_dir.mkdir(parents=True, exist_ok=True)

    train_file = base_dir / "train_accessions.txt"
    test_file = base_dir / "test_accessions.txt"

    config_file = base_dir / "split_config.json"
    config_file.write_text(json.dumps({"one_sample": one_sample, "test_size": test_size}, indent=2))


    if train_file.exists() and test_file.exists():
        print(f"[INFO] Using a existing train/test split for {msl_tag}")
        train_acc = train_file.read_text().splitlines()
        test_acc = test_file.read_text().splitlines()
    else:
        print(f"[INFO] Creating a new train/test splits for {msl_tag}")

        print(f"[WARNING] {df[stratify_col].isna().sum()} sequences have Genus missing values. Removing them...")
        df_valid = df[df[stratify_col].notna()].copy() # should not have missing values

        genus_counts = df_valid[stratify_col].value_counts()
        valid_taxa = genus_counts[genus_counts >= 2].index
        print(f"[INFO] Substracting Genus represented by more than one sample...")
        print(f"Genus in {msl_tag}: {len(list(set(df_valid[stratify_col].values)))}")
        print(f"Genus with more than one sample in {msl_tag}: {len(valid_taxa)}")

        filtered_data = df_valid[df_valid[stratify_col].isin(values=valid_taxa)].copy()
        complementary = df_valid.loc[~df_valid[stratify_col].isin(valid_taxa)].copy()

        train_acc, test_acc = train_test_split(
            filtered_data["Accession"],
            test_size=test_size,
            stratify=filtered_data[stratify_col],
            random_state=2000
        )

        if one_sample:
            train_acc = pd.concat([train_acc,complementary["Accession"]], axis=0, ignore_index=True)

        train_file.write_text("\n".join(map(str, train_acc)))
        test_file.write_text("\n".join(map(str, test_acc)))

    train_df = df[df["Accession"].isin(train_acc)].copy()
    test_df = df[df["Accession"].isin(test_acc)].copy()

    print(f"[INFO]: \n {len(train_df)} training sequences \n {len(test_df)} test sequences")

    return train_df, test_df



import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit

def stratified_split(
    df: pd.DataFrame,
    test_size: float = 0.2,
    random_state: int = 42,
    strategy: int = 1,
    taxo_col: str = "Genus"
):
    """
    Crea splits train/test para tres escenarios:
      1) Estratificado SOLO con géneros con >=2 muestras (los de 1 muestra se excluyen del split).
      2) Como (1), pero los géneros con 1 muestra van a TRAIN.
      3) Como (1), pero los géneros con 1 muestra van a TEST (útiles para OSR).

    Notas:
    - Las filas con Genus NaN no se usan para estratificar; se añaden a TRAIN por defecto.
    - Devuelve (train_df, test_df) con índices preservados del df original.

    Parámetros:
      df: DataFrame con columna 'Genus'
      test_size: proporción de test para el estrato de géneros con >=2 muestras
      strategy: 1 | 2 | 3 (ver descripción)
      taxo_col: nombre de la columna de etiqueta

    """
    assert strategy in (1, 2, 3), "strategy debe ser 1, 2 o 3"

    df = df.copy()

    # Separar etiquetados vs no etiquetados
    df_labeled = df[~df[taxo_col].isna()].copy()
    df_unlabeled = df[df[taxo_col].isna()].copy()  # los mandaremos a train

    # Contar por género
    counts = df_labeled[taxo_col].value_counts()
    singletons = counts[counts == 1].index
    multis     = counts[counts >= 2].index

    df_multi = df_labeled[df_labeled[taxo_col].isin(multis)].copy()
    df_single = df_labeled[df_labeled[taxo_col].isin(singletons)].copy()

    # Stratified split en el subconjunto con >=2
    if len(df_multi) > 0:
        sss = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
        idx = np.arange(len(df_multi))
        y = df_multi[taxo_col].values
        train_idx, test_idx = next(sss.split(idx, y))
        multi_train = df_multi.iloc[train_idx]
        multi_test  = df_multi.iloc[test_idx]
    else:
        # Raro, pero por si acaso
        multi_train = df_multi.iloc[:0]
        multi_test  = df_multi.iloc[:0]

    # Reubicar singletons según estrategia
    if strategy == 1:
        # Se excluyen del split (no van ni a train ni a test)
        single_train = df_single.iloc[:0]
        single_test  = df_single.iloc[:0]
    elif strategy == 2:
        # Todos los singletons a TRAIN
        single_train = df_single
        single_test  = df_single.iloc[:0]
    else:  # strategy == 3
        # Todos los singletons a TEST (para evaluar OSR)
        single_train = df_single.iloc[:0]
        single_test  = df_single

    # Añadir no etiquetados a TRAIN (opción conservadora)
    train_df = pd.concat([multi_train, single_train, df_unlabeled], axis=0).sort_index()
    test_df  = pd.concat([multi_test,  single_test], axis=0).sort_index()

    # Info útil
    def summary(name, d):
        total = len(d)
        uniq = d[taxo_col].dropna().nunique()
        singles = (d[taxo_col].value_counts() == 1).sum()
        print(f"{name}: N={total}, géneros distintos (etiquetados)={uniq}, singletons={singles}")
    print("=== Resumen del split ===")
    print(f"Estrategia: {strategy}")
    print(f"Singletons totales: {len(df_single)}")
    summary("TRAIN", train_df)
    summary("TEST ", test_df)

    # Comprobar intersección de géneros entre train/test (solo etiquetados)
    train_g = set(train_df[taxo_col].dropna().unique())
    test_g  = set(test_df[taxo_col].dropna().unique())
    inter   = train_g.intersection(test_g)
    print(f"Géneros en común entre TRAIN y TEST: {len(inter)}")

    if strategy == 3:
        # En 3, esperamos que los singletons estén sólo en TEST (no en TRAIN)
        only_test_singletons = set(singletons) - train_g
        print(f"Singletons ubicados exclusivamente en TEST: {len(only_test_singletons)}")

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


