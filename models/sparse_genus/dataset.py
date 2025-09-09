
from torch.utils.data import Dataset
import pandas as pd
from typing import Dict, Any
from sparse_genus.utils import csr_to_sparse_tensor  


class GenusDataset(Dataset):
    """
    Dataset for virus-level VPF vectors and genus labels.
    Expects a DataFrame with:
        - 'hmms_conteos': sparse matrix rows (per virus)
        - 'Genus': true genus label (string)
        - 'Accession': virus ID
    """
    def __init__(self, df: pd.DataFrame, genus_to_idx: Dict[str, int]):
        self.df = df
        self.genus_to_idx = genus_to_idx

        required_cols = {'hmms_conteos', 'Genus', 'Accession'}
        if not required_cols.issubset(df.columns):
            raise ValueError(f"Missing required columns in input dataframe: {required_cols - set(df.columns)}")

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, index: int) -> Any:
        row = self.df.iloc[index]

        x = csr_to_sparse_tensor(row['hmms_conteos'])  # scipy csr → torch sparse
        genus = row['Genus']
        accession = row['Accession']

        y_genus = self.genus_to_idx.get(genus, -1) if pd.notna(genus) else -1
        return x, y_genus, accession
