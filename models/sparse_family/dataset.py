from torch.utils.data import Dataset
import pandas as pd
from typing import Dict, Any
from sparse_family.utils import csr_to_sparse_tensor


class FamilyDataset(Dataset):
    """
    Dataset for virus-level VPF vectors and family labels.
    Expects a DataFrame with:
        - 'hmms_conteos': sparse matrix rows (per virus)
        - 'Family': true family label (string)
        - 'Accession': virus ID
    """
    def __init__(self, df: pd.DataFrame, family_to_idx: Dict[str, int]):
        self.df = df
        self.family_to_idx = family_to_idx

        required_cols = {'hmms_conteos', 'Family', 'Accession'}
        if not required_cols.issubset(df.columns):
            missing = required_cols - set(df.columns)
            raise ValueError(f"Missing required columns in input dataframe: {missing}")

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, index: int) -> Any:
        row = self.df.iloc[index]

        # scipy csr → torch sparse
        x = csr_to_sparse_tensor(row['hmms_conteos'])
        family = row['Family']
        accession = row['Accession']

        y_family = self.family_to_idx.get(family, -1) if pd.notna(family) else -1
        return x, y_family, accession
