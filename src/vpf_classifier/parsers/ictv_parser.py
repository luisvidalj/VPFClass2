import pandas as pd
from vpf_classifier.utils.utils import clean_and_split_accessions


def clean_ictv_csv(ictv_df: pd.DataFrame, msl_tag: str) -> pd.DataFrame:
    
    # 1. Genome Coverage --> Genome (for some releases)
    if 'Genome Coverage' in ictv_df.columns:
        ictv_df = ictv_df.rename(columns={"Genome Coverage": "Genome"})

    if 'Virus GENBANK accession' not in ictv_df.columns:
        raise ValueError("Virus GENBANK accession not found")
        
    # 2. Creating a 'base accession' columns to merge the MSLXX.fasta
    expanded_rows = []

    for _, row in ictv_df.iterrows():
        accessions = clean_and_split_accessions(row['Virus GENBANK accession'])
        for i in range(len(accessions)):
            if '.' in accessions[i]:
                new = accessions[i].split('.')[0]
                accessions[i] = new
                print(new)
        for acc in accessions:
            new_row = row.copy()
            new_row['base_accession'] = acc
            expanded_rows.append(new_row)


    ictv_df_expanded = pd.DataFrame(expanded_rows)

    return ictv_df_expanded


def merge_vpf_with_ictv(vpf_database: pd.DataFrame, ictv_df: pd.DataFrame, out_of_use: bool = False) -> pd.DataFrame:
    """
    Merge VPF dataframe with cleaned ICTV annotations using base accessions.

    Args:
        vpf_database (pd.DataFrame): DataFrame from VPF_parser, must contain 'Accession' column.
        ictv_df (pd.DataFrame): Raw ICTV dataframe loaded from .csv.
        out_of_use: (internal use) -> set true to use with generic df

    Returns:
        pd.DataFrame: Merged DataFrame ready for training, with taxonomic and genome metadata.
    """
    # 1. Normalize accession in VPF to base format (remove .X if present)
    vpf_database = vpf_database.copy()
    vpf_database['base_accession'] = vpf_database['Accession'].str.split('.').str[0]

    # 2. Expand the ICTV DataFrame (duplicate rows for each accession if needed)
    expanded_rows = []
    for _, row in ictv_df.iterrows():
        accessions = clean_and_split_accessions(row['Virus GENBANK accession'])
        for i in range(len(accessions)):
            if '.' in accessions[i]:
                accessions[i] = accessions[i].split('.')[0]
        for acc in accessions:
            new_row = row.copy()
            new_row['base_accession'] = acc
            expanded_rows.append(new_row)
    df_ictv_expanded = pd.DataFrame(expanded_rows)

    # 3. Restrict columns of ICTV to relevant ones
    taxonomy_cols = [
        'Isolate ID','Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum',
        'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily',
        'Genus', 'Subgenus', 'Species'
    ]
    other_cols = ['base_accession', 'Genome', 'Genome composition', 'Host source']
    cols_to_keep = taxonomy_cols + other_cols
    df_ictv_expanded = df_ictv_expanded[[col for col in cols_to_keep if col in df_ictv_expanded.columns]]
    # vpf_database = vpf_database[[col for col in cols_to_keep if col in vpf_database.columns]]

    # 4. Merge on base_accession
    df_merged = pd.merge(
        vpf_database,
        df_ictv_expanded,
        on='base_accession',
        how='inner',
        suffixes=('', '_drop')
    )

    # Check duplicates
    duplicated = df_merged[df_merged.duplicated(subset='Accession', keep=False)]

    if not duplicated.empty:
        inconsistent = duplicated.groupby('Accession')['Genus'].nunique()
        problematic = inconsistent[inconsistent > 1]

        if not problematic.empty:
            # raise ValueError(
            #     f"Found duplicated accessions with inconsistent Genus labels:\n{problematic.index.tolist()}"
            # )
            print(f"[IMPORTANT] The following Accession are problematic due to MSL release conflict: \n{problematic.index.tolist()}")
        else:
            print(f"[INFO] Duplicated accessions found but with consistent Genus labels: {duplicated['Accession'].nunique()} entries")


    # 5. Add count of protein hits (based on unique proteins)
    if out_of_use:
        pass
    else:
        df_merged['prots_hits'] = df_merged['protein_accessions'].apply(lambda x: len(set(x)))

    # 6. Remove duplicate virus entries (same Accession)
    print("[WARNING] Duplicated accessions after merge:", df_merged[df_merged.duplicated(subset='Accession', keep=False)]['Accession'].unique().tolist())
    df_merged = df_merged.drop_duplicates(subset='Accession', keep='first').reset_index(drop=True)

    # 7. Fill missing Family and Subfamily with 'Unknown'
    df_merged['Family'] = df_merged['Family'].fillna('Unknown')
    df_merged['Subfamily'] = df_merged['Subfamily'].fillna('Unknown')

    return df_merged
