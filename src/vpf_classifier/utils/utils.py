import pandas as pd
import re

def clean_and_split_accessions(raw_string):
    """
    Given an accession string from the ICTV MSL file, splits it into clean GenBank accessions.

    Handles cases like:
    - DNA-A: HM585445; DNA-B: HM585446
    - AE006468 (2844298.2877981)
    - U65738; U65739; U65737

    Returns:
        List of clean accession strings.
    """
    if pd.isna(raw_string):
        return []

    split_items = raw_string.split(';')
    clean_accessions = []

    for item in split_items:
        # Remove segment labels like "DNA-A:", "RNA2:", etc.
        if ':' in item:
            item = item.split(':')[-1]

        # Remove any parenthesis content
        item = re.sub(r'\(.*?\)', '', item)

        # Strip whitespace
        item = item.strip()

        # Remove trailing GenBank version numbers (e.g., ".1")
        if '.' in item:
            item = item.split('.')[0]

        # Validate basic format
        if re.match(r'^[A-Za-z0-9_]+$', item):
            clean_accessions.append(item)

    return clean_accessions
