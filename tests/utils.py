import hashlib
import logging
import subprocess
from typing import List, Dict

import pandas as pd


def get_filtered_dataframe(file: str, cols_to_compare: Dict[str, type]) -> pd.DataFrame:
    return pd.read_csv(
        file,
        sep="\t",
        dtype=cols_to_compare,
        usecols=cols_to_compare.keys(),
    )


def hash_large_file(
    file: str, identity_cols: List[str], score_cols: List[str], round_decimals: int = 6
) -> (int, str):
    """Row count + MD5 of identity columns and rounded score columns.

    Used for files too large to git-commit as full expected-output fixtures.
    Scores are rounded before hashing so cross-node floating-point noise
    (differs at ~15 significant digits between BLAS/hardware) doesn't produce
    false mismatches, while still catching real regressions.
    """
    df = pd.read_csv(file, sep="\t", usecols=identity_cols + score_cols)
    for col in score_cols:
        df[col] = df[col].round(round_decimals)
    df = df.sort_values(by=identity_cols).reset_index(drop=True)
    row_count = len(df)
    digest = hashlib.md5(df.to_csv(index=False).encode()).hexdigest()
    return row_count, digest


def run_cmd(cmd: str, raise_ex: bool = True) -> bool:
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error: {e}")
        if raise_ex:
            raise
        return False
    return True


def get_biosample_names(biosamples_tsv: str) -> List[str]:
    biosample_data = pd.read_csv(biosamples_tsv, sep="\t")
    return biosample_data.iloc[:, 0].tolist()
