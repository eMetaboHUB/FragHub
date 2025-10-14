import pandas as pd
import numpy as np

def normalize_to_not_found(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes empty string values in a DataFrame by replacing them with "NOT FOUND".

    This function creates a copy of the input DataFrame, then replaces all empty string values
    within the DataFrame with the string "NOT FOUND". The modified DataFrame is returned as a result.

    Parameters:
    df: pd.DataFrame
        The input DataFrame to process, with potential empty string values.

    Returns:
    pd.DataFrame
        A copy of the input DataFrame with all empty string values replaced by "NOT FOUND".
    """
    df_processed = df.copy()
    df_processed = df_processed.replace('', 'NOT FOUND')

    return df_processed