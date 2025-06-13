import pandas as pd
import scripts.globals_vars
import os


def ontologies_completion(spectrum_list, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    The `ontologies_completion` function enriches a DataFrame containing spectral data with ontological information.

    Parameters:
    - spectrum_list (pd.DataFrame): The DataFrame containing spectral data with 'INCHIKEY' and initial ontology columns.
    - progress_callback (callable, optional): Function for updating progress during the process.
    - total_items_callback (callable, optional): Function for setting the total number of items to process.
    - prefix_callback (callable, optional): Function to describe the task being performed.
    - item_type_callback (callable, optional): Function to specify the type of items being processed.

    Returns:
    - pd.DataFrame: The enriched DataFrame with the completed ontology information.
    """
    # Add columns with default values 'UNKNOWN'
    spectrum_list['CLASSYFIRE_SUPERCLASS'] = "UNKNOWN"
    spectrum_list['CLASSYFIRE_CLASS'] = "UNKNOWN"
    spectrum_list['CLASSYFIRE_SUBCLASS'] = "UNKNOWN"
    spectrum_list['NPCLASS_PATHWAY'] = "UNKNOWN"
    spectrum_list['NPCLASS_SUPERCLASS'] = "UNKNOWN"
    spectrum_list['NPCLASS_CLASS'] = "UNKNOWN"

    # Count the number of unique INCHIKEYs in the spectrum_list
    num_keys = spectrum_list['INCHIKEY'].nunique()

    # Initialize tracking with the callbacks
    if prefix_callback:
        prefix_callback("updating ontologies:")

    if item_type_callback:
        item_type_callback("rows")

    if total_items_callback:
        total_items_callback(num_keys, 0)  # Set the total at the start

    # Merge spectrum_list with ontologies_df on 'INCHIKEY'
    completed_df = pd.merge(
        spectrum_list,
        scripts.globals_vars.ontologies_df[
            ["INCHIKEY", "CLASSYFIRE_SUPERCLASS", "CLASSYFIRE_CLASS", "CLASSYFIRE_SUBCLASS", "NPCLASS_PATHWAY",
             "NPCLASS_SUPERCLASS", "NPCLASS_CLASS"]
        ],
        on='INCHIKEY',
        how='left'
    )

    # Simulate updating key-by-key for progress tracking
    processed_keys = 0
    for _ in range(num_keys):
        processed_keys += 1
        if progress_callback:
            progress_callback(processed_keys)

    # Replace initial values with the merged ones
    completed_df['CLASSYFIRE_SUPERCLASS'] = completed_df['CLASSYFIRE_SUPERCLASS_y'].combine_first(
        completed_df['CLASSYFIRE_SUPERCLASS_x'])
    completed_df['CLASSYFIRE_CLASS'] = completed_df['CLASSYFIRE_CLASS_y'].combine_first(
        completed_df['CLASSYFIRE_CLASS_x'])
    completed_df['CLASSYFIRE_SUBCLASS'] = completed_df['CLASSYFIRE_SUBCLASS_y'].combine_first(
        completed_df['CLASSYFIRE_SUBCLASS_x'])
    completed_df['NPCLASS_PATHWAY'] = completed_df['NPCLASS_PATHWAY_y'].combine_first(completed_df['NPCLASS_PATHWAY_x'])
    completed_df['NPCLASS_SUPERCLASS'] = completed_df['NPCLASS_SUPERCLASS_y'].combine_first(
        completed_df['NPCLASS_SUPERCLASS_x'])
    completed_df['NPCLASS_CLASS'] = completed_df['NPCLASS_CLASS_y'].combine_first(completed_df['NPCLASS_CLASS_x'])

    # Remove temporary columns
    completed_df.drop(
        columns=[
            'CLASSYFIRE_SUPERCLASS_x', 'CLASSYFIRE_SUPERCLASS_y',
            'CLASSYFIRE_CLASS_x', 'CLASSYFIRE_CLASS_y',
            'CLASSYFIRE_SUBCLASS_x', 'CLASSYFIRE_SUBCLASS_y',
            'NPCLASS_PATHWAY_x', 'NPCLASS_PATHWAY_y',
            'NPCLASS_SUPERCLASS_x', 'NPCLASS_SUPERCLASS_y',
            'NPCLASS_CLASS_x', 'NPCLASS_CLASS_y'
        ],
        inplace=True
    )

    return completed_df