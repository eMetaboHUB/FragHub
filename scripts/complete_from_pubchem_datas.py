import pandas as pd
import globals_vars
import os


def complete_from_pubchem_datas(CONCATENATE_DF, progress_callback=None, total_items_callback=None, prefix_callback=None,
                                item_type_callback=None):
    """
    Enriches the data from CONCATENATE_DF with data from the PubChem DataFrame (pubchem_datas),
    by matching on the 'INCHIKEY' column and updating specific columns to include missing information.

    :param CONCATENATE_DF: DataFrame to enrich
    :param progress_callback: A function to update progress during processing (optional)
    :param total_items_callback: A function to set the total number of items (optional)
    :param prefix_callback: A function to describe the task in progress (optional)
    :param item_type_callback: A function to describe the type of items being processed (optional)
    :return: Enriched DataFrame
    """
    # Define the prefix for the task and type of items if provided
    if prefix_callback:
        prefix_callback("enriching data from PubChem:")
    if item_type_callback:
        item_type_callback("rows")

    # Create a copy of CONCATENATE_DF to avoid modifying the original DataFrame
    concatenate_df_copy = CONCATENATE_DF.copy()

    # Set the total items for tracking (rows in CONCATENATE_DF)
    if total_items_callback:
        total_items_callback(len(concatenate_df_copy), 0)  # Total to process, initially 0 completed

    # Enrich the DataFrame by merging with pubchem_datas
    enriched_df = concatenate_df_copy.merge(
        globals_vars.pubchem_datas,
        on='INCHIKEY',
        suffixes=('', '_pubchem'),
        how='left'
    )

    # Initialize counters for tracking progress
    processed_items = 0
    total_items = len(enriched_df) if progress_callback else 0

    # Define the columns to update with PubChem data
    columns_to_update = ['INCHI', 'SMILES', 'FORMULA', 'NAME', 'EXACTMASS', 'AVERAGEMASS']

    # Update each column in enriched_df with data from PubChem (if available)
    for col in columns_to_update:
        # Update cells in columns with missing data
        enriched_df[col] = enriched_df[col + '_pubchem'].combine_first(enriched_df[col])

        # Update progress after processing each column
        processed_items += 1
        if progress_callback:
            progress_callback(processed_items)

    # Remove the additional columns created by the PubChem merge
    enriched_df.drop(columns=[col + '_pubchem' for col in columns_to_update], inplace=True)

    if progress_callback:
        progress_callback(len(concatenate_df_copy))

    return enriched_df

