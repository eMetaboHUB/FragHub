from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from rdkit import RDLogger, Chem
import deletion_report
import globals_vars
import pandas as pd
import re
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages


def apply_transformations(inchi_smiles):
    """
    Apply transformations to a given InChI or SMILES string.
    :param inchi_smiles: The InChI or SMILES string.
    :return: A dictionary containing the transformed values of the input string.
    """
    # Initialize an empty dictionary for transformed values
    transforms = {}

    # Check if the input string does not contain 'InChI='. If true, some corrections will be applied to the string
    if 'InChI=' not in inchi_smiles:
        inchi_smiles = re.sub(globals_vars.indigo_smiles_correction_pattern, "", inchi_smiles)

    # Check if the modified input string is indeed a string
    if isinstance(inchi_smiles, str):
        # Depending on whether it is in InChI format or SMILES format, appropriate conversion will be applied
        mol = Chem.MolFromInchi(inchi_smiles) if 'InChI=' in inchi_smiles else Chem.MolFromSmiles(inchi_smiles)

        # Checking if the conversion produced a valid result. If so, chemical information is extracted
        if mol is not None:
            transforms = {
                'INCHI': Chem.MolToInchi(mol),
                'INCHIKEY': Chem.MolToInchiKey(mol),
                'SMILES': Chem.MolToSmiles(mol),
                'FORMULA': CalcMolFormula(mol),
            }
        else:
            transforms = {
                'INCHI': '',
                'INCHIKEY': '',
                'SMILES': '',
                'FORMULA': '',
            }

        # Calculating the masses. If any errors occur during this process, default to blank strings
        if transforms:
            mol = Chem.MolFromInchi(transforms['INCHI']) if 'InChI=' in inchi_smiles else Chem.MolFromSmiles(transforms['SMILES'])
            if mol is not None:
                try:
                    transforms['EXACTMASS'] = ExactMolWt(mol)
                    transforms['AVERAGEMASS'] = MolWt(mol)
                except:
                    transforms['EXACTMASS'] = ''
                    transforms['AVERAGEMASS'] = ''
                    return transforms
            else:
                transforms['EXACTMASS'] = ''
                transforms['AVERAGEMASS'] = ''

    # Return the transformations dictionary
    return transforms

def map_transformations(row, unique_transforms):
    """
    Transforms the given row based on the unique transformations specified.
    :param row: A dictionary representing a row of data.
    :param unique_transforms: A dictionary of unique transformations.
    :return: The transformed row.
    """
    # Set original_inchi to the value of 'INCHI' in the row or None if it isn't available.
    original_inchi = row['INCHI'] if pd.notna(row['INCHI']) else None

    # Set original_smiles to the value of 'SMILES' in the row or None if it isn't available.
    original_smiles = row['SMILES'] if pd.notna(row['SMILES']) else None

    # If original_inchi is present and also there in unique_transforms, apply the transformations in the row.
    if original_inchi and original_inchi in unique_transforms:
        for key, value in unique_transforms[original_inchi].items():
            row[key] = value

    # If original_smiles is present and also there in unique_transforms, apply the transformations in the row.
    elif original_smiles and original_smiles in unique_transforms:
        for key, value in unique_transforms[original_smiles].items():
            row[key] = value

    # Return the transformed row
    return row


def mols_derivation_and_calculation(CONCATENATE_DF, output_directory, progress_callback=None, total_items_callback=None,
                                    prefix_callback=None, item_type_callback=None):
    """
    Derives and calculates molecular properties based on unique INCHI and SMILES in the given DataFrame with progress reporting using callbacks.

    :param CONCATENATE_DF: DataFrame containing the INCHI and SMILES columns.
    :param progress_callback: A function to update progress as items are processed.
    :param total_items_callback: A function to set the total number of items to process.
    :param prefix_callback: A function to dynamically set the prefix for the operation (if needed).
    :param item_type_callback: A function to specify the type of items (e.g., "molecules").
    :return: Tuple where the first element is the filtered DataFrame with calculated properties,
             and the second element is another DataFrame with dropped rows.
    """

    # Set the progress-related prefix, if provided
    if prefix_callback:
        prefix_callback("derivation and calculation:")

    # Specify the type of items being processed, if provided
    if item_type_callback:
        item_type_callback("rows")

    # Step 1: Concatenate INCHI and SMILES columns, drop nulls, and get unique values
    unique_inchi_smiles = pd.concat([CONCATENATE_DF['INCHI'], CONCATENATE_DF['SMILES']]).dropna().unique()

    # Set the total number of unique INCHI/SMILES for progress tracking
    if total_items_callback:
        total_items_callback(len(unique_inchi_smiles), 0)  # Total = number of unique molecules, Completed = 0

    # Step 2: Apply transformations for each unique INCHI/SMILES string
    processed_items = 0
    unique_transforms = {}

    for inchi_smiles in unique_inchi_smiles:
        # Apply transformation
        unique_transforms[inchi_smiles] = apply_transformations(inchi_smiles)

        # Update progress after each transformation
        processed_items += 1
        if progress_callback:
            progress_callback(processed_items)

    # Step 3: Map transformations back to the DataFrame
    if prefix_callback:
        prefix_callback("updating rows:")
    results_processed = 0
    total_items = len(CONCATENATE_DF)  # Total rows in the DataFrame
    if total_items_callback:
        total_items_callback(total_items, 0)

    def apply_row_mapping(row):
        nonlocal results_processed
        results_processed += 1

        # Report progress
        if progress_callback:
            progress_callback(results_processed)

        # Apply transformations using the dictionary of unique transformations
        return map_transformations(row, unique_transforms)

    CONCATENATE_DF = CONCATENATE_DF.apply(apply_row_mapping, axis=1)

    # Step 4: Validate the 'INCHIKEY' column using a predefined pattern
    mask = CONCATENATE_DF['INCHIKEY'].str.fullmatch(globals_vars.inchikey_pattern, na=False)

    # Apply the mask to retain only valid rows
    CONCATENATE_DF = CONCATENATE_DF[mask]

    # Store the length of the DataFrame before dropping nulls
    before = len(CONCATENATE_DF)

    # Step 5: Drop null values in critical columns
    critical_columns = ['EXACTMASS', 'AVERAGEMASS', 'SMILES', 'INCHI', 'INCHIKEY']
    rows_to_drop = CONCATENATE_DF[CONCATENATE_DF[critical_columns].isnull().any(axis=1)]  # Rows to drop
    CONCATENATE_DF = CONCATENATE_DF.dropna(subset=critical_columns)  # Filtered DataFrame

    rows_to_drop['DELETION_REASON'] = "spectrum deleted because it has neither inchi nor smiles nor inchikey, even after re calculation"

    # Step 6: Write dropped rows to a CSV file
    # Ensure the target directory exists
    deletion_dir = os.path.join(output_directory, "DELETED_SPECTRUMS")

    # Define file name and path
    deleted_file_path = os.path.join(deletion_dir, "deleted_no_inchi_smiles_inchikey_after_re_calculation.csv")

    # Write the rows_to_drop DataFrame to the CSV file
    rows_to_drop.to_csv(deleted_file_path, index=False, sep='\t', encoding='utf-8')

    del rows_to_drop

    # Calculate the number of missing rows
    after = len(CONCATENATE_DF)
    missing = before - after  # Number of deleted rows
    deletion_report.no_smiles_no_inchi_no_inchikey += missing

    # Return both the filtered DataFrame and the dropped rows
    return CONCATENATE_DF
