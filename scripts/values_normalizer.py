
def normalize_values(metadata_dict):
    metadata_dict = normalize_empties(metadata_dict)
    
    metadata_dict = normalize_adduct(metadata_dict)
    metadata_dict = normalize_ionmode(metadata_dict)
    metadata_dict = normalize_retention_time(metadata_dict)
    metadata_dict = normalize_ms_level(metadata_dict)
    metadata_dict = normalize_synonymes(metadata_dict)
    metadata_dict = normalize_formula(metadata_dict)
    metadata_dict = normalize_predicted(metadata_dict)
    metadata_dict = normalize_db_informations(metadata_dict)

    return metadata_dict