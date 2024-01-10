import pandas as pd
def csv_to_json_processing(FINAL_CSV):
    """
    Convert a CSV file to a list of JSON objects.

    :param FINAL_CSV: pandas DataFrame representing the CSV file
    :type FINAL_CSV: pandas.DataFrame
    :return: list of JSON objects representing each row in the CSV file
    :rtype: list

    """
    json_list = FINAL_CSV.to_dict('records')

    return json_list
