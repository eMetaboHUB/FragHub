import concurrent.futures
from tqdm import tqdm
import pandas as pd
import subprocess
import os
import re

def cfm_id_test(smile):
    completed = []
    cmd = "docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cd /cfmid/public/; cfm-predict \'" + smile + "\' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 stdout\""
    completed.append(subprocess.run(["powershell", "-Command", cmd], capture_output=True).stdout.decode("UTF-8"))  # Running cfm-id docker into container via Windows Powershell

    cmd = "docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cd /cfmid/public/; cfm-predict \'" + smile + "\' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 stdout\""
    completed.append(subprocess.run(["powershell", "-Command", cmd], capture_output=True).stdout.decode("UTF-8"))  # Running cfm-id docker into container via Windows Powershell

    return completed


def charging_datas():
    dataframe = pd.read_csv(r"C:\Users\Axel\Documents\PYTHON\MSP_V3\FragHub\cmf_id_test\report.csv",sep=',', quotechar="\"", encoding="UTF-8")
    liste_smiles = dataframe['SMILES'].tolist()

    return liste_smiles[0:10]


def CMF_id_test_processing():
    liste_smiles = charging_datas()
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
        results = list(tqdm(executor.map(cfm_id_test, liste_smiles), total=len(liste_smiles), unit=" smiles", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.

if __name__ == "__main__":
    completed = CMF_id_test_processing()
    print(completed)

