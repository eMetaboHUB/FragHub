import concurrent.futures

import pandas as pd
from tqdm import tqdm
import os
import re

def cfm_id_test(mode):
    if mode == "pos":

def charging_datas():
    dataframe = pd.read_csv(r"C:\Users\Axel\Documents\PYTHON\MSP_V3\FragBank\cmf_id_test\report.csv",sep=',', quotechar="\"", encoding="UTF-8")


def CMF_id_test_processing():
    for files in os.listdir(r"..\OUTPUT\CSV"):
        if os.path.isdir(os.path.join(r"..\OUTPUT\CSV",files)):
            if files == "POS":
                mode = "pos"
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    results = list(tqdm(executor.map(cfm_id_test, FINAL_XML), total=len(FINAL_XML), unit=" spectrums", colour="green", desc="\t  converting"))

                final = [res for res in results if res is not None]

                return final # returns the list of different worker executions.

