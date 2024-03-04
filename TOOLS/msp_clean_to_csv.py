import pandas as pd
import re


def msp_clean_to_csv(MSP):
    dictionary = {}
    CSV_df = pd.DataFrame()
    first = True
    empty = False
    if len(MSP) == 0:
        empty = True

    for spectrum in MSP:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        CSV_df = pd.DataFrame.from_dict(dictionary)

    CSV_df.to_csv(r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\COSIN\TEST.csv", sep=";", encoding="UTF-8", index=False)

MSP =  r"C:\Users\Axel\Documents\Présentations\MSP\datas diagram\COSIN\TEST.msp"

with open(MSP, "r", encoding="UTF-8") as buffer:
    content = buffer.read()

content = content.split("\n\n")

msp_clean_to_csv(content)