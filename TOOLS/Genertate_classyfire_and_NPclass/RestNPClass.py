from requests import get
from tqdm import tqdm
import pandas as pd
import json

def checkRequest(r):
    """
    :param r: The response object from an HTTP request
    :return: True if the response has a status code of 200, False otherwise
    """
    if r.status_code !=200 :
        print(r.status_code)
        print(r)
        return False
    return True

def getNPClassifFromSmile(smile):
    """
    Retrieves the NPClassifier structure from the given SMILE.

    :param smile: The SMILE string of the desired compound.
    :type smile: str
    :return: The NPClassifier structure if successfully retrieved, False otherwise.
    :rtype: dict or bool
    """
    r = get('https://npclassifier.ucsd.edu/classify?smiles=%s'%(smile))
    if checkRequest(r):
        NPClassif=r.json()
        return NPClassif
    else:
        print("Error : impossible to load NP Classifier stucture from smile %s"%(smile))
        return False 

def readNPClassif(NPCresults):
    """

    Method: readNPClassif

    This method takes in one parameter, NPCresults, and returns a dictionary.

    :param NPCresults: A dictionary containing the classification results for NPClassif
    :return: A dictionary with the updated classification results

    """
    if NPCresults == False:
        return {'class_results': None, 'superclass_results': None, 'pathway_results': None, 'isglycoside': None}
    res = {
        "ClassNP":convertlist2singlestring(convertstr2list(NPCresults['class_results'])),
        "SuperClassNP":convertlist2singlestring(convertstr2list(NPCresults['superclass_results'])),
        "PathwayNP":convertlist2singlestring(convertstr2list(NPCresults['pathway_results'])),
        "IsGlycoNP":convertlist2singlestring(convertstr2list(NPCresults['isglycoside']))
    }
    return res

def addNPClassif(dfcpd, SMILES_col_name="First(Smiles (canonical SMILES))"):
    """
    :param dfcpd: The Pandas DataFrame object containing the data to be classified
    :param SMILES_col_name: The column name in the DataFrame that contains the SMILES string. Default value is "First(Smiles (canonical SMILES))"
    :return: The modified DataFrame with additional columns containing the NPClassif results for each SMILES string

    """
    for ind in dfcpd.index:
        NPCresults=readNPClassif(getNPClassifFromSmile(dfcpd[SMILES_col_name][ind]))
        for key in NPCresults.keys():
            dfcpd[key]=NPCresults[key]
    return dfcpd

def convertstr2list(string):
    """
    Convert a string representation of a list to a list.

    :param string: A string representation of a list.
    :return: A list.
    """
    if type(string)!=str:
        string = str(string).replace("nan","")
    if type(string)==str:
        res=string.replace("[","")
        res=res.replace("]","")
        res=res.replace("'",'')
        res=res.split(", ")
    while '' in res:
        res.remove('')
    return res

def convertlist2singlestring(list):
    """
    Convert a list of strings into a single string separated by a '|'.

    :param list: A list of strings.
    :return: A single string consisting of the elements of the list separated by '|'. Returns None if the input list is None.
    """
    if str(list)!="None":
        return '|'.join(list)
    else :
        return None


## test script
dfcpd = pd.read_excel(r"D:\Axel\DATAS\ISTD_ALL\filtered_DB\Agromix_C18_1.5RT_Jul2023_neg_filtered.xlsx")
res = []
SMILES_col_name = "SMILES"
lenghtdf = len(dfcpd)

for ind in tqdm(range(lenghtdf), total=lenghtdf, colour="green"):
    r = readNPClassif(getNPClassifFromSmile(dfcpd[SMILES_col_name][ind]))
    # print(r)
    res.append(r)

with open(r"D:\Axel\DATAS\ISTD_ALL\filtered_DB\ISTD_neg_smiles_NPclass_complete.json","w") as jsonFile:
    json.dump(res,jsonFile,indent=4)
    jsonFile.close