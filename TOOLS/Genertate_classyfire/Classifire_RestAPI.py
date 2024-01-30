from requests import get
from tqdm import tqdm
import pandas as pd
import json
import time


def checkRequest(r):
    """
    Check the status code of a request.

    :param r: The response object returned from the request.
    :return: True if the status code is 200, False otherwise.
    """
    if r.status_code !=200 :
        print(r.status_code)
        print(r)
        return False
    return True

def get_urlRequest(url):
    """
    Makes a GET request to the specified URL and returns the response content as JSON if the status code is 200.
    If the status code is 429, it sleeps for 60 seconds and retries the request.
    If the status code is neither 200 nor 429, it returns False.

    :param url: The URL to make the GET request to.
    :return: The response content as JSON if the status code is 200, False otherwise.
    """
    r = get(url)
    if r.status_code == 429:
        time.sleep(60)
        r = get_urlRequest(url)
    elif r.status_code == 200:
        Classyfire = r.json()
        return Classyfire
    else:
        return False

def readlvlname(reslvl):
    """
    Extracts the level name from a dictionary.

    :param reslvl: A dictionary containing the level's information.
    :return: The name of the level, or None if not found.
    """
    try:
        res=reslvl["name"]
    except:
        res=None
    return res

def getClassyfireFromInChIKey(InChIKey):
    """
    :param InChIKey: The InChIKey is a unique identifier for chemical substances. It is used to retrieve classification information from the Classyfire database.
    :return: A dictionary containing the classification information for the given InChIKey. The dictionary includes the following keys:
        - "Classyfire_kingdom": The kingdom classification of the substance.
        - "Classyfire_Superclass": The superclass classification of the substance.
        - "Classyfire_class": The class classification of the substance.
        - "Classyfire_Subclass": The subclass classification of the substance.
        - "direct_parent": The direct parent classification of the substance.
    """
    r = get_urlRequest("http://classyfire.wishartlab.com/entities/%s.json?"%(InChIKey))
    print(r)
    if r ==False or r==None:
        return {"Classyfire_kingdom":None, "Classyfire_Superclass":None, "Classyfire_class":None, "Classyfire_Subclass":None, "direct_parent":None}

    res ={}

    if "kingdom" in r:
        res["Classyfire_kingdom"] = readlvlname(r["kingdom"])
    if "superclass" in r:
        res["Classyfire_Superclass"] = readlvlname(r["superclass"])
    if "class" in r:
        res["Classyfire_class"] = readlvlname(r["class"])
    if "subclass" in r:
        res["Classyfire_Subclass"] = readlvlname(r["subclass"])
    if "direct_parent" in r:
        res["direct_parent"] = readlvlname(r["direct_parent"])

    return res

def addClassyfire(dfcpd,InChIKey_col_name):
    """
    Add Classyfire annotations to a dataframe.

    :param dfcpd: A pandas DataFrame containing chemical compounds.
    :param InChIKey_col_name: The name of the column that contains the InChIKeys of the compounds.
    :return: The input dataframe with additional columns containing Classyfire annotations.
    """
    for ind in dfcpd.index:
        ClassyfireResults=getClassyfireFromInChIKey(dfcpd[InChIKey_col_name][ind])
        for key in ClassyfireResults.keys():
            dfcpd[key]=ClassyfireResults[key]
    return dfcpd

res = []

dfcpd = pd.read_excel(r"C:\Users\Axel\Documents\Présentations\FragBank\Publication\Diagrammes\Sunburst\MSP_inchikeys_classyfire.xlsx")
InChIKey_col_name = "INCHIKEY"
pathfinalfile=r"C:\Users\Axel\Documents\Présentations\FragBank\Publication\Diagrammes\Sunburst\MSP_inchikeys_classyfire_complete.json"

lenghtdf = len(dfcpd)
print(lenghtdf)

for ind in tqdm(range(lenghtdf), total=lenghtdf, colour="green"):
    r = getClassyfireFromInChIKey(dfcpd[InChIKey_col_name][ind])
    
    res.append(r)

with open(pathfinalfile,"w") as jsonFile:
    json.dump(res,jsonFile,indent=4)
    jsonFile.close