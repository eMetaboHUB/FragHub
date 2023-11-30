from requests import get
import json
import pandas as pd
import time


def checkRequest(r):
    """[Check if API response is OK (status code 200)]

    Args:
        r ([type]): API response

    Returns:
        [bool]: [False if status code is different than 200 else True]
    """
    if r.status_code !=200 :
        print(r.status_code)
        print(r)
        return False
    return True

def get_urlRequest(url):
    r = get(url)
    if r.status_code == 429:
        time.sleep(30)
        r = get_urlRequest(url)
    elif r.status_code == 200:
        Classyfire = r.json()
        return Classyfire
    else:
        return False

def readlvlname(reslvl):
    try:
        res=reslvl["name"]
    except:
        res=None
    return res

def getClassyfireFromInChIKey(InChIKey):
    r = get_urlRequest("http://classyfire.wishartlab.com/entities/%s.json?"%(InChIKey))
    if r ==False or r==None:
        return {"Classyfire_kingdom":None, "Classyfire_Superclass":None, "Classyfire_class":None, "Classyfire_Subclass":None, "direct_parent":None}
    res ={
        "Classyfire_kingdom":readlvlname(r["kingdom"]), 
        "Classyfire_Superclass":readlvlname(r["superclass"]), 
        "Classyfire_class":readlvlname(r["class"]), 
        "Classyfire_Subclass":readlvlname(r["subclass"]),
        "direct_parent":readlvlname(r["direct_parent"])
    }
    return res

def addClassyfire(dfcpd,InChIKey_col_name):
    for ind in dfcpd.index:
        ClassyfireResults=getClassyfireFromInChIKey(dfcpd[InChIKey_col_name][ind])
        for key in ClassyfireResults.keys():
            dfcpd[key]=ClassyfireResults[key]
    return dfcpd


## test script
#dfcpd = pd.read_csv("/home/solweig/Thèse/chemomaps/pharmakon/input/testPharmakon0706.csv",sep =";")
res=[]
#InChIKey_col_name = "Smiles (InChiKey)"
#pathfinalfile="/home/solweig/Thèse/chemomaps/pharmakon/input/Classyfire/testres.json"


dfcpd = pd.read_csv("/home/solweig/Thèse/chemomaps/mspConvert/mspData/inchikeyUniq.csv",sep =";")
InChIKey_col_name = "INCHIKEY"
pathfinalfile="/home/solweig/Thèse/chemomaps/mspConvert/mspData/inchikeyUniq.csv"
lenghtdf = len(dfcpd)
print(lenghtdf)
for ind in range(lenghtdf):
    r=getClassyfireFromInChIKey(dfcpd[InChIKey_col_name][ind])
    
    res.append(r)
    if ind%50==0:
        print(ind,"/",lenghtdf)
with open(pathfinalfile,"w") as jsonFile:
    json.dump(res,jsonFile,indent=4)
    jsonFile.close