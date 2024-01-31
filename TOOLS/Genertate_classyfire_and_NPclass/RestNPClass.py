from requests import get
import json
import pandas as pd


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

def getNPClassifFromSmile(smile):
    r = get('https://npclassifier.ucsd.edu/classify?smiles=%s'%(smile))
    if checkRequest(r):
        NPClassif=r.json()
        return NPClassif
    else:
        print("Error : impossible to load NP Classifier stucture from smile %s"%(smile))
        return False 

def readNPClassif(NPCresults):
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
    for ind in dfcpd.index:
        NPCresults=readNPClassif(getNPClassifFromSmile(dfcpd[SMILES_col_name][ind]))
        for key in NPCresults.keys():
            dfcpd[key]=NPCresults[key]
    return dfcpd
        
"""
def addNPclassif(dfcpd):
    for ind in dfcpd.index:
        NPCresults=getNPClassifFromSmile(dfcpd["First(Smiles (canonical SMILES))"][ind])
        if NPCresults == False:
            NPCresults = {'class_results': None, 'superclass_results': None, 'pathway_results': None, 'isglycoside': None}
            #convertlist2singlestring(convertstr2list(NPCresults['class_results']))
        dfcpd["ClassNP"][ind]=convertlist2singlestring(convertstr2list(NPCresults['class_results']))
        dfcpd["SuperClassNP"][ind]=convertlist2singlestring(convertstr2list(NPCresults['superclass_results']))
        dfcpd["PathwayNP"][ind]=convertlist2singlestring(convertstr2list(NPCresults['pathway_results']))
        dfcpd["IsGlycoNP"][ind]=convertlist2singlestring(convertstr2list(NPCresults['isglycoside']))
        print(ind,"/",len(dfcpd.index))
    return dfcpd
"""
def convertstr2list(string):
    
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
    if str(list)!="None":
        return '|'.join(list)
    else :
        return None


## test script
dfcpd = pd.read_csv("/home/solweig/Thèse/chemomaps/pharmakon/input/testPharmakon0706.csv",sep =";")
res=[]
SMILES_col_name = "canonical SMILES"
for ind in range(100):
    r=readNPClassif(getNPClassifFromSmile(dfcpd[SMILES_col_name][ind]))

    res.append(r)
with open("/home/solweig/Thèse/chemomaps/pharmakon/input/NPClassifier/testres.json","w") as jsonFile:
    json.dump(res,jsonFile,indent=4)
    jsonFile.close