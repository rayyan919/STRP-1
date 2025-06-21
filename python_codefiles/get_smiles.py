import requests
from glob import glob
import pandas as pd
import time

def get_CIDs(filename = "../ligands_3D/*.sdf"):
    listCID = []

    sdf_files = glob(filename)
    for path in sdf_files:
        filename = path.split('.')[0]
        listCID.append(filename.split('/')[-1].split('_')[-1])
    return listCID

def get_smiles(cid):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/TXT"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            return response.text.strip()
        else:
            return None
    except:
        return None

def get_all_smiles():
    results = []
    listCID = get_CIDs()
    for cid in listCID:
        smiles = get_smiles(cid)
        print(f"{cid}: {smiles}")
        results.append({'CID': cid, 'SMILES': smiles})
        time.sleep(0.2)  # API limits

    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv("../output/ligands_smiles.csv", index=False)

