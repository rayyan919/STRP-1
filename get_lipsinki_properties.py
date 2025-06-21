import requests
from get_smiles import get_CIDs
import json
import time

def get_lipsinki_properties(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,MolecularWeight,XLogP,RotatableBondCount,TPSA,HBondDonorCount,HBondAcceptorCount/JSON'
    response = requests.get(url)
    if response.status_code == 200:
        data = json.loads(response.text.strip())
        # get the first compound's properties
        compound_props = data["PropertyTable"]["Properties"][0]
        return compound_props
    else:
        print(f"Error fetching properties for CID {cid}")
        return None

def fetch_all_compound_properties():
    listCID = get_CIDs()
    all_props = {}
    for cid in listCID:
        print(f"Fetching properties for CID {cid}...")
        props = get_lipsinki_properties(cid)
        if props:
            all_props[cid] = props
        else:
            print(f"Failed to fetch properties for CID {cid}")
        time.sleep(0.4)
    # Save detailed JSON
    with open("ligands_lipsinki_properties.json", "w") as f:
        json.dump(all_props, f, indent=2)

def num_violations():
    '''
    Filterations will be done based on the following guidelines:
    mw =< 500
    hacceptor =< 10
    hdonor =< 5
    rotatable =< 10
    20=<tpsa=<130
    xlogp =< 5
    '''
    guides = {
        "MolecularWeight": 500,
        "HBondAcceptorCount": 10,
        "HBondDonorCount": 5,
        "RotatableBondCount": 10,
        "TPSA": (20, 130),
        "XLogP": 5
    }
    with open('ligands_lipsinki_properties.json', 'r') as f:
        all_props = json.load(f)
    violation_counts = []
    for cid, props in all_props.items():
        failed_properties = []
        violations = 0
        title = props.get('Title')
        mw = float(props.get('MolecularWeight'))
        hacceptor = props.get('HBondAcceptorCount')
        hdonor = props.get("HBondDonorCount")
        rotatable = props.get("RotatableBondCount")
        tpsa = props.get("TPSA")
        xlogp = props.get("XLogP")
        if mw > guides["MolecularWeight"]:
            violations += 1
            failed_properties.append("MolecularWeight")
        if hacceptor > guides["HBondAcceptorCount"]:
            violations += 1 
            failed_properties.append("HBondAcceptorCount")
        if hdonor > guides["HBondDonorCount"]:
            violations += 1
            failed_properties.append("HBondDonorCount")
        if rotatable > guides["RotatableBondCount"]:
            violations += 1
            failed_properties.append("RotatableBondCount")
        if tpsa < guides["TPSA"][0] or tpsa > guides["TPSA"][1]:
            violations += 1
            failed_properties.append("TPSA")
        if xlogp > guides["XLogP"]:
            violations += 1
            failed_properties.append("XLogP")
        violation_counts.append((title, cid, violations, failed_properties))
    return violation_counts


def early_drug_filteration(allowed=3):
    num_violations_list = num_violations()
    discarded_cids = []
    for title, cid, violations, failed_properties in num_violations_list:
        if violations > allowed:
            discarded_cids.append((title, cid, violations, failed_properties))
    print(f"Discarded {len(discarded_cids)} compounds with more than {allowed} violations.")
    with open("primary_lipsinki_discarded_compounds.json", "w") as f:
        json.dump(discarded_cids, f, indent=2)
    with open('FilteredCompounds.txt', 'w') as f:
        listCIDs = get_CIDs()
        for cid in listCIDs:
            if cid not in [item[1] for item in discarded_cids]:
                f.write(f"CID:{cid}\n")
            

if __name__ == "__main__":
    # fetch_all_compound_properties()
    # print("All properties fetched and saved to ligands_lipinski_properties.json")
    early_drug_filteration()