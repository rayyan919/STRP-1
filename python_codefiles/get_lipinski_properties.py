import requests
import json
import time
from glob import glob
import os
from get_smiles import get_CIDs, get_smiles

def get_lipinski_properties(cid):
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
        props = get_lipinski_properties(cid)
        if props:
            all_props[cid] = props
        else:
            print(f"Failed to fetch properties for CID {cid}")
        time.sleep(0.4)
    # Save detailed JSON
    with open("../output/ligands_lipinski_properties.json", "w") as f:
        json.dump(all_props, f, indent=2)

def num_violations():
    '''
    Filterations will be done based on the following guidelines:
    mw =< 500
    hacceptor =< 10
    hdonor =< 5
    rotatable =< 10
    20=<tpsa=<130
    -1.0 <= xlogp =< 5.6
    '''
    guides = {
        "MolecularWeight": 500,
        "HBondAcceptorCount": 10,
        "HBondDonorCount": 5,
        "RotatableBondCount": 10,
        "TPSA": (20, 130),
        "XLogP": (-1.0, 5.6)
    }
    with open('../output/ligands_lipinski_properties.json', 'r') as f:
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
        if xlogp < guides["XLogP"][0] or xlogp > guides["XLogP"][1]:
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
    with open("../output/primary_lipinski_discarded_compounds.json", "w") as f:
        json.dump(discarded_cids, f, indent=2)
    with open('../output/FilteredCompounds.txt', 'w') as f:
        listCIDs = get_CIDs()
        for cid in listCIDs:
            if cid not in [item[1] for item in discarded_cids]:
                f.write(f"CID:{cid}\n")


def write_filtered_smiles():
    with open('../output/FilteredCompounds.txt', 'r') as f:
        filtered_cids = [line.strip().split(':')[1] for line in f.readlines()]
    
    smiles_file = '../output/ligands_smiles.csv'
    filtered_smiles = []
    
    with open(smiles_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 2:
                cid = parts[0].strip()
                smiles = parts[1].strip()
                if cid in filtered_cids:
                    filtered_smiles.append(smiles)
    with open('../output/FilteredSmiles.txt', 'w') as f:
        for smiles in filtered_smiles:
            f.write(f"{smiles}\n")

def print_all_smiles_only(filename="../input/InitialCompounds.txt"):

    # check whether it is a text file or a csv file
    if filename.endswith('.csv'):
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    smiles = parts[1].strip()
                    print(smiles)
    elif filename.endswith(".txt"):
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split(':')
                print(get_smiles(parts[1]))
                time.sleep(0.2)
    


if __name__ == "__main__":
    # fetch_all_compound_properties()
    # print("All properties fetched and saved to ligands_lipinski_properties.json")
    # early_drug_filteration()
    # write_filtered_smiles()
    print_all_smiles_only()