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
    with open("ligands_lipinski_properties.json", "w") as f:
        json.dump(all_props, f, indent=2)

if __name__ == "__main__":
    fetch_all_compound_properties()
    print("All properties fetched and saved to ligands_lipinski_properties.json")
