import requests
from glob import glob
import time
from get_smiles import get_CIDs


def fetch_ligand_CIDs(filename):
    listCID = []
    
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                listCID.append(line.strip().split(':')[1].strip())
    return set(listCID)

def get_3d_conformers(cid):
    existing_CID = get_CIDs()

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    downloaded = False
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            if cid not in existing_CID:
                downloaded = True
                with open(f"fetched/Conformer3D_COMPOUND_CID_{cid}.sdf", 'w') as file:
                    file.write(response.text)
                print(f"Downloaded 3D conformer for CID: {cid}")
            else:
                print(f"3D conformer for CID: {cid} already exists.")
        else:
            print(f"Failed to fetch 3D conformer for CID: {cid}, status code: {response.status_code}")
    except requests.RequestException as e:
        print(f"Error fetching 3D conformer for CID: {cid}, error: {e}")
    return downloaded
    
def get_3d_conformers_from_file(filename="Compounds.txt"):
    listCID = fetch_ligand_CIDs(filename)
    count = 0
    for cid in listCID:
        if get_3d_conformers(cid):
            count += 1
            time.sleep(0.6)  # Respect API limits
    print(f"Total 3D conformers downloaded: {count}")

