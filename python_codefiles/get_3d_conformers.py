import requests
import os
import time
from glob import glob
from get_smiles import get_CIDs


def fetch_ligand_CIDs(filename):
    """
    Fetch the list of CIDs from the provided text file.
    The CID format assumed is 'SomeText: CID_Number'.
    
    :param filename: Path to the text file containing CIDs.
    :return: A set of CIDs to avoid duplicates.
    """
    listCID = []
    
    # Read the file and extract the CIDs
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                try:
                    listCID.append(line.strip().split(':')[1].strip())
                except IndexError:
                    print(f"Warning: Skipping invalid line: {line.strip()}")
    
    return set(listCID)


def get_3d_conformers(cid, retries=3):
    """
    Fetch the 3D conformer SDF file for a given CID from PubChem.
    Retries in case of failure or transient errors.

    :param cid: PubChem CID of the compound.
    :param retries: Number of retry attempts in case of failure.
    :return: True if the file was downloaded successfully, False otherwise.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    downloaded = False
    currentCIDs = get_CIDs()
    # Retry logic
    if cid in currentCIDs:
        print(f"3D conformer for CID: {cid} already exists in the current CIDs list.")
        return
    
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                output_dir = '../fetched'
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                file_path = os.path.join(output_dir, f"Conformer3D_COMPOUND_CID_{cid}.sdf")
                if not os.path.exists(file_path):
                    with open(file_path, 'w') as file:
                        file.write(response.text)
                    print(f"Downloaded 3D conformer for CID: {cid}")
                    downloaded = True
                else:
                    print(f"3D conformer for CID: {cid} already exists.")
                break
            else:
                print(f"Failed to fetch 3D conformer for CID: {cid}, status code: {response.status_code}")
        except requests.RequestException as e:
            print(f"Error fetching 3D conformer for CID: {cid}, attempt {attempt + 1}/{retries}. Error: {e}")
            time.sleep(5)  # Wait before retrying

    return downloaded


def get_3d_conformers_from_file(filename="../input/InitialCompounds.txt"):
    """
    Fetch 3D conformers for all CIDs listed in the provided file.
    
    :param filename: Path to the text file containing the CIDs.
    :return: None
    """
    listCID = fetch_ligand_CIDs(filename)
    count = 0
    failed_CIDs = []  # List to track failed downloads

    for cid in listCID:
        if get_3d_conformers(cid):
            count += 1
            time.sleep(0.6)  # Respect PubChem's API rate limits
        else:
            failed_CIDs.append(cid)

    print(f"Total 3D conformers downloaded: {count}")
    if failed_CIDs:
        print(f"Failed to download 3D conformers for the following CIDs: {failed_CIDs}")
        log_failed_CIDs(failed_CIDs)


def log_failed_CIDs(failed_CIDs):
    """
    Logs failed CIDs into a text file for later review.

    :param failed_CIDs: List of CIDs that failed to download.
    :return: None
    """
    with open('failed_downloads.txt', 'w') as log_file:
        for cid in failed_CIDs:
            log_file.write(f"{cid}\n")
    print("Failed CIDs have been logged in 'failed_downloads.txt'")


if __name__ == "__main__":
    get_3d_conformers_from_file()

