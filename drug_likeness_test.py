import requests
import numpy as np

def get_molecular_weight(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight/TXT'
    response = requests.get(url)
    if response.status_code == 200:
        return float(response.text.strip())
    else:
        raise Exception(f"Error fetching molecular weight for CID {cid}: {response.status_code}")
    
def get_logp(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/LogP/TXT'
    response = requests.get(url)
    if response.status_code == 200:
        return float(response.text.strip())
    else:
        raise Exception(f"Error fetching LogP for CID {cid}: {response.status_code}")

def get_num_h_bond_donors(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/HBondDonorCount/TXT'
    response = requests.get(url)
    if response.status_code == 200:
        return int(response.text.strip())
    else:
        raise Exception(f"Error fetching number of H-bond donors for CID {cid}: {response.status_code}")

def get_num_h_bond_acceptors(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/HBondAcceptorCount/TXT'
    response = requests.get(url)
    if response.status_code == 200:
        return int(response.text.strip())
    else:
        raise Exception(f"Error fetching number of H-bond acceptors for CID {cid}: {response.status_code}")


def lipsinki_rule_of_five(cid):
    mw = get_molecular_weight(cid)                          # Flag if Molecular Weight > 500
    logp = get_logp(cid)                                    # Flag if LogP value > 5 
    h_bond_donors = get_num_h_bond_donors(cid)              # Flag if H-bond donors > 5
    h_bond_acceptors = get_num_h_bond_acceptors(cid)        # Flag if H-bond acceptors > 10





        
