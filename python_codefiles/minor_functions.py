from get_smiles import get_CIDs
import requests
import os
import time
import json

def write_all_CIDs(filename="../input/InitialCompounds.txt"):
    # Check whether the file is empty or not if not empty then empty it first and then write all CIDs
    try:
        with open(filename, 'r') as file:
            if file.read().strip():
                print(f"{filename} is not empty. Emptying the file before writing all CIDs.")
                with open(filename, 'w') as empty_file:
                    empty_file.write("")
    except FileNotFoundError:
        print(f"{filename} does not exist. A new file will be created.")

    existing_CID = get_CIDs()
    with open(filename, 'w') as file:
        for cid in existing_CID:
            file.write(f"CID:{cid}\n")
    print(f"All CIDs written to {filename}")

def fetch_smiles(CID):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/property/SMILES/txt"
    print(f"Fetching SMILES for {CID}")

    response = requests.get(url, timeout=5)
    return response


def write_selected_smiles(file_dir, res_file):

    input = os.path.abspath(f"../input/{res_file}")
    files = sorted(os.listdir(file_dir))
    with open(input, 'w') as f:
        for filename in files:
            CID, _ = os.path.splitext(filename)
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/property/SMILES/txt"
            print(f"Fetching SMILES for {CID}")

            try:
                response = requests.get(url, timeout=5)
                if response.status_code == 200:
                    smiles = response.text.strip()
                else:
                    print(f"failed to fetch for {CID}")
                    continue
            except:
                print(f"failed to fetch for {CID}")
                continue
            f.write(f"{smiles} {CID}\n")
            time.sleep(0.2)
        
    print(f"Wrote file {res_file} in input directory")


def create_analysis_json(input_dir, output_dir):

    from plip.structure.preparation import PDBComplex
    from openbabel import pybel

    output_file = os.path.join(output_dir, f"{os.path.basename(input_dir)}.json")

    data = {}

    for pdb_file in os.listdir(input_dir):
        root = os.path.abspath(f"../output/plip_analysis/complexes/{os.path.basename(input_dir)}")
        pdb_path = os.path.join(root, pdb_file)
        prot = PDBComplex()
        prot.load_pdb(pdb_path)
        name = os.path.basename(pdb_path).strip().split('_')[0]
        try:
            prot.analyze()
        except Exception as e:
            print(f"Error analyzing {name}: {e}")
            return
        
    
        for lig in prot.ligands:
            lig_id = f"{lig.hetid}:{lig.chain}:{lig.position}"
            inter = prot.interaction_sets[lig_id]
    
            data[name] = {}
            data[name]["H_ldon"] = []
            data[name]["H_pdon"] = []
            data[name]["Salt"] = []
            data[name]["Hydrophobic"] = []
            data[name]["Pi"] = []
            data[name]["Pication"] = []
            data[name]["Metal"] = []
            data[name]["Water"] = []
            data[name]["Halogen"] = []
            data[name]["HBA"] = []
            data[name]["HBD"] = []

            if inter.hbonds_ldon or inter.hbonds_pdon:
                for hbond in inter.hbonds_ldon:
                    data[name]["H_ldon"].append({
                        "res": hbond.restype,
                        "resnr": f"{hbond.resnr}{hbond.reschain}",
                        "distance": float(f"{hbond.distance_ah:.2f}"),
                        "angle": float(f"{hbond.angle:.1f}")
                    })
                for hbond in inter.hbonds_pdon:
                    data[name]["H_pdon"].append({
                        "res": hbond.restype,
                        "resnr": f"{hbond.resnr}{hbond.reschain}",
                        "distance": float(f"{hbond.distance_ah:.2f}"),
                        "angle": float(f"{hbond.angle:.1f}")
                    })
            
            if inter.saltbridge_lneg or inter.saltbridge_pneg:
                for sb in inter.saltbridge_lneg + inter.saltbridge_pneg:
                    data[name]["Salt"].append({
                        "res": sb.restype,
                        "resnr": f"{sb.resnr}{sb.reschain}",
                        "distance": float(f"{sb.distance:.2f}")
                    })
            
            if inter.hydrophobic_contacts:
                for contact in inter.hydrophobic_contacts:
                    data[name]["Hydrophobic"].append({
                        "res": contact.restype,
                        "resnr": f"{contact.resnr}{contact.reschain}",
                        "distance": float(f"{contact.distance:.2f}")
                    })
            
            if inter.pistacking:
                for pi in inter.pistacking:
                    data[name]["Pi"].append({
                        "res": pi.restype,
                        "resnr": f"{pi.resnr}{pi.reschain}",
                        "offset": float(f"{pi.offset:.2f}")
                    })
            
            if inter.pication_laro or inter.pication_paro:
                for pc in inter.pication_laro + inter.pication_paro:
                    data[name]["Pication"].append({
                        "res": pc.restype,
                        "resnr": f"{pc.resnr}{pc.reschain}",
                        "distance": float(f"{pc.distance:.2f}")
                    })
            
            if inter.metal_complexes:
                for mc in inter.metal_complexes:
                    data[name]["Metal"].append({
                        "res": mc.restype,
                        "resnr": f"{mc.resnr}{mc.reschain}"
                    })

            if inter.water_bridges:
                for wb in inter.water_bridges:
                    data[name]["Water"].append({
                        "res" : wb.restype,
                        "resnr": f"{wb.restype}{wb.resnr}",
                        "idx": wb.water.idx
                    })
            
            if inter.halogen_bonds:
                for hb in inter.halogen_bonds:
                    data[name]["Halogen"].append({
                        "res":hb.restype,
                        "resnr": f"{hb.restype}{hb.resnr}",
                        "distance": f"{hb.distance:.2f}"
                    })
            
            if inter.unpaired_hba:
                for atom in inter.unpaired_hba:
                    data[name]["HBA"].append({
                        "idx": atom.idx,
                        "num": atom.atomicnum,
                        "type": atom.type,
                        "coords": atom.coords,
                        "pc": f"{atom.partialcharge:.3f}"
                    })
            
            if inter.unpaired_hbd:
                for atom in inter.unpaired_hbd:
                    data[name]["HBD"].append({
                        "idx": atom.idx,
                        "num": atom.atomicnum,
                        "type": atom.type,
                        "coords": atom.coords,
                        "pc": f"{atom.partialcharge:.3f}"
                    })
    with open(output_file, "w") as f:    
        json.dump(data, f, indent=2)
        print(f"Wrote all information to {output_file}")

                
import pandas as pd

def compare_ligands_vs_controls(compounds_file, controls_file, output_csv="top7_final_hits.csv", top_n=7):
    # Load the two CSVs
    df_hits = pd.read_csv(compounds_file)
    df_controls = pd.read_csv(controls_file)

    # Compute a composite score (weighted combination)
    # Normalization (z-score like ranking using min-max scaling)
    def normalize(series):
        return (series - series.min()) / (series.max() - series.min())

    df_hits["Affinity Norm"] = normalize(-df_hits["Binding Affinity"])
    df_hits["Interaction Norm"] = normalize(df_hits["Interaction Score"])
    df_hits["Composite Score"] = df_hits["Affinity Norm"] * 0.4 + df_hits["Interaction Norm"] * 0.6

    df_controls["Affinity Norm"] = normalize(-df_controls["Binding Affinity"])
    df_controls["Interaction Norm"] = normalize(df_controls["Interaction Score"])
    df_controls["Composite Score"] = df_controls["Affinity Norm"] * 0.4 + df_controls["Interaction Norm"] * 0.6

    # Select top N ligands only from hits, sorted by composite score
    top_ligands = df_hits.sort_values(by="Composite Score", ascending=False).head(top_n)

    # Output
    top_ligands[[
        "Ligand Name", "CID", "Binding Affinity", "Interaction Score",
        "Key Interactions", "Composite Score"
    ]].to_csv(output_csv, index=False)

    return top_ligands.reset_index(drop=True)

            

if __name__ == "__main__":
    # file_dir1 = os.path.abspath("../output/plip_analysis/complexes/antibiotic_controls")
    # file_dir2 = os.path.abspath("../output/plip_analysis/complexes/top_test_ligands_complexes")
    # output_dir = os.path.abspath("../output/plip_analysis/json_data")
    # create_analysis_json(file_dir1, output_dir)
    compare_ligands_vs_controls("ranked_compounds.csv", "ranked_antibiotics.csv")

    # print(os.listdir(file_dir1))
    # res_file1 = "moderate_smiles.txt"
    # res_file2 = "top_smiles.txt"
    # write_selected_smiles(file_dir1, res_file1)
    # write_selected_smiles(file_dir2, res_file2)


