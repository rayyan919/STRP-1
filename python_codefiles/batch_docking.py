import os
import subprocess
import csv

vina_path = r"C:\vina\vina_1.2.7_win.exe"
receptor = os.path.abspath("../proteins_pdb/prepared/7rpb.pdbqt")
ligand_dir = os.path.abspath("../bound_ligand/ligands_pdbqt_minimal")
output_dir = os.path.abspath("../docking_results")
log_dir = os.path.abspath("../docking_results/logs")
results_csv = os.path.abspath("../docking_results/scores/docking_scores.csv")
config_file = os.path.abspath("../input/config.txt")


center_x, center_y, center_z = 0.75, 34.5, -11.01
size_x, size_y, size_z = 25, 25, 25

os.makedirs(output_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(os.path.dirname(results_csv), exist_ok=True)

# Function to extract score from PDBQT
def extract_score_from_pdbqt(pdbqt_file):
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except:
        pass
    return None

def batch_docking():
    results = []

    

 
    cmd = [
        vina_path,
        "--config", config_file,
        "--out", os.path.join(log_dir, "docking.log"),
        "--exhaustiveness=8",
    ] 

    subprocess.call(cmd, cwd=output_dir)
            # Run docking and capture log
            # with open(log_path, "w") as log_file:
            #     result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)

    #         # Extract binding score
    #         score = extract_score_from_pdbqt(output_path)
    #         results.append([ligand_name, score])

    # # Save to CSV
    # with open(results_csv, "w", newline='') as csvfile:
    #     writer = csv.writer(csvfile)
    #     writer.writerow(["Ligand Name", "Binding Score (kcal/mol)"])
    #     writer.writerows(results)

    # print("Batch docking completed.")
    # print(f"Docked files saved in: {output_dir}")
    # print(f"Logs saved in: {log_dir}")
    # print(f"CSV results saved in: {results_csv}")

if __name__ == "__main__":
    batch_docking()
