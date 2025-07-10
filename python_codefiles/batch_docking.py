import os
import subprocess
import csv


vina_path = r"C:\vina\vina_1.2.7_win.exe"
receptor = os.path.abspath("../proteins_pdb/prepared/7rpb.pdbqt")
ligand_dir = os.path.abspath("../antibiotics/antibiotics_pdbqt_minimal")
output_dir = os.path.abspath("../antibiotics/docking_results/docked")
log_dir = os.path.abspath("../antibiotics/docking_results/logs")
results_csv = os.path.abspath("../antibiotics/docking_results/docking_scores.csv")


center_x, center_y, center_z = 0.75, 34.5, -11.01
size_x, size_y, size_z = 30, 30, 30

# Make output directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)


def extract_score_from_pdbqt(pdbqt_file):
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except:
        pass
    return None

# Main docking function
def batch_docking():
    results = []

    for ligand_file in os.listdir(ligand_dir):
        if ligand_file.endswith(".pdbqt"):
            ligand_name = ligand_file[:-6]
            ligand_path = os.path.join(ligand_dir, ligand_file)
            output_path = os.path.join(output_dir, f"{ligand_name}_out.pdbqt")
            log_path = os.path.join(log_dir, f"{ligand_name}.log")

            # Build Vina command
            cmd = [
                vina_path,
                "--receptor", receptor,
                "--ligand", ligand_path,
                "--center_x", str(center_x),
                "--center_y", str(center_y),
                "--center_z", str(center_z),
                "--size_x", str(size_x),
                "--size_y", str(size_y),
                "--size_z", str(size_z),
                "--out", output_path,
                "--exhaustiveness", "8"
            ]

            # Run docking and capture log
            with open(log_path, "w") as log_file:
                result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)

            # Extract binding score
            score = extract_score_from_pdbqt(output_path)
            results.append([ligand_name, score])

    # Save to CSV
    with open(results_csv, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Ligand Name", "Binding Score (kcal/mol)"])
        writer.writerows(results)

    print("✅ Batch docking completed.")
    print(f"📁 Docked files saved in: {output_dir}")
    print(f"📄 Logs saved in: {log_dir}")
    print(f"📊 CSV results saved in: {results_csv}")

if __name__ == "__main__":
    batch_docking()
