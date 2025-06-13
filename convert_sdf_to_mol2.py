import os
import glob
import subprocess

input_dir = "ligands_3D"
output_folder = "ligands_mol2"

if not os.path.exists(output_folder):
    print("Creating output folder:", output_folder)
    os.makedirs(output_folder)

sdf_files = glob.glob(os.path.join(input_dir, "*.sdf"))

for sdf_file in sdf_files:
    base = os.path.splitext(os.path.basename(sdf_file))[0]
    mol2_file = os.path.join(output_folder, base + ".mol2")
    
    command = ["obabel", "-i", "sdf", sdf_file, "-o", "mol2", "-O", mol2_file]
    print("Converting:", " ".join(command))
    
    subprocess.call(command)

print("Batch conversion complete. MOL2 files saved in:", output_folder)
