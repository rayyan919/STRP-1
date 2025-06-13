import os
import glob
import subprocess

input_dir = "ligands_3D"
output_folder = "ligands_mol2_obabel_minimal"
log_folder = "sdf_to_mol2_conversion_logs"

if not os.path.exists(log_folder):
    print("Creating log folder:", log_folder)
    os.makedirs(log_folder)

if not os.path.exists(output_folder):
    print("Creating output folder:", output_folder)
    os.makedirs(output_folder)

sdf_files = glob.glob(os.path.join(input_dir, "*.sdf"))

for sdf_file in sdf_files:
    base = os.path.splitext(os.path.basename(sdf_file))[0]
    mol2_file = os.path.join(output_folder, base + ".mol2")
    
    # Define the log file path for each molecule
    log_file = os.path.join(log_folder, base + "_minimization.log")
    
    # Open the log file for writing
    with open(log_file, 'w') as log:
        
        log.write(f"Minimization log for {base}:\n")

        # Minimal command
        command = [
            "obabel", 
            sdf_file, 
            "-O", mol2_file, 
            "-h"  # Add explicit hydrogens
        ]
        
        # Optimized command
        # command = [
        #     "obabel", 
        #     sdf_file, 
        #     "-O", mol2_file, 
        #     "-h",          # Add explicit hydrogens (important for docking)
        #     "--minimize",  # Energy minimization to ensure correct geometry
        #     "--ff", "MMFF94",  # Use MMFF94 force field for optimization
        #     "--steps", "2500",  # Number of minimization steps
        #     "--log", log_file  # Log the energy minimization process to the log file
        # ]
        
        # Write the command to the log file
        log.write("Command used for conversion:\n")
        log.write(" ".join(command) + "\n")
        
        print("Converting:", " ".join(command))
        
        subprocess.call(command)

print("Batch conversion complete. MOL2 files saved in:", output_folder)
