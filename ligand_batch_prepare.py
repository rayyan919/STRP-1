import os
import glob
import subprocess

mol2_dir = r"C:\Users\rayya\OneDrive\Documents\ligands\ligands_mol2"
output_dir = r"C:\Users\rayya\OneDrive\Documents\ligands\ligands_pdbqt"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

mgltools_python = r"C:\Program Files (x86)\MGLTools-1.5.7\python.exe"
prepare_script = r"C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py"

mol2_files = glob.glob(os.path.join(mol2_dir, "*.mol2"))
error_files = []

for mol2_file in mol2_files:
    base = os.path.splitext(os.path.basename(mol2_file))[0]
    pdbqt_file = os.path.join(output_dir, base + ".pdbqt")
    
    if not os.path.isfile(mol2_file):
        error_files.append(mol2_file)
        continue

    cmd = [
        mgltools_python,
        prepare_script,
        "-l", os.path.basename(mol2_file),  # pass just filename, not full path
        "-o", pdbqt_file,
        "-A", "hydrogens"
    ]
    print("Processing:", mol2_file)
    
    # Run with cwd set to mol2_file directory so prepare_ligand4.py finds ligand file by name
    subprocess.call(cmd, cwd=mol2_dir)
if error_files:
    print("The following files were not processed due to errors:")
    for error_file in error_files:
        print(error_file)
else:
    print("Batch ligand preparation complete.")
