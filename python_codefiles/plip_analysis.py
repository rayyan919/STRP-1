from plip.structure.preparation import PDBComplex
from openbabel import pybel
import os

def pdbqt_to_pdb(input_file, output_file):
    mol_gen = pybel.readfile("pdbqt", input_file)
    mol = next(mol_gen, None)
    if mol is None:
        raise ValueError(f"Could not read ligand file: {input_file}")
    mol.OBMol.AddHydrogens()
    mol.OBMol.PerceiveBondOrders()
    mol.write("pdb", output_file, overwrite=True)

def create_complex(protein_pdb, ligand_pdbqt, name, output_dir):
    # Create required directories
    os.makedirs(output_dir, exist_ok=True)
    temp_dir = os.path.join(output_dir, "temp_files")
    os.makedirs(temp_dir, exist_ok=True)

    output_pdb = os.path.join(output_dir, f"{name}_complex.pdb")
    ligand_pdb = os.path.join(temp_dir, f"{name}_temp_ligand.pdb")
    
    pdbqt_to_pdb(ligand_pdbqt, ligand_pdb)

    with open(output_pdb, "w") as out_f:
        with open(protein_pdb, "r") as p_f:
            for line in p_f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    out_f.write(line)
        with open(ligand_pdb, "r") as l_f:
            for line in l_f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    out_f.write(line)
        out_f.write("END\n")

    os.remove(ligand_pdb)
    return output_pdb


def summarize_unpaired_atoms(atom_list, label):
    logs = []
    if not atom_list:
        return [f"No unpaired {label}."]
    else:
        logs.append(f"\nUnpaired {label} atoms in ligand:")
        for atom in atom_list:
            logs.append(f"  Atom index: {atom.idx}")
            logs.append(f"  Atomic number: {atom.atomicnum}")
            logs.append(f"  Type: {atom.type}")
            logs.append(f"  Coordinates: {atom.coords}")
            logs.append(f"  Partial charge: {atom.partialcharge:.3f}")
            logs.append("")
        return logs

def analyze_complex(pdb_path, name, output_dir):
    prot = PDBComplex()
    prot.load_pdb(pdb_path)
    try:
        prot.analyze()
    except Exception as e:
        print(f"Error analyzing {name}: {e}")
        return

    output_file = os.path.join(output_dir, f"{name}.txt")

    with open(output_file, 'w', encoding='utf-8') as file:
        for lig in prot.ligands:
            lig_id = f"{lig.hetid}:{lig.chain}:{lig.position}"
            inter = prot.interaction_sets[lig_id]


            file.write(f"\n=== Ligand: {name} {lig.hetid} (Chain {lig.chain}, Res {lig.position}) ===\n")
            file.write("Residues involved in binding:\n-----------------------------\n")

            # 1. Hydrogen Bonds
            if inter.hbonds_ldon or inter.hbonds_pdon:
                file.write("\nHydrogen Bonds:\n")
                for hbond in inter.hbonds_ldon:
                    file.write(f"  Ligand -> Protein | {hbond.restype} {hbond.resnr} (Chain {hbond.reschain}) | Dist A-H: {hbond.distance_ah:.2f}Å | Angle: {hbond.angle:.1f}°\n")
                for hbond in inter.hbonds_pdon:
                    file.write(f"  Protein -> Ligand | {hbond.restype} {hbond.resnr} (Chain {hbond.reschain}) | Dist A-H: {hbond.distance_ah:.2f}Å | Angle: {hbond.angle:.1f}°\n")
            else:
                file.write("No hydrogen bonds found.\n")

            # 2. Salt Bridges
            if inter.saltbridge_lneg or inter.saltbridge_pneg:
                file.write("\nSalt Bridges:\n")
                for sb in inter.saltbridge_lneg + inter.saltbridge_pneg:
                    file.write(f"  {sb.restype} {sb.resnr} (Chain {sb.reschain}) <-> Ligand | Distance: {sb.distance:.2f}Å\n")
            else:
                file.write("No salt bridges found.\n")

            # 3. Hydrophobic Contacts
            if inter.hydrophobic_contacts:
                file.write("\nHydrophobic Interactions:\n")
                for contact in inter.hydrophobic_contacts:
                    file.write(f"  {contact.restype} {contact.resnr} (Chain {contact.reschain}) <-> Ligand | Dist: {contact.distance:.2f}Å\n")
            else:
                file.write("No hydrophobic contacts found.\n")

            # 4. Pi-Stacking
            if inter.pistacking:
                file.write("\nPi-Stacking Interactions:\n")
                for pi in inter.pistacking:
                    file.write(f"  {pi.restype} {pi.resnr} (Chain {pi.reschain}) <-> Ligand | Offset: {pi.offset:.2f}Å\n")
            else:
                file.write("No pi-stacking interactions found.\n")

            # 5. Pi-Cation
            if inter.pication_laro or inter.pication_paro:
                file.write("\nPi-Cation Interactions:\n")
                for pc in inter.pication_laro + inter.pication_paro:
                    file.write(f"  {pc.restype} {pc.resnr} (Chain {pc.reschain}) <-> Ligand | Distance: {pc.distance:.2f}Å\n")
            else:
                file.write("No pi-cation interactions found.\n")

            # 6. Metal Complexes
            if inter.metal_complexes:
                file.write("\nMetal Complexes:\n")
                for mc in inter.metal_complexes:
                    file.write(f"  {mc.restype} {mc.resnr} (Chain {mc.reschain}) <-> Ligand\n")
            else:
                file.write("No metal complexes found.\n")

            # 7. Water Bridges
            if inter.water_bridges:
                file.write("\nWater Bridges:\n")
                for wb in inter.water_bridges:
                    file.write(f"  {wb.restype} {wb.resnr} (Chain {wb.reschain}) <-> Ligand | Water idx: {wb.water.idx}\n")
            else:
                file.write("No water bridges found.\n")

            # 8. Halogen Bonds
            if inter.halogen_bonds:
                file.write("\nHalogen Bonds:\n")
                for hb in inter.halogen_bonds:
                    file.write(f"  {hb.restype} {hb.resnr} (Chain {hb.reschain}) <-> Ligand | Distance: {hb.distance:.2f}Å\n")
            else:
                file.write("No halogen bonds found.\n")

            # 9. Unpaired Donors and Acceptors
            file.write("\n==============================\n")
            hba_logs = summarize_unpaired_atoms(inter.unpaired_hba, "HBA (acceptors)")
            for log in hba_logs:
                file.write(log+"\n")
            hbd_logs = summarize_unpaired_atoms(inter.unpaired_hbd, "HBD (donors)")
            for log in hbd_logs:
                file.write(log+"\n")
            file.write("==============================\n")

def analyze_all(ligands_dir, protein_file, output_dir):
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "analysis_logs"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "temp_files"), exist_ok=True)

    for ligand_file in os.listdir(ligands_dir):
        if ligand_file.endswith(".pdbqt"):
            if "antibiotic" in ligands_dir:
                parts = ligand_file.strip().split('_')
                name = parts[0]+'_'+parts[1]
                output_dir_new = output_dir_new = os.path.join(output_dir, 'complexes', 'antibiotic_controls')
                log_dir = os.path.join(output_dir, "analysis_logs", "antibiotic_controls")
            else:
                parts = ligand_file.strip().split('_')
                name = parts[-2].split('.')[0]
                rank = os.path.basename(ligands_dir)
                output_dir_new = os.path.join(output_dir, "complexes", f"{rank}_test_ligands_complexes")
                log_dir = os.path.join(output_dir, "analysis_logs", f"{rank}_test_ligands_complexes")

            ligand_path = os.path.join(ligands_dir, ligand_file)
            complex_pdb = create_complex(protein_file, ligand_path, name, output_dir_new)
            os.makedirs(log_dir, exist_ok=True)

            analyze_complex(complex_pdb, name, log_dir)

if __name__ == "__main__":
    # protein_file = os.path.abspath("../proteins_pdb/prepared/7rpb_complex.pdb")
#     antibiotic_dir = os.path.abspath("../docking_results/antibiotic_controls/docked")
#     meropenem_dir = os.path.abspath("../bound_ligand_meropenem/reddocking_results/docked")
#     top_dir = os.path.abspath("../docking_results/tests_docked/docked/top")
#     moderate_dir = os.path.abspath("../docking_results/tests_docked/docked/moderate")
#     output_dir = os.path.abspath("../output/plip_analysis")

    # analyze_all(antibiotic_dir, protein_file, output_dir)
    # analyze_all(meropenem_dir, protein_file, output_dir)
    # analyze_all(top_dir, protein_file, output_dir)
    # analyze_all(moderate_dir, protein_file, output_dir)
    # print("All logs and complexes created and organized!")

