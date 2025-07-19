import json
import csv
import pandas as pd

# === CONFIGURATION === #

# Key active-site or functional residues (resname + resnumber, e.g. ARG261)
KEY_RESIDUES = {'SER81', 'TRP221', 'ARG261', 'TYR112', 'MET223', 'VAL130', 'LYS84'}

# Weights for different interaction types (base importance)
INTERACTION_WEIGHTS = {
    "H_pdon": 3.0,
    "H_ldon": 3.0,
    "Salt": 4.0,
    "Hydrophobic": 1.5,
    "Pi": 2.0,
    "Pication": 2.5
}

# Distance and angle thresholds for optimal hydrogen bonding
MIN_DIST = 1.5
MAX_DIST = 4.0
IDEAL_ANGLE = 180
ANGLE_TOL = 60

# === FUNCTIONS === #

def compute_interaction_score(interaction, key_residues, weights):
    score = 0
    matched_residues = set()

    for itype, items in interaction.items():
        if itype not in weights:
            continue
        for i in items:
            res = i.get('res')
            resnr = i.get('resnr')
            res_id = f"{res.upper()}{resnr[:-1]}" if res and resnr else None

            if res_id not in key_residues:
                continue

            base_weight = weights[itype]
            dist = i.get('distance')
            angle = i.get('angle')

            # Distance-based contribution
            if dist and isinstance(dist, (float, int)):
                dist_score = max(0, 1 - ((dist - MIN_DIST) / (MAX_DIST - MIN_DIST)))
            else:
                dist_score = 1

            # Angle-based contribution (if applicable)
            if angle and isinstance(angle, (float, int)):
                angle_score = max(0, 1 - abs(angle - IDEAL_ANGLE) / ANGLE_TOL)
            else:
                angle_score = 1

            interaction_score = base_weight * dist_score * angle_score
            score += interaction_score

            descriptor = f"{res_id} ({itype}"
            if dist:
                descriptor += f", {dist:.2f} Å"
            if angle:
                descriptor += f", {angle:.1f}°"
            descriptor += ")"

            matched_residues.add(descriptor)

    return score, matched_residues


def rank_top_compounds(energy_csv_path, interaction_json_path, output_csv_path, top_n=5):
    # Load data
    df = pd.read_csv(energy_csv_path)
    with open(interaction_json_path, 'r') as f:
        interaction_data = json.load(f)

    records = []

    for _, row in df.iterrows():
        name = row['Ligand Name']
        affinity = float(row['Binding Affinity (kcal/mol)'])
        cid = name.split("_")[-1]

        if cid not in interaction_data:
            continue

        score, residues = compute_interaction_score(
            interaction=interaction_data[cid],
            key_residues=KEY_RESIDUES,
            weights=INTERACTION_WEIGHTS
        )

        records.append({
            "Ligand Name": name,
            "CID": cid,
            "Binding Affinity": affinity,
            "Interaction Score": round(score, 3),
            "Key Interactions": "; ".join(sorted(residues))
        })

    # Sort by interaction score (then affinity)
    top_records = sorted(records, key=lambda x: (-x['Interaction Score'], x['Binding Affinity']))[:top_n]
    output_df = pd.DataFrame(top_records)
    output_df.to_csv(output_csv_path, index=False)
    return output_df

# === USAGE EXAMPLE (you can modify this part in your script) === #
top_df = rank_top_compounds("../docking_results/antibiotic_controls/docking_scores.csv", "../output/plip_analysis/json_data/antibiotic_controls.json", "ranked_antibiotics.csv", 7)
print(top_df)

