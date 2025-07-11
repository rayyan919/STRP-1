import os
import json

output_dir = os.path.abspath('../output/screening_rankings')

def create_primary_rankings(control='../docking_results/antibiotic_controls/docking_scores.csv', 
                            test='../docking_results/tests_docked/docking_scores.csv'):
    threshold = []
    props_file = os.path.abspath('../output/ligands_lipinski_properties.json')
    os.makedirs(output_dir, exist_ok=True)

    with open(control, 'r') as controlfile:
        next(controlfile)
        for line in controlfile:
            parts = line.strip().split(',')
            threshold.append(float(parts[1]))
    threshold = sorted(threshold)
    
    top_threshold = threshold[1]
    moderate_threshold = threshold[4]
    bottom_threshold = threshold[-1]

    top_file = output_dir+'/top_hits.csv'
    moderate_file = output_dir+'/moderate_hits.csv'
    bottom_file = output_dir+'/bottom_hits.csv'
    eliminate_file = output_dir+'/eliminated.csv'

    top_hits = []
    moderate_hits = []
    bottom_hits = []
    eliminated = []

    with open(props_file, 'r') as props:
        properties = json.load(props)

    with open(test, 'r') as file:
        next(file)
        for line in file:
            ID, affinity = line.strip().split(',')
            affinity = float(affinity)
            CID = ID.strip().split('_')[-1]

            if CID not in properties:
                print(f"{CID} not found in file, skipping it...")
                continue
            if ',' in properties[CID]["Title"]:
                name = '"'+ properties[CID]["Title"]+'_'+CID+'"'
            else:
                name = properties[CID]["Title"]+'_'+CID

            if affinity <= top_threshold:
                top_hits.append((name, affinity))
            elif affinity <= moderate_threshold:
                moderate_hits.append((name, affinity))
            elif affinity <= bottom_threshold:
                bottom_hits.append((name, affinity))
            else:
                eliminated.append((name, affinity))

    with open(top_file, 'w') as top, open(moderate_file, 'w') as moderate, open(bottom_file, 'w') as bottom, open(eliminate_file, 'w') as eliminate:

        top.write("Ligand Name, Binding Affinity (kcal/mol)\n")
        moderate.write("Ligand Name, Binding Affinity (kcal/mol)\n")
        bottom.write("Ligand Name, Binding Affinity (kcal/mol)\n")
        eliminate.write("Ligand Name, Binding Affinity (kcal/mol)\n")

        for th in top_hits:
            top.write(f"{th[0]},{th[1]}\n")

        for mh in moderate_hits:
            moderate.write(f"{mh[0]},{mh[1]}\n")

        for bh in bottom_hits:
            bottom.write(f"{bh[0]},{bh[1]}\n")

        for e in eliminated:
            eliminate.write(f"{e[0]},{e[1]}\n")

    print("All ranks assigned")

if __name__ == "__main__":
    create_primary_rankings()

