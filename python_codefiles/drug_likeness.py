import pandas as pd
import json

# Define the property rules
PROPERTY_RULES = [
    ('Muegge', 'Muegge #violations', '<=', 2),
    ('BA', 'Bioavailability Score', 'range', (0.55, 0.90)),
    ('SA', 'Synthetic Accessibility', '<=', 7),
    ('LogP', 'Consensus Log P', '<=', 4.5),
    ('Sw', 'ESOL Solubility (mg/ml)', '>=', 0.010),
    ('LogS', 'ESOL Log S', 'range', (-6.0, 0.5)),
    ('Fsp³', 'Fraction Csp3', '>', 0.36)
]

def check_property(value, rule):
    """Check if a property value satisfies the given rule"""
    operator, threshold = rule[2], rule[3]
    
    if pd.isna(value):
        return False
    
    if operator == '<':
        return value < threshold
    elif operator == '<=':
        return value <= threshold
    elif operator == '>':
        return value > threshold
    elif operator == '>=':
        return value >= threshold
    elif operator == 'range':
        return threshold[0] <= value <= threshold[1]
    
    return False

def analyze_druglikeness(molecule_ids, csv_file, discarded_file, output_file=None, max_violations=2):
    # Read discarded molecule IDs from JSON
    try:
        with open(discarded_file) as f:
            discarded = [item[1] for item in json.load(f)]  
    except FileNotFoundError:
        print(f"Error: JSON file '{discarded_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return

    # Read CSV data
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: CSV file '{csv_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    # Identify molecule column
    molecule_col = 'Molecule'
    molecule_ids = [str(m) for m in molecule_ids]

    # Filter dataframe for molecules in molecule_ids and not in discarded
    filtered_df = df[(df[molecule_col].astype(str).isin(molecule_ids)) & (~df[molecule_col].astype(str).isin(discarded))]
    if filtered_df.empty:
        print("No molecules found matching criteria (in molecule_ids and not in discarded).")
        return

    # Prepare results
    results = []
    for mol_id in molecule_ids:
        # Skip if molecule is in discarded list
        if mol_id in discarded:
            continue

        # Filter molecule data
        molecule_data = filtered_df[filtered_df[molecule_col].astype(str) == mol_id]
        if molecule_data.empty:
            results.append({
                'Molecule': mol_id,
                'Violations': 'Not found in CSV',
                'Violation Count': None,
                'Percentage': None,
                'Test_in_vitro': 'No'
            })
            continue

        row = molecule_data.iloc[0]
        violations = []
        violation_count = 0

        # Check each property
        for prop in PROPERTY_RULES:
            csv_col = prop[1]
            if csv_col not in df.columns:
                violations.append(f"{prop[0]}: Column '{csv_col}' missing")
                violation_count += 1
                continue

            value = row[csv_col]
            if not check_property(value, prop):
                violations.append(prop[0])
                violation_count += 1

        # Calculate violation percentage
        violation_percentage = (violation_count / len(PROPERTY_RULES)) * 100

        # Determine if should be tested in vitro
        test_in_vitro = "Yes" if violation_count <= max_violations else "No"

        results.append({
            'Molecule': mol_id,
            'Violations': ', '.join(violations) if violations else 'None',
            'Violation Count': violation_count,
            'Percentage': f"{violation_percentage:.2f}%",
            'Test_in_vitro': test_in_vitro
        })

    # Create results dataframe
    results_df = pd.DataFrame(results)
    
    # Save or display results
    if output_file:
        results_df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    else:
        print("\nDrug-likeness Analysis Results:")
        print(results_df.to_string(index=False))

def final_screening(file):
    try:
        df = pd.read_csv(file)
    except FileNotFoundError:
        print(f"Error: CSV file '{file}' not found.")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    # Extract molecule IDs (assuming 'CID' or alternatives)
    molecule_col = 'CID'

    cids = df[molecule_col].astype(str).tolist()

    # Analyze drug-likeness with discarded list from JSON
    analyze_druglikeness(cids, '../input/top_swissadme.csv', '../output/primary_lipinski_discarded_compounds.json', '../output/final_invitro_compounds.csv')

if __name__ == "__main__":
    final_screening('../output/top7_final_hits.csv')