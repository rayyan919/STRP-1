from get_smiles import get_CIDs

def write_all_CIDs(filename="InitialCompounds.txt"):
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



