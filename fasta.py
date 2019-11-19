'''
Will remove duplicates in the fasta files and the dictionary
Duplicates are based on protein sequence to reduce redundancy during alignment
'''
import json

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

def remove(proteins):
    # Will hold the protein sequences as the key and its header line as the value
    # Used to remove duplicate sequences
    tempdict = {}
    removed = []

    with open("data/ecod.latest.fasta.txt", "r") as f:
        # Look at the current and previous line
        current = f.readline()
        for line in f:
            previous = current
            current = line
            # If it's the sequence, then add it as a key to the dictionary and the previous line as a value
            if not current.startswith(">") and current != "\n" and current not in tempdict:
                tempdict[current] = previous

            elif not current.startswith(">") and current != "\n" and current in tempdict:
                removed.append((((previous.split('|')[0])[1:]), (previous.split('|')[2])))

    # Write the dictionary values to a fasta file
    with open("data/ecod.latest.fasta.txt", 'w') as f:
         for key, value in tempdict.items():
             line = value + key
             f.write(line)

    # Split the id so we can parse through the proteins dictionary
    for value in removed:
        x, h, t = value[1].split(".")[0], value[1].split(".")[1], value[1].split(".")[2]
        f = value[1].split(".")[3] if len(value[1].split(".")) == 4 else "None"
        hierarchy = proteins[x][h][t][f]

        # Remove the duplicate proteins in the json file
        id = value[0]
        for i in range(len(hierarchy)):
            if hierarchy[i]["u_id"] == id:
                del hierarchy[i]
                break

    # Output the updated file with pretty JSON
    with open("proteins.json", "w") as f:
        json.dump(proteins, f, ensure_ascii=False, indent=4)

def main():
    proteins = openFile()
    remove(proteins)


if __name__ == "__main__":
    main()
