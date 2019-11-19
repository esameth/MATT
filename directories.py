'''
Goes through the hierarchy of each protein and will create a directory of all proteins within the level.
Each folder has a sequence folder and the bottom most child has a fasta file of all the proteins belonging to it
'''

import json
import os
import tempfile

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

# Create the directories based on f_id
def makeDir(path):
    if not os.path.exists(path):
        os.makedirs(os.path.join(path, 'sequences'))
        os.makedirs(os.path.join(path, 'alignments'))

# Uses the grep command to pull from the fasta file those that belong to the hierarchy
# Pipes it to a temporary file
def makeFasta(path, f_id, proteins):
    # Will hold the protein sequences as the key and its header line as the value
    # Used to remove duplicate sequences
    tempdict = {}
    removed = []

    # Create a temporary file for the grep output
    with tempfile.NamedTemporaryFile(mode="w+t") as temp:
        cmd = "LC_ALL=C fgrep '|" + f_id + "|' -A 1 data/ecod.latest.F40.fasta.txt | grep -v '^--' >> " + temp.name
        os.system(cmd)
        temp.seek(0)

        # Look at the current and previous line
        current = temp.readline()
        for line in temp:
            previous = current
            current = line
            # If it's the sequence, then add it as a key to the dictionary and the previous line as a value
            if not current.startswith(">") and current not in tempdict:
                tempdict[current] = previous
            elif not current.startswith(">") and current in tempdict:
                removed.append((previous.split('|')[0])[1:])

    # Write the dictionary values to a fasta file
    with open(os.path.join(path, 'sequences', 'sequences.fasta'), 'w+') as f:
        for key, value in tempdict.items():
            line = value + key
            f.write(line)

    # Split the id so we can parse through the proteins dictionary
    x_level, h_level, t_level = f_id.split(".")[0], f_id.split(".")[1], f_id.split(".")[2]
    f_level = f_id.split(".")[3] if len(f_id.split(".")) == 4 else "None"

    hierarchy = proteins[x_level][h_level][t_level][f_level]

    # Remove the duplicate proteins in the json file
    for id in removed:
        for i in range(len(hierarchy)):
            if hierarchy[i]["u_id"] == id:
                del hierarchy[i]
                break

    # Output the updated file with pretty JSON
    with open("proteins.json", "w") as f:
        json.dump(proteins, f, ensure_ascii=False, indent=4)

# Go through the dictionary
def parseDict(proteins):
    for x, x_dict in proteins.items():
        makeDir(os.path.join(x))
        for h, h_dict in x_dict.items():
            makeDir(os.path.join(x, h))
            for t, t_dict in h_dict.items():
                makeDir(os.path.join(x, h, t))
                for f, f_dict in t_dict.items():
                    path = os.path.join(x, h, t, f)
                    makeDir(path)
                    # Replace the path with . so we can look at the f_id (ex: 1/1/1/2 --> 1.1.1.2)
                    f_id = path.replace('/', '.')
                    # The dictionary had "None" as a placeholder for those without an f-level (ex: 1.1.3)
                    if f_id.rsplit('.', 1)[1] == "None":
                        f_id = f_id.rsplit('.', 1)[0]
                    makeFasta(path, f_id, proteins)

def main():
    proteins = openFile()
    parseDict(proteins)


if __name__ == "__main__":
    main()
