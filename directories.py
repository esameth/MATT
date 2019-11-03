'''
Goes through the hierarchy of each protein and will create a directory of all proteins within the level.
Each folder has a sequence folder and the bottom most child has a fasta file of all the proteins belonging to it
'''

import json
import os

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

# Create the directories based on f_id
def makeDir(path):
    if not os.path.exists(path):
        os.makedirs(os.path.join(path, 'sequences'))

# Uses the grep command to pull from the fasta file those that belong to the hierarchy
# Pipes it to a multifasta file in the parent/sequences directory
def makeFasta(path, f_id):
    with open(os.path.join(path, 'sequences', 'multifasta.txt'), 'a+'):
        cmd = "LC_ALL=C fgrep '|" + f_id + "|' -A 1 data/ecod.latest.fasta.txt | grep -v '^--' >> " + os.path.join(path, 'sequences', 'multifasta.txt')
        os.system(cmd)

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
                    makeFasta(path, f_id)

def main():
    proteins = openFile()
    parseDict(proteins)


if __name__ == "__main__":
    main()
