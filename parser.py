'''
Parse through the fasta file and returns a json file that includes each protein's pdb ID, chain, and hierarchy
Creates a directory of the hierarchy with the fasta file of all proteins in that level
'''

import json
import re
import sys
import os

fastaFile = sys.argv[1]

import fasta
records = fasta.remove()
proteins = {}

# Each protein has a pdb, id and a chain
class Protein:
    def __init__(self, pdb, chain):
        self.pdb = pdb
        self.chain = chain

def addToDict(x, h, t, f, protein):
    # Create X-level (is not dependent on other levels)
    if x not in proteins:
        proteins[x] = {}

    # Create H level which is dependent on X level
    if h not in proteins[x]:
        proteins[x][h] = {}

    # Create T level which is dependent on X level and H level
    if t not in proteins[x][h]:
        proteins[x][h][t] = {}

    # Create F level which is dependent on X level, H level, and T level
    if f not in proteins[x][h][t]:
        proteins[x][h][t][f] = []

    # Add each protein's information to the dict
    proteins[x][h][t][f].append(protein)

# Create the directories
def makeDir(dir):
    # Holds the folder path
    path = ""
    # Holds the ID to grep for the fasta file
    f_id = "|"

    # Go through every single level and create a directory for it
    for level in dir:
        path = os.path.join(path, level)
        f_id += level
        f_id += "." if f_id.count(".") != 3 else "|"

        if not os.path.exists(path):
            os.makedirs(os.path.join(path, 'sequences'))
            os.makedirs(os.path.join(path, 'alignments'))
            # Uses the grep command to pull from the fasta file those that belong to the hierarchy
            # Pipes it to a multifasta file in the sequences directory
            cmd = "LC_ALL=C fgrep '" + f_id + "' -A 1 " + fastaFile + " | grep -v '^--' >> " + os.path.join(path, 'sequences', 'multifasta.fasta')
            os.system(cmd)
            
# Get information needed for json file and creation of directories
def data():
    for record in records:
        record = str(record.id).split('|')
        # Get the values in specific columns
        pdb, f_id, chain = record[1][1:5], record[2].split("."), re.split(':|,', record[3])

        # Split chain if more than one and get only chain names
        chain = ','.join([chain[i] for i in range(len(chain)) if i % 2 == 0])

        # Get the levels
        x, h, t = f_id[0], f_id[1], f_id[2]
        # Some f_id do not have a F group
        f = f_id[3] if len(f_id) == 4 else "None"

        # Add values to dictionary for json file
        protein = Protein(pdb, chain).__dict__
        addToDict(x, h, t, f, protein)

        # Make directory folder for it and make the fasta file
        makeDir([x, h, t, f])

if __name__ == "__main__":
    data()
    # Creates the json file for easy opening of the data
    with open('proteins.json', 'w+', encoding='utf-8') as f:
        json.dump(proteins, f, ensure_ascii=False, indent=4)
