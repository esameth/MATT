'''
Will remove duplicates and similar sequences in the fasta files and the dictionary
Duplicates are based on protein sequence to reduce redundancy during alignment
Sequence similarity is based on Levenshtein Edit Distance between strings
'''

import json
import Levenshtein
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class Protein:
    def __init__(self, level, u_id):
        self.level = level
        self.u_id = u_id

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

# Will determine what values need to be removed
def remove():
    # Used to remove sequences
    removed = {}
    # Key = level, Value = [Sequences]
    temp = {}
    # Holds SeqRecords to write to a fasta file
    records = []

    with open("data/ecod.latest.F40.fasta.copy.fasta", "rt") as input:
        for reads in SeqIO.parse(input, 'fasta'):
            sequence = str(reads.seq)
            ID = str(reads.id)

            # Get the protein hierarchy level
            level = ID.split("|")[2]

            # Used as boolean for if they are different enough to add to the dictionary
            OK = True

            # Increases run time - we only get the Levenshtein ratio of those of the same hierarchy
            if level not in temp:
                temp[level] = []

            if sequence not in temp[level] and sequence not in removed:
                # Check the keys in the dictionary of the same level to see the number of differences
                for seq in temp[level]:
                    # If they are similar, then we will not add it to dictionary (90% sequence identity)
                    if Levenshtein.ratio(sequence, seq) >= 0.9:
                        OK = False
                        break
                # Enough differences so add it to the dictionary
                if OK == True:
                    temp[level].append(sequence)
                    records.append(SeqRecord(Seq(sequence, IUPAC.protein), id=ID, name=ID, description=ID))

            if sequence in temp or OK == False:
                if sequence not in removed:
                    removed[sequence] = []
                removed[sequence].append(Protein(level, ID.split("|")[0]).__dict__)

    # Rewrites the fasta file so that there are no duplicates or similar sequences (< 3 differences)
    with open("fastacopy_test.fasta", 'w') as output:
        SeqIO.write(records, output, 'fasta')

    return removed


# Rewrites the json file so that there are no duplicates or similar sequences (< 3 differences)
def removeDict(proteins, removed):
    for protein in removed.values():
        for value in protein:
            # Split the id so we can parse through the proteins dictionary
            x, h, t = value["level"].split(".")[0], value["level"].split(".")[1], value["level"].split(".")[2]
            f = value["level"].split(".")[3] if len(value["level"].split(".")) == 4 else "None"
            hierarchy = proteins[x][h][t][f]
            # Remove the duplicate proteins in the json file
            id = value["u_id"]
            for i in range(len(hierarchy)):
                if hierarchy[i]["u_id"] == id:
                    del hierarchy[i]
                    break

    # Output the updated file
    with open("proteins.json", "w") as f:
        json.dump(proteins, f, ensure_ascii=False, indent=4)

def main():
    proteins = openFile()
    removed = remove()
    removeDict(proteins, removed)

if __name__ == "__main__":
    main()
