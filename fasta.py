'''
Will remove duplicates and similar sequences in the fasta files and the dictionary
Duplicates are based on protein sequence to reduce redundancy during alignment
Sequence similarity is based on Levenshtein Edit Distance between strings
'''

import json
from Bio import SeqIO

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

# Will determine what values need to be removed
def remove():
    # Will hold the protein sequences as the key and its header line as the value
    tempdict = {}
    # Used to remove sequences
    removed = {}

    with open("data/ecod.latest.fasta.txt", "rt") as f:
        for reads in SeqIO.parse(f, 'fasta'):
            sequence = str(reads.seq)
            ID = str(reads.id)
            # Remove the entry if the sequence is in the dictionary already
            if sequence in tempdict:
                if sequence not in removed:
                    removed[sequence] = []
                removed[sequence].append((((ID.split('|')[0])[1:]), (ID.split('|')[2])))

            # Add sequence as a key to the dictionary and the ID as a value
            elif sequence not in tempdict:
                # Used as boolean for if they are different enough to add to the dictionary
                OK = True

                if sequence in removed:
                    OK = False

                else:
                    # Get the protein hierarchy level
                    level = "|" + ID.split("|")[2] + "|"
                    # Check the keys in the dictionary of the same level to see the number of differences
                    for seq, id in tempdict.items():
                        # If they are similar, then we will not add it to dictionary (90% sequence identity)
                        if level in id and identity(sequence, seq) >= 0.9:
                            OK = False
                            break

                # Enough differences so add it to the dictionary
                if OK == True:
                    tempdict[reads.seq] = reads.id

                # Not enough differences so we will remove it from the file
                else:
                    if sequence not in removed:
                        removed[sequence] = []
                    removed[sequence].append((((ID.split('|')[0])[1:]), (ID.split('|')[2])))
    return tempdict, removed

# Rewrites the fasta file so that there are no duplicates or similar sequences (< 3 differences)
def removeFasta(tempdict):
    # Write the dictionary values to a fasta file
    with open("data/ecod.latest.fasta.txt", 'w') as f:
         SeqIO.write(tempdict.values(), f, 'fasta')

# Rewrites the json file so that there are no duplicates or similar sequences (< 3 differences)
def removeDict(proteins, removed):
    for protein in removed.values():
        for value in protein:
            # Split the id so we can parse through the proteins dictionary
            x, h, t = value[1].split(".")[0], value[1].split(".")[1], value[1].split(".")[2]
            f = value[1].split(".")[3] if len(value[1].split(".")) == 4 else "None"
            hierarchy = proteins[x][h][t][f]

            # Remove the duplicate proteins in the json file
            id = value[0]
            for i in range(len(hierarchy)):
                if hierarchy[i]["u_id"] == id:
                    del hierarchy[i]
                    break

    # Output the updated file
    with open("proteins.json", "w") as f:
        json.dump(proteins, f, ensure_ascii=False, indent=4)

# Calculates the number of differences between two strings (the distance)
def identity(s1,s2):
    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2, char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1, char1 in enumerate(s1):
            if char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
        distances = newDistances
    return round((len(s2) - distances[-1])/len(s2), 1)

def main():
    proteins = openFile()
    tempdict, removed = remove()
    removeFasta(tempdict)
    removeDict(proteins, removed)
   
if __name__ == "__main__":
    main()
