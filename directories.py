'''
Goes through the hierarchy of each protein and will create a directory of all proteins within the level
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

def getID(f_dict, path):
    path = os.path.join(path, "sequences")
    u_idList = []
    for protein in f_dict:
        u_idList.append(protein['u_id'])

    makeFasta(u_idList, path)

def makeFasta(u_idList, path):
    with open('../../Downloads/ecod.latest.fasta.txt', 'r') as infile, open(os.path.join(path, 'multifasta.txt'), 'w') as outfile:
        currentLine = ""
        for line in infile:
            if line.startswith('>'):
                for u_id in u_idList:
                    if u_id in currentLine:
                        outfile.write(currentLine)
                currentLine = ""
            currentLine += line

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
                    getID(f_dict, path)

def main():
    proteins = openFile()
    parseDict(proteins)


if __name__ == "__main__":
    main()
