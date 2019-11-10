'''
Goes through the json file and parses the pdb id to get which folder the pdb file would lie in.
Creates a temp file of the pdb file paths that will be used as input to MATT to create alignments for each level
The temp file is deleted when done
'''

import json
import os
import tempfile

# Path to pdb files
path = "/data/pdb"

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

# Parse through the dictionary and pull the pdb id
def parse(proteins):
    for x, x_dict in proteins.items():
        XpdbList = []
        for h, h_dict in x_dict.items():
            HpdbList = []
            for t, t_dict in h_dict.items():
                TpdbList = []
                for f, f_dict in t_dict.items():
                    FpdbList = []
                    for protein in f_dict:
                        FpdbList.append(protein["pdb"])
                        TpdbList.append(protein["pdb"])
                        HpdbList.append(protein["pdb"])
                    makeTempFile(FpdbList)
                makeTempFile(TpdbList)
            makeTempFile(HpdbList)
        makeTempFile(XpdbList)

# Create a temporary file that will be passed into MATT
def makeTempFile(pdbList):
    with tempfile.TemporaryFile(mode = 'w+') as fp:
        for pdb in pdbList:
            fp.write(os.path.join(path, pdb[1:3], "pdb" + pdb + ".ent.gz") + "\n")
            print(os.path.join(path, pdb[1:3], "pdb" + pdb + ".ent.gz"))

def main():
    proteins = openFile()
    parse(proteins)

if __name__ == "__main__":
    main()
