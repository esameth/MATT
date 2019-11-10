'''
Goes through the json file and parses the pdb id to get which folder the pdb file would lie in.
Creates a temp file of the pdb file paths that will be used as input to MATT to create alignments for each level
The temp file is deleted when done
'''

import json
import os
import tempfile

path = "/data/pdb"

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

def parse(proteins):
    for x, x_dict in proteins.items():
        Hpath = []
        for h, h_dict in x_dict.items():
            Tpath = []
            Hpath.append(os.path.join(x, h,"sequences", "multifasta.txt"))
            HpathP = os.path.join(x, "sequences", "multifasta.txt")
            for t, t_dict in h_dict.items():
                pdbList = []
                Tpath.append(os.path.join(x, h, t, "sequences", "multifasta.txt"))
                TpathP = os.path.join(x, h, "sequences", "multifasta.txt")
                for f, f_dict in t_dict.items():
                    for protein in f_dict:
                        pdbList.append(protein["pdb"])
                print(pdbList)
                makeFile(pdbList)

def makeFile(pdbList):
    with open("./temp.txt", "w") as f:
        for pdb in pdbList:
            f.write(os.path.join(path, pdb[1:3], "pdb" + pdb + ".ent.gz") + "\n")

'''def makeTempFile():
    with tempfile.TemporaryFile() as fp:
        fp.write(b"Hi")
        print(tempfile.gettempprefix())'''


def main():
    proteins = openFile()
    parse(proteins)

if __name__ == "__main__":
    main()