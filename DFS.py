'''
Goes through the json file to retrieve file paths.
These paths will be used to merge files from the child directories to the parent directory.
Ex: 1/1/1/4, 1/1/1/5, 1/1/1/10 fasta files will be merged into 1/1/1
Ex: 1/1/1, 1/1/2, 1/1/3 fasta files will be merged into 1/1
'''

import os
import json
import shutil

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json


# Merges the fasta files of the subdirectory into the parent directory
# paths: list of file paths to be merged
# parent: the parent path that the files will be merged into
def merge(paths, parent):
    with open(parent, 'wb') as outfile:
        for file in paths:
            with open(file, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

# Go through the dictionary and create lists of file paths to be merged
# __path is the list of files to be merged
# __pathP is the parent path to merge the file into
def parseDict(proteins):
    for x, x_dict in proteins.items():
        Hpath = []
        for h, h_dict in x_dict.items():
            Tpath = []
            Hpath.append(os.path.join(x, h,"sequences", "sequences.fasta"))
            HpathP = os.path.join(x, "sequences", "sequences.fasta")
            for t, t_dict in h_dict.items():
                Fpath = []
                Tpath.append(os.path.join(x, h, t, "sequences", "sequences.fasta"))
                TpathP = os.path.join(x, h, "sequences", "sequences.fasta")
                for f, f_dict in t_dict.items():
                    Fpath.append(os.path.join(x, h, t, f, "sequences", "sequences.fasta"))
                    FpathP = os.path.join(x, h, t, "sequences", "sequences.fasta")
                merge(Fpath, FpathP)
            merge(Tpath, TpathP)
        merge(Hpath, HpathP)

def main():
    proteins = openFile()
    parseDict(proteins)


if __name__ == "__main__":
    main()
