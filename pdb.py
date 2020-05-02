'''
Goes through the json file and parses the pdb id to get which folder the pdb file would lie in.
Creates a temp file of the pdb file paths that will be used as input to MATT to create alignments for each level
The temp file is deleted when done
'''

import json
import os
import tempfile
import subprocess

# Path to pdb files
pdbPath = "/data/pdb"
home = os.getcwd()

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

# Parse through the dictionary and pull the pdb id
# Each _pdbList will hold the pdb ids that will be in the alignment
# Each _path will have the path so we can change to that directory for outputting the alignment
def parse(proteins):
    for x, x_dict in proteins.items():
        XpdbList = []
        Xpath = os.path.join(x, "alignments")
        for h, h_dict in x_dict.items():
            HpdbList = []
            Hpath = os.path.join(x, h, "alignments")
            for t, t_dict in h_dict.items():
                TpdbList = []
                Tpath = os.path.join(x, h, t, "alignments")
                for f, f_dict in t_dict.items():
                    FpdbList = []
                    Fpath = os.path.join(x, h, t, f, "alignments")
                    for protein in f_dict:
                        # Add the pdb ids to each list
                        FpdbList.append((protein["pdb"], protein["chain"]))
                        TpdbList.append((protein["pdb"], protein["chain"]))
                        HpdbList.append((protein["pdb"], protein["chain"]))
                        XpdbList.append((protein["pdb"], protein["chain"]))

                    # Create the temp file of pdb paths for aligning
                    makeTempFile(FpdbList, Fpath)
                makeTempFile(TpdbList, Tpath)
            makeTempFile(HpdbList, Hpath)
        makeTempFile(XpdbList, Xpath)

# Create a temporary file that will be passed into MATT
def makeTempFile(pdbList, pathway):
    # Used so we can see the folder we are writing to
    fpath = pathway.replace("/", "_")
    # Change to the directory that will hold the alignments
    pathway = os.path.join(home, pathway)
    os.chdir(pathway)

    # If alignment folder is empty, run the alignment
    if len(os.listdir(pathway)) == 0:
        # Create a tem file that will hold the pdb paths
        with tempfile.NamedTemporaryFile(mode = "w+t", prefix = fpath) as temp:
            for pdb, chain in pdbList:
                temp.write(os.path.join(pdbPath, pdb[1:3], "pdb" + pdb + ".ent.gz") + ":" + chain + "\n")
            # Run MATT with the output names as 'alignment' and the list file as temp file name
            temp.seek(0)
            cmd = ["/usr/local/bin/matt", "-o", "alignment", "-t", "28", "-L", temp.name]
            proc = subprocess.Popen(cmd)
            # Timeout after an hour - print in realtime
            try:
                print(proc.communicate(timeout=3600))
            # If timed out, kill the subprocess
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.communicate()
            # If killed, output folder name to a file to later run
            if proc.returncode != 0:
                with open(os.path.join(home, "timedout.txt"), "a+") as f:
                    f.write(fpath + "\n")

    # Change back to the home directory
    os.chdir(home)

if __name__ == "__main__":
    proteins = openFile()
    parse(proteins)
