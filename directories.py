'''
Parse through the levels dictionary and creates directories and the fasta file for it
'''

import sys
import os

fastaFile = sys.argv[1]

import parser
levelDict, proteins = parser.parse(fastaFile)

# Create the directories with the fasta files
def makeDir():
    for dir in levelDict:
        dir = dir.split(".")

        # Holds the folder path
        path = ""
        # Holds the ID to grep for the fasta file
        f_id = "|"

        if len(dir) != 4:
            dir.append("None")
        # Go through every single level and create a directory for it
        for i in range(len(dir)):
            level = dir[i]
            path = os.path.join(path, level)
            f_id += level
            f_id += "." if i < len(dir) - 1 else "|"

            if not os.path.exists(path):
                os.makedirs(os.path.join(path, 'sequences'))
                os.makedirs(os.path.join(path, 'alignments'))
                # Uses the grep command to pull from the fasta file those that belong to the hierarchy
                # Pipes it to a multifasta file in the sequences directory
                cmd = "LC_ALL=C fgrep '" + f_id + "' -A 1 " + fastaFile + " | grep -v '^--' >> " + os.path.join(path, 'sequences', 'multifasta.fasta')
                os.system(cmd)

if __name__ == "__main__":
    makeDir()
