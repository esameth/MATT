'''
Will remove duplicates and similar sequences in the fasta file
Duplicates are based on protein sequence to reduce redundancy during alignment
Sequence similarity is based on Levenshtein Edit Distance between strings
'''

import re
import Levenshtein
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from collections import defaultdict

# Nested dictionary to hold our levels and the proteins in it
proteins = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

# Gets Levenshtein distance of the sequence and all values in the same level
# Returns False if >= 90% sequence similarity
def calcDist(seq, dict):
    for sequence in dict:
        if Levenshtein.ratio(seq, sequence) >= 0.9:
            return False
    return True

# Removes duplicates in the fasta file, and creates the dictionaries necessary for us to create directories
# and to determine which proteins will be aligned together
def parse(fasta):
    # Holds records that will be added to fasta file
    records = []
    levelDict = {}
    removed = set()

    for reads in SeqIO.parse(fasta, "fasta"):
        id = str(reads.id).split("|")
        seq, pdb, level, chain = str(reads.seq), id[1][1:5], id[2], re.split(':|,', id[3])

        # Split chain if more than one and get only chain names
        chain = ','.join([chain[i] for i in range(len(chain)) if i % 2 == 0])

        # Get the levels
        f_id = level.split(".")
        x, h, t = f_id[0], f_id[1], f_id[2]
        # Some f_id do not have a F group
        f = f_id[3] if len(f_id) == 4 else "None"

        if level not in levelDict:
            levelDict[level] = []

        # If it has not been looked at already and it has not been removed, get the distance
        if seq not in levelDict[level] and seq not in removed:
            # Add it to the list of removed if there's more than 90% sequence identity
            if not calcDist(seq, levelDict[level]):
                removed.add(seq)
            # Enough similarities so add it to our dict and fasta file
            else:
                levelDict[level].append(seq)
                records.append(SeqRecord(Seq(seq, IUPAC.protein), id=reads.id, name=reads.id, description=reads.id))
                # Add each protein's information to the dict
                proteins[x][h][t][f].append((pdb, chain))

    # Rewrites the fasta file so that there are no duplicates or similar sequences (< 3 differences)
    with open(fasta, 'w') as output:
        SeqIO.write(records, output, 'fasta')

    return levelDict.keys(), proteins
