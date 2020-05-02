'''
Will remove duplicates and similar sequences in the fasta file
Duplicates are based on protein sequence to reduce redundancy during alignment
Sequence similarity is based on Levenshtein Edit Distance between strings
'''

import Levenshtein
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Gets Levenshtein distance of the sequence and all values in the same level
# Returns False if >= 90% sequence similarity
def calcDist(seq, dict):
    for sequence in dict:
        if Levenshtein.ratio(seq, sequence) >= 0.9:
            return False
    return True

# Removes duplicates in the fasta file
def remove(fasta):
    # Holds records that will be added to fasta file
    records = []
    levelDict = {}
    removed = set()

    for reads in SeqIO.parse(fasta, "fasta"):
        seq, level = str(reads.seq), str(reads.id).split("|")[2]
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
                records.append(SeqRecord(Seq(seq, IUPAC.protein), id=reads.id))
    # Rewrites the fasta file so that there are no duplicates or similar sequences (< 3 differences)
    with open(fasta, 'w') as output:
        SeqIO.write(records, output, 'fasta')

    return records
