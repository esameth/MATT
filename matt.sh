#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: [ECOD fasta]"
    exit
else
    if [ ! -f "$1" ]; then
        echo "'$1' does not exist"
        echo "Usage: [ECOD fasta]"
        exit
    fi
fi

time python3 parser.py "$1"
echo "Removed fasta duplicates and created directories"
time python3 pdb.py
echo "Finished all alignments"
