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

time python3 directories.py "$1"
time python3 pdb.py "$1"
