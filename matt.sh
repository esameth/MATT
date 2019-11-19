#!/bin/bash

python3 parser.py
echo Finished creating dictionary
python3 fasta.py
echo Finished removing duplicates
python3 directories.py
echo Finished creating directories
python3 DFS.py
echo Finished merging fasta files
python3 pdb.py
echo Finished alignment
