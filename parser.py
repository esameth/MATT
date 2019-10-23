'''
Parse through the ECOD F40 Representative Domain text file and return a json file that includes
each protein's pdb ID, chain, and hierarchy
'''

import json
import re
from itertools import islice

# Each protein has a pdb id and a chain
class Protein:
    def __init__(self, u_id, pdb, chain):
        self.u_id = u_id
        self.pdb = pdb
        self.chain = chain

def data():
    proteins = {}

    with open('../../Downloads/ecod.latest.F40.domains.copy.txt') as f:
        # Skip the first 5 lines because they are headers
        for line in islice(f, 5, None):
            line = (line.replace('"', '')).split("\t")

            # Get the values in specific columns
            u_id, f_id, pdb, chain = line[0], line[3], line[4], line[5]

            # Some proteins have multiple chains designated as a "." that must be read from the pdb range column
            if chain == '.':
                chain = line[6]
                chain = re.split(':|,', chain)
                chain = ','.join([chain[i] for i in range(len(chain)) if i%2 == 0])

            x_level, h_level, t_level = f_id.split(".")[0], f_id.split(".")[1], f_id.split(".")[2]

            # Some f_id do not have a F group
            f_level = f_id.split(".")[3] if len(f_id.split(".")) == 4 else "None"

            # Create X-level (is not dependent on other levels)
            if x_level not in proteins:
                proteins[x_level] = {}

            # Create H level which is dependent on X level
            if h_level not in proteins[x_level]:
                proteins[x_level][h_level] = {}

            # Create T level which is dependent on X level and H level
            if t_level not in proteins[x_level][h_level]:
                proteins[x_level][h_level][t_level] = {}

            # Create F level which is dependent on X level, H level, and T level
            if f_level not in proteins[x_level][h_level][t_level]:
                proteins[x_level][h_level][t_level][f_level] = []

            # Add each protein's information to the dict
            proteins[x_level][h_level][t_level][f_level].append(Protein(u_id, pdb, chain).__dict__)

    return proteins


def main():
    proteins = data()

    # Creates the json file for easy opening of the data
    with open('proteins.json', 'w+', encoding='utf-8') as f:
        json.dump(proteins, f, ensure_ascii=False, indent=4)


if __name__ == "__main__":
    main()
