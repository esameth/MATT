'''
Parse through the ECOD F40 Representative Domain text file and return a json file that includes
each protein's pdb ID, chain, and hierarchy
'''

import json
from itertools import islice

# Each protein has a pdb id and a chain
class Protein:
    def __init__(self, pdb, chain):
        self.pdb = pdb
        self.chain = chain

def data():
    proteins = {}

    with open('../../Downloads/ecod.latest.F40.domains.copy.txt') as f:
        # Skip the first 5 lines because they are headers
        for line in islice(f, 5, None):
            line = (line.replace('"', '')).split("\t")

            # Get the values in specific columns
            f_id, pdb, chain = line[3], line[4], line[5]
            arch_name, x_name, h_name, t_name, f_name = line[8], line[9], line[10], line[11], line[12]
            x_level, h_level, t_level = f_id.split(".")[0], f_id.split(".")[1], f_id.split(".")[2]

            # Some f_id do not have a F group
            if len(f_id.split(".")) == 4:
                f_level = f_id.split(".")[3]
            else:
                f_level = "None"

            if arch_name not in proteins:
                proteins[arch_name] = {}

            # Each key has 2 values: the name of the family and a dictionary of the other levels
            # Create X-level (is not dependent on other levels)
            if x_level not in proteins[arch_name]:
                proteins[arch_name][x_level] = [x_name, {}]

            # Create H level which is dependent on X level
            if h_level not in proteins[arch_name][x_level][1]:
                proteins[arch_name][x_level][1][h_level] = [h_name, {}]

            # Create T level which is dependent on X level and H level
            if t_level not in proteins[arch_name][x_level][1][h_level][1]:
                proteins[arch_name][x_level][1][h_level][1][t_level] = [t_name, {}]

            # Create F level which is dependent on X level, H level, and T level
            if f_level not in proteins[arch_name][x_level][1][h_level][1][t_level][1]:
                proteins[arch_name][x_level][1][h_level][1][t_level][1][f_level] = [f_name, []]

            # Add each protein's information to the dict
            proteins[arch_name][x_level][1][h_level][1][t_level][1][f_level][1].append(Protein(pdb, chain).__dict__)

    return proteins


def main():
    proteins = data()

    # Creates the json file for easy opening of the data
    with open('proteins.json', 'w+', encoding='utf-8') as f:
        json.dump(proteins, f, ensure_ascii=False, indent=4)


if __name__ == "__main__":
    main()
