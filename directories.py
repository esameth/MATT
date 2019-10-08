'''
Goes through the hierarchy of each protein and will create a multifasta file of all proteins within the level
'''

import json
import os

# Open the json file
def openFile():
    with open("proteins.json") as f:
        loaded_json = json.load(f)
    return loaded_json

# Go through the dictionary to create the directories based on hierarchy level name
def getLevel(proteins):
    for arch, arch_value in proteins.items():
        for x, x_value in arch_value.items():
            for h, h_value in x_value[1].items():
                for t, t_value in h_value[1].items():
                    for f, f_value in t_value[1].items():
                        if not os.path.exists(os.path.join(arch, x_value[0], h_value[0], t_value[0], f_value[0])):
                            os.makedirs(os.path.join(arch, x_value[0], h_value[0], t_value[0], f_value[0]))

def main():
    proteins = openFile()
    getLevel(proteins)


if __name__ == "__main__":
    main()
