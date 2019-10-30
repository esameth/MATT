import os

'''def walk(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for file in files:
            path = os.path.join(root.rsplit('/', 2) [0], "sequences", "multifasta.txt")
            if file == "multifasta.txt":
                fasta = os.path.join(root,file)
            print(fasta, "----", path)

            #create(fasta, path)


def merge(arr):
    for path in arr:
        with open(path, 'r') as infile:
            for line in infile:
                outfile.write(line)

    with open(fasta, 'r') as infile, open(path, "a+") as outfile:
        for line in infile:
            outfile.write(line)
    "return path to parent directory"'''


def dfs(path):
    for f in os.scandir(path):
        if f.is_dir():
            dfs(f)
            print(f.path)


def main():
    dfs("1")


if __name__ == "__main__":
    main()
