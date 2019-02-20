import argparse

ap = argparse.ArgumentParser(description = "Removes everything except first ID from FASTA ID line")
ap.add_argument("-f", help = "FASTA file name", required = True)

opts = ap.parse_args()
FASTAFileName = opts.f

for Line in open(FASTAFileName):
    if ">" in Line:
        FullID = Line[1:-1].split(" ")[0]
        ShortID = FullID
        if "|" in ShortID:
            ShortID = ShortID.split("|")[1]

        print(">" + ShortID)
    else:
        print(Line[:-1])