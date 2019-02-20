import argparse

ap = argparse.ArgumentParser(description = "Convert mmseqs2 tsv clustering output to format where all cluster members represented by one line")
ap.add_argument("-f", help = "Cluster file name", required = True)

opts = ap.parse_args()
ClusterFileName = opts.f

RepresentativesToMembersDict = dict()
for Line in open(ClusterFileName):
    LineValues = Line[:-1].split("\t")

    LineValues[0] = LineValues[0].replace(">", "")
    LineValues[1] = LineValues[1].replace(">", "")

    if LineValues[0] not in RepresentativesToMembersDict:
        RepresentativesToMembersDict[LineValues[0]] = []

    RepresentativesToMembersDict[LineValues[0]].append(LineValues[1])

for Repr in RepresentativesToMembersDict:
    print(Repr + "\t" + " ".join(RepresentativesToMembersDict[Repr]))