import argparse
import os
import subprocess

import statistics
import linecache

HIT_Id = 0
HIT_Score = 1
HIT_Start = 2
HIT_Stop = 3
HIT_Sequence = 4
HIT_ClusterId = 5

HIT_IsInVicinity = 7
HIT_DistToBait = 10

StrictThreshold = 0.9


def GetEffectiveNoFromAlignment(TargetHits, PrefixFileName, FilterField = -1):
    ClusterCount = -1
    Dists = []

    FAFileName = PrefixFileName + ".tmp"
    ClustFileName = PrefixFileName + ".clust"

    ContentExists = False

    DistsDict = dict()
    DistField = HIT_DistToBait

    with open(FAFileName, "w+") as TempFile:
        for Line in TargetHits:
            if (FilterField == -1) or (FilterField != -1 and Line[FilterField] == "1"):
                TempFile.write(">gi|" + Line[HIT_Id] + "\n")
                TempFile.write(Line[HIT_Sequence] + "\n")
                ContentExists = True

                DistsDict[Line[HIT_Id]] = int(Line[DistField])

    if ContentExists:
        devnull = open(os.devnull, 'wb')
        subprocess.call(
            "bash RunClust.sh " + FAFileName + " " + str(StrictThreshold) + " " + ClustFileName, shell=True, stdout=devnull, stderr=devnull)

        ClusterCount = 0
        for Line in open(ClustFileName, "r"):
            ClusterCount += 1
            Dists.append(DistsDict[Line[:-1].split("\t")[0]])

    if os.path.exists(FAFileName):
        os.remove(FAFileName)
    if os.path.exists(ClustFileName):
        os.remove(ClustFileName)

    MedianDist = -1
    if len(Dists) > 0:
        MedianDist = statistics.median_low(Dists)

    return ClusterCount, MedianDist

def MakeGIListFile(ClusterHits, TempGIListFileName):
    GIsFileName = TempGIListFileName

    GIs = set([x[HIT_Id] for x in ClusterHits])
    with open(GIsFileName, mode = "w+") as GIsFile:
        for GI in GIs:
            GIsFile.write(GI + "\n")

def GetSelectedClusterIDs(ClustersFileName, LineNo):
    ClusterLine = linecache.getline(ClustersFileName, LineNo)
    ClusterLine = ClusterLine[:-1]

    ClusterLine = ClusterLine[ClusterLine.find("\t") + 1:]
    IDs = ClusterLine.split(" ")

    return IDs

def FilterBlastHits(ClusterId, ClusterHits, Database, ClustersFileName):
    ClusterGIListFileName = ClusterId + ".gis"
    ClusterFASTA = ClusterId + ".faa"
    ClusterSubclustersFileName = ClusterId + ".clust_permissive"

    MakeGIListFile(ClusterHits, ClusterGIListFileName)

    subprocess.call(
        "blastdbcmd -db " + Database + " -entry_batch " + ClusterGIListFileName + " -long_seqids > " + ClusterFASTA, shell = True)

    devnull = open(os.devnull, 'wb')
    subprocess.call(
        "bash RunClust.sh " + ClusterFASTA + " 0.3 " + ClusterSubclustersFileName, shell = True, stdout = devnull, stderr = devnull)

    ClusterLineNo = int(ClusterId.split("_")[1])
    OriginalClusterIDs = set(GetSelectedClusterIDs(ClustersFileName, ClusterLineNo))

    FinalClusterSet = set()

    for Line in open(ClusterSubclustersFileName, "r"):
        Clusters = set(Line[:-1].split("\t")[1].split(" "))
        if len(Clusters.intersection(OriginalClusterIDs)) > 0:
            FinalClusterSet = FinalClusterSet.union(Clusters)

    if os.path.exists(ClusterGIListFileName):
        os.remove(ClusterGIListFileName)
    if os.path.exists(ClusterFASTA):
        os.remove(ClusterFASTA)
    if os.path.exists(ClusterSubclustersFileName):
        os.remove(ClusterSubclustersFileName)

    return [x for x in ClusterHits if x[0] in FinalClusterSet]

def LoadHits(SortedHitsFileName):
    ClusterHits = []
    GIs = set()
    with open(SortedHitsFileName, "r") as SortedHits:
        for Line in SortedHits:
            LineValues = Line[:-1].split("\t")
            GIs.add(LineValues[HIT_Id])
            ClusterHits.append(LineValues)

    return ClusterHits


ap = argparse.ArgumentParser(description = "Calculates number of different proteins at the baits, in all genomic database and calculates median distance to the baits using sorted PSIBLAST search results")
ap.add_argument("-f", help = "Sorted PSIBLAST hits file name", required = True)
ap.add_argument("-o", help = "Result tsv file", required = True)
ap.add_argument("-d", help = "Genomic database", required = True)
ap.add_argument("-c", help = "Permissive clusters file name", required = True)

opts = ap.parse_args()
SortedHitsFileName = opts.f
ResultFileName = opts.o
Database = opts.d
ClustersFileName = opts.c

ClusterId = os.path.splitext(os.path.basename(SortedHitsFileName))[0]
ClusterHits = LoadHits(SortedHitsFileName)

# filter out all protein clusters that don't contain representatives of original cluster
ClusterHits = FilterBlastHits(ClusterId, ClusterHits, Database, ClustersFileName)

# get effective size and distance
AllClusterCount, Tmp = GetEffectiveNoFromAlignment(ClusterHits, SortedHitsFileName, FilterField= -1)
BaitClusterCount, DistToBait = GetEffectiveNoFromAlignment(ClusterHits, SortedHitsFileName, FilterField= HIT_IsInVicinity)

with open(ResultFileName, "w") as ResultFile:
    ResultFile.write(ClusterId + "\t" + str(BaitClusterCount) + "\t" + str(AllClusterCount) + "\t" +
                     str(DistToBait) + "\t" + str(BaitClusterCount / AllClusterCount) + "\n")


