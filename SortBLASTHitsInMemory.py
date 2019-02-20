import os
import argparse

HIT_Score = 1
HIT_Cluster = 5
HIT_ALIGNMENT_START = 2
HIT_ALIGNMENT_STOP = 3
HIT_ALIGNMENT_Sequence = 4
HIT_Protein_ID = 0

INFO_CONTIG = 0
INFO_BAITED = 1
INFO_ORF_START = 2
INFO_ORF_STOP = 3
INFO_ORF_DISTANCE_TO_BAIT = 4

def AddHitToDict(TargetHitsDict, Hit):
    BetterHitExists = False

    # check for better hit, remove lowest score hit if overlapping
    for DictHit in TargetHitsDict[Hit[HIT_Protein_ID]]:
        Intersect = min(Hit[HIT_ALIGNMENT_STOP], DictHit[HIT_ALIGNMENT_STOP]) - max(Hit[HIT_ALIGNMENT_START], DictHit[HIT_ALIGNMENT_START])
        MinLength = min(Hit[HIT_ALIGNMENT_STOP] - Hit[HIT_ALIGNMENT_START], DictHit[HIT_ALIGNMENT_STOP] - DictHit[HIT_ALIGNMENT_START])
        if Intersect/MinLength > MIN_COVERAGE_FILTER:
            if Hit[HIT_Score] > DictHit[HIT_Score]:
                TargetHitsDict[Hit[HIT_Protein_ID]].remove(DictHit)
            else:
                BetterHitExists = True
                break # better hit exists

    if not BetterHitExists:
        TargetHitsDict[Hit[HIT_Protein_ID]].append(Hit)

    return TargetHitsDict

def AddTargetHitsValues(TargetHitsFileName, TargetHitsDict, ProteinInfoDict):#TargetHitsDict):
    BLAST_FORMAT_ALIGNMENT_TARGET_PROTEIN_ID = 1
    BLAST_FORMAT_ALIGNMENT_START = 8
    BLAST_FORMAT_ALIGNMENT_STOP = 9
    BLAST_FORMAT_ALIGNMENT_SCORE = 10
    BLAST_FORMAT_ALIGNMENT_SEQUENCE = 7

    TargetValues = []
    BestScoreLine = []
    MaxScore = 0
    ClusterId = os.path.splitext(os.path.basename(TargetHitsFileName))[0]

    with open(TargetHitsFileName, "r") as TargetHitsFile:
        for Line in TargetHitsFile:
            if Line[0] == "#":
                continue

            LineValues = Line.split("\t")
            LineValues[-1] = LineValues[-1][:-1]

            #Alignment Data
            Start = int(LineValues[BLAST_FORMAT_ALIGNMENT_START])
            Stop = int(LineValues[BLAST_FORMAT_ALIGNMENT_STOP])
            Score = int(LineValues[BLAST_FORMAT_ALIGNMENT_SCORE])
            ProteinID = LineValues[BLAST_FORMAT_ALIGNMENT_TARGET_PROTEIN_ID].split("|")[1] # gi|ID| format expected

            HitLine = [ProteinID, Score, Start, Stop, LineValues[BLAST_FORMAT_ALIGNMENT_SEQUENCE], ClusterId]
            HitLine.extend(ProteinInfoDict[ProteinID])

            if Score > MaxScore:
                BestScoreLine = HitLine
                MaxScore = Score

            TargetValues.append(HitLine)

    if len(TargetValues) > 0: # no hits
        # max length of alignment is taken as selflength to filter by coverage
        # this value is used as representative length of alignment
        BestScoreLength = len(BestScoreLine[HIT_ALIGNMENT_Sequence].replace("-", ""))

        for Hit in TargetValues:
            if len(Hit[HIT_ALIGNMENT_Sequence].replace("-", "")) / BestScoreLength >= MIN_OVERLAP_FILTER: # filtering by coverage
                if Hit[HIT_Protein_ID] in TargetHitsDict:
                    TargetHitsDict = AddHitToDict(TargetHitsDict, Hit)
                else:
                    TargetHitsDict[Hit[HIT_Protein_ID]] = [Hit]

    return ClusterId

def GetDistanceToSeed(Start, Stop, Seed):
    if len(Seed) == 0:
        return 100000

    return min(abs(Start - Seed[1]), abs(Seed[0] - Stop))

def GetDistToSeedDict(AnnotationFileName, SeedsDict):
    LOCI_PROTEIN_ID = 0
    LOCI_CONTIG_ID = 4
    LOCI_COORDINATES = 1

    DistToSeed = dict()

    GIs = dict()
    IslandSeeds = []
    for Line in open(AnnotationFileName, "r"):
        LineValues = Line[:-1].split("\t")

        if LineValues[LOCI_PROTEIN_ID] == "===":
            if len(GIs) > 0:
                if len(IslandSeeds) == 0:
                    raise "NoSeedFound"
                else:
                    for GI in GIs:
                        DistToSeed[GI] = min([abs(x - GIs[GI]) for x in IslandSeeds])

            GIs = dict()
            IslandSeeds = []
        else:
            GIs[LineValues[LOCI_PROTEIN_ID]] = len(GIs) + 1

            StartStop = LineValues[LOCI_COORDINATES].split("..")
            if IsInSeeds(LineValues[LOCI_CONTIG_ID], int(StartStop[0]), int(StartStop[1]), SeedsDict):
                IslandSeeds.append(len(GIs))

    if len(GIs) > 0:
        if len(IslandSeeds) == 0:
            raise "NoSeedFound"
        else:
            for GI in GIs:
                DistToSeed[GI] = min([abs(x - GIs[GI]) for x in IslandSeeds])

    return DistToSeed

def LoadProteinInfoDict(PTYFileName, VicinityProteinIDs, LociDists):
    PTY_PROTEIN_ID = 6
    PTY_CONTIG_ID = 4
    PTY_COORDINATES = 1

    ProteinInfoDict = dict()

    for Line in open(PTYFileName):
        if Line[0] == "#":
            continue
        LineValues = Line[:-1].split("\t")

        StartStop = LineValues[PTY_COORDINATES].split("..")

        if LineValues[PTY_PROTEIN_ID] in VicinityProteinIDs:
            InVicinity = "1"
        else:
            InVicinity = "0"

        if LineValues[PTY_PROTEIN_ID] in LociDists:
            Dist = LociDists[LineValues[PTY_PROTEIN_ID]]
        else:
            Dist = 10000

        # in assumption that Protein IDs are unique
        ProteinInfoDict[LineValues[PTY_PROTEIN_ID]] = [LineValues[PTY_CONTIG_ID], InVicinity,
                                                       StartStop[0], StartStop[1], Dist]

    return ProteinInfoDict

def LoadList(FileName):
    ListValues = []
    for Line in open(FileName):
        ListValues.append(Line[:-1])
    return ListValues

def LoadSeedsDict(SeedsFileName):
    Seeds = dict()
    ContigField = 1
    StartField = 2
    EndField = 3
    Offset = 0

    return AddSeedsDictByFields(SeedsFileName, Seeds, ContigField, StartField, EndField, Offset)

def AddSeedsDictByFields(FileName, SeedDict, ContigField, StartField, EndField, Offset):
    with open(FileName, "r") as File:
        for Line in File:
            LineValues = Line.split("\t")
            LineValues[-1] = LineValues[-1][:-1]
            ID = LineValues[0]

            Start, End = minmax(int(LineValues[StartField]), int(LineValues[EndField]))
            if LineValues[ContigField] in SeedDict:
                SeedDict[LineValues[ContigField]].append([max(Start - Offset, 0), End + Offset, ID, LineValues])
            else:
                SeedDict[LineValues[ContigField]] = [[max(Start - Offset, 0), End + Offset, ID, LineValues]]

    for Key in SeedDict:
        SeedDict[Key] = sorted(SeedDict[Key], key = lambda e: e[0])

    return SeedDict

def minmax(X, Y):
    return min(X, Y), max(X, Y)

def IsInSeeds(ContigID, Start, Stop, SeedsDict):
    Seed = GetSeed(ContigID, Start, Stop, SeedsDict)
    if len(Seed) == 0:
        return False
    return True

def GetSeed(ContigID, Start, Stop, SeedsDict):
    if not ContigID in SeedsDict:
        return []

    for ORF in SeedsDict[ContigID]:
        if (Start <= ORF[1]) and (Stop >= ORF[0]):
            return ORF

    return []


ap = argparse.ArgumentParser(description = "Collecting all PSIBLAST hits in specified folder and sorting them between clusters. Marking hits that are nearby seeds and storing distance to the seeds.")
ap.add_argument("-c", help = "Folder name where hits stored", required = True)
ap.add_argument("-o", help = "Folder name where sorted result hits stored", required = True)
ap.add_argument("-p", help = "PTY file containing coordinates of ORFs of genomic database", required = True)
ap.add_argument("-i", help = "List of protein IDs in vicinity of baits", required = True)
ap.add_argument("-s", help = "Baits coordinates", required = True)
ap.add_argument("-v", help = "Vicinity around the baits, needed to calculate distances to the baits", required = True)
ap.add_argument("-z", help = "Overlap threshold, hits are subject to sort between two profiles if they are overlapping for more than threshold value", required = True)
ap.add_argument("-x", help = "Coverage threshold, hits are stored if they cover original profile for more than threshold value", required = True)

opts = ap.parse_args()
ClustersHitsFolder = opts.c
SortedHitsFolder = opts.o
PTYFileName = opts.p
VicinityProteinIDsFileName = opts.i
SeedsFileName = opts.s
VicinityFileName = opts.v
MIN_COVERAGE_FILTER = float(opts.x)
MIN_OVERLAP_FILTER = float(opts.z)


VicinityProteinIDs = set(LoadList(VicinityProteinIDsFileName))
Seeds = LoadSeedsDict(SeedsFileName)
LociDists = GetDistToSeedDict(VicinityFileName, Seeds)
ProteinInfoDict = LoadProteinInfoDict(PTYFileName, VicinityProteinIDs, LociDists)

if not os.path.exists(SortedHitsFolder):
    os.makedirs(SortedHitsFolder)

ProteinHitsDict = dict()
Clusters = []
for FileName in os.listdir(ClustersHitsFolder):
    if not FileName.endswith(".hits"):
        continue

    ClusterId = AddTargetHitsValues(ClustersHitsFolder + FileName, ProteinHitsDict, ProteinInfoDict)
    Clusters.append(ClusterId)

print("Loading BLAST hits complete.")

ClusterHitsDict = dict()
for ProteinID in ProteinHitsDict:
    for Hit in ProteinHitsDict[ProteinID]:
        if Hit[HIT_Cluster] in ClusterHitsDict:
            ClusterHitsDict[Hit[HIT_Cluster]].append(Hit)
        else:
            ClusterHitsDict[Hit[HIT_Cluster]] = [Hit]

print("Moving sorting hits by clusters complete.")

for ClusterId in Clusters:
    if ClusterId in ClusterHitsDict:
        with open(SortedHitsFolder + str(ClusterId) + ".hits_sorted", "w") as SortedFile:
            for Hit in ClusterHitsDict[ClusterId]:
                ResLine = "\t".join([str(x) for x in Hit])

                SortedFile.write(ResLine + "\n")
    else: # empty file
        with open(SortedHitsFolder + str(ClusterId) + ".hits_sorted", "w") as SortedFile:
            SortedFile.write("")
