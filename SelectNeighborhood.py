import argparse

def WriteSeedIslands(SeedIslandPTYLines, ResultFile):
    if len(SeedIslandPTYLines) > 0:
        ResultFile.write("===\n")
        for LineValues in SeedIslandPTYLines:
            ResultFile.write(LineValues[0] + "\t" + LineValues[1] + "\t" + LineValues[2] + "\t" +
                             LineValues[3] + "\t" + LineValues[4] + "\n")

def GetSeedIslandCoordinates(SeedCoordinates, SeedIslandOffset):
    if not SeedCoordinates:
        SeedIslandStart = -1
        SeedIslandEnd = -1
    else:
        IslandFound = False
        SeedIslandStart = int(SeedCoordinates[0][0]) - SeedIslandOffset

        while not IslandFound:
            SeedIslandEnd = int(SeedCoordinates[0][1]) + SeedIslandOffset
            del SeedCoordinates[0]

            if len(SeedCoordinates) == 0:
                IslandFound = True
            elif int(SeedCoordinates[0][0]) - SeedIslandOffset > SeedIslandEnd:
                IslandFound = True

    return SeedIslandStart, SeedIslandEnd

def WriteIslands(ContigLines, ResultFile, ContigID, SeedsDict, SeedIslandOffset):
    if ContigID not in SeedsDict:
        return

    if len(SeedsDict[ContigID]) == 0:
        return

    ContigLines = sorted(ContigLines, key = lambda e: e[0])

    SeedIslandPTYLines = []
    for Line in ContigLines:
        LineValues = Line[1]

        Coordinates = LineValues[PTY_Coordinates].split("..")
        Start = int(Coordinates[0])
        End = int(Coordinates[1])

        if IsInSeeds(ContigID, Start - SeedIslandOffset, End + SeedIslandOffset, SeedsDict):
            SeedIslandPTYLines.append(
                [LineValues[PTY_GI], LineValues[PTY_Coordinates], LineValues[PTY_Strand],
                 LineValues[PTY_SpecieName], LineValues[PTY_ContigId]])
        else:
            if len(SeedIslandPTYLines) > 0:
                WriteSeedIslands(SeedIslandPTYLines, ResultFile)

            SeedIslandPTYLines = []

    if len(SeedIslandPTYLines) > 0:
        WriteSeedIslands(SeedIslandPTYLines, ResultFile)

def LoadSeedsDict(SeedsFileName, NoID = False):
    Seeds = dict()
    ContigField = 1
    StartField = 2
    EndField = 3
    Offset = 0

    return AddSeedsDictByFields(SeedsFileName, Seeds, ContigField, StartField, EndField, Offset, NoID)

def AddSeedsDictByFields(FileName, SeedDict, ContigField, StartField, EndField, Offset, NoID = False):
    with open(FileName, "r") as File:
        for Line in File:
            LineValues = Line.split("\t")
            LineValues[-1] = LineValues[-1][:-1]
            if not NoID:
                ID = LineValues[0]
            else:
                ID = "Seed"

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

# PTY tsv file columns
PTY_LocusTag = 0
PTY_Coordinates = 1
PTY_Strand = 2
PTY_SpecieName = 3
PTY_ContigId = 4
PTY_AccessionNo = 5
PTY_GI = 6

# Reading arguments from command line
ap = argparse.ArgumentParser(description = "Selects loci around seeds from pty file.")
ap.add_argument("-p", help = "PTYDataFileName, complete pty for contigs", required = True)
ap.add_argument("-s", help = "SeedsFileName, seeds tsv file", required = True)
ap.add_argument("-o", help = "ResultFileName, output pty file", required = True)
ap.add_argument("-d", help = "Offset around seed (base pairs)", required = True)
opts = ap.parse_args()

PTYDataFileName = opts.p
SeedsFileName = opts.s
ResultFileName = opts.o
SeedIslandOffset = int(opts.d)

# Loading seeds file
SeedsDict = LoadSeedsDict(SeedsFileName)

# Run through PTY and cut loci
ContigId = ""
ContigLines = []
with open(ResultFileName, "w") as ResultFile:
    with open(PTYDataFileName, "r") as PTYData:
        for Line in PTYData:
            Line = Line[:-1]
            LineValues = Line.split("\t")

            if LineValues[0] == "#":
                continue

            if ContigId != LineValues[PTY_ContigId]: # write island on change of a contig
                WriteIslands(ContigLines, ResultFile, ContigId, SeedsDict, SeedIslandOffset)
                ContigLines = []

            ContigId = LineValues[PTY_ContigId]
            Coordinates = LineValues[PTY_Coordinates].split("..")
            Start = int(Coordinates[0])

            ContigLines.append([Start, LineValues])

        WriteIslands(ContigLines, ResultFile, ContigId, SeedsDict, SeedIslandOffset)