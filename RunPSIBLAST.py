import os
import subprocess
import argparse

def GetNumberOfSequences(FASTAFileName):
    Count = 0
    for Line in open(FASTAFileName):
        if ">" in Line:
            Count += 1

    return Count

ap = argparse.ArgumentParser(description = "Running PSIBLAST for each protein profile in specified folder and saving results in same folder")
ap.add_argument("-c", help = "Folder name where profiles stored", required = True)
ap.add_argument("-d", help = "Path to protein database", required = True)

opts = ap.parse_args()
ClustersFolderName = opts.c
Database = opts.d

Evalue = 0.000001

for FileName in os.listdir(ClustersFolderName):
    if not FileName.endswith("ali"):
        continue

    BLASTHitsFileName = ClustersFolderName + "/" + os.path.splitext(FileName)[0] + ".hits"

    if GetNumberOfSequences(ClustersFolderName + "/" + FileName) == 1:
        Query = "-query"
    else:
        Query = "-in_msa"
    Query += " " + ClustersFolderName + "/" + FileName

    # to parallelize
    print("Processing " + FileName)
    subprocess.call(
        "psiblast -db " + Database + " " +
        "-outfmt \"7 qseqid sseqid slen sstart send evalue qseq sseq qstart qend score\"" +
        " -seg no -evalue " + str(Evalue) +
        " -dbsize 20000000 -max_target_seqs 10000 -comp_based_stats no " +
        Query + " > " + BLASTHitsFileName, shell = True)
