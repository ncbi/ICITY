import argparse
import subprocess
import os

ap = argparse.ArgumentParser(description = "Creates protein profiles in specified folder for clusters in given file")
ap.add_argument("-f", help = "Clusters file name", required = True)
ap.add_argument("-c", help = "Folder name where profiles will be saved", required = True)
ap.add_argument("-d", help = "Path to protein database", required = True)

opts = ap.parse_args()
ClustersFileName = opts.f
ClustersFolderName = opts.c
Database = opts.d

ClusterNo = 0
TmpIDsFileName = "Tmp_IDs.lst"
TmpFASTAFileName = "Tmp_FASTA.faa"

if not os.path.exists(ClustersFolderName):
    subprocess.call("mkdir " + ClustersFolderName, shell = True)

for Line in open(ClustersFileName):
    ClusterNo += 1
    ClusterIDs = Line[:-1].split("\t")[1].split(" ")
    ClusterProfileFileName = "CLUSTER_" + str(ClusterNo) + ".ali"

    with open(TmpIDsFileName, "w") as IDsFile:
        IDsFile.write("\n".join(ClusterIDs))

    subprocess.call("blastdbcmd -db " + Database + " -entry_batch " + TmpIDsFileName + " -long_seqids > " + TmpFASTAFileName,
                    shell = True)

    subprocess.call("muscle -in " + TmpFASTAFileName + " -out " + ClustersFolderName + "/" + ClusterProfileFileName,
                    shell = True)

subprocess.call("rm " + TmpIDsFileName, shell = True)
subprocess.call("rm " + TmpFASTAFileName, shell = True)