ICITY_CONFIG_INPUT = {
    "PTYFile": "Database/CDS.pty",
    "SeedsFile": "Seeds.tsv",
    "NeighborhoodVicinitySize": 10000,
    "PathToDatabase": "Database/ProteinsDB",
    "PermissiveClusteringThreshold": 0.3,
    "SortingOverlapThreshold": 0.4,
    "SortingCoverageThresold": 0.25
}

ICITY_CONFIG_OUTPUT = {
    "ICITYFileName": "Relevance.tsv",
    "VicinityClustersFileName": "VicinityPermissiveClustsLinear.tsv"
}

ICITY_CONFIG_TEMPORARYFILES = {
    "VicinityFileName": "Vicinity.tsv",
    "VicinityIDsFileName": "VicinityIDs.lst",
    "VicinityFASTAFileName": "Vicinity.faa",
    "VicinityClustersFileName": "VicinityPermissiveClustsLinear.tsv",
    "ProfilesFolder": "CLUSTERS/",
    "SortedBLASTHitsFolder": "CLUSTERS/Sorted/",
}