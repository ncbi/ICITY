#!/bin/bash

if [ "$1" == "" ]
then
    echo "Folder name with sorted clusters as first parameter needed"
	exit 1
fi

FolderName="$1"

if [ "$2" == "" ]
then
    echo "Path to protein database as second parameter needed"
	exit 1
fi

DBPath="$2"


if [ "$3" == "" ]
then
    echo "Linear clusters file name is needed as second parameter"
	exit 1
fi

ClustersFileName="$3"


for f in $FolderName/*.hits_sorted;           
do
	echo $f
	python GetIcityForBLASTHits.py -f $f -o $FolderName/$(basename "$f" .hits_sorted).tsv -d $DBPath -c $ClustersFileName	
done