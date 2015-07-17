#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./FastCommunityMH"
METRIC="coverage"		  	   # change

# extensions
EXT_IN=".pairs"

# input file info
LOC="../../input/"			   # change
GRAPH_NAME=(				   # change
	"karate"
	"dolphins"
	"polbooks"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_INPUT="${LOC}${name}${EXT_IN}"
	$PROGRAM -f $FILE_INPUT -l $METRIC
	# ./FastCommunityMH -f ../../input/{name}.pairs -l {METRIC}
done