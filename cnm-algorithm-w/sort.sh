#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
ALG="-cnm"
METRIC="-modularity"		   # change

# extensions
EXT_IN=".groups"
EXT_OUT=".graph"

# input file info
LOC="../../cnm/"			   # change
GRAPH_NAME=(				   # change
	"adjnoun"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_INPUT="${LOC}${name}${ALG}${METRIC}${EXT_IN}"
	FILE_OUTPUT="${LOC}${name}${ALG}${METRIC}${EXT_OUT}"
	sort -n $FILE_INPUT -o $FILE_OUTPUT
done