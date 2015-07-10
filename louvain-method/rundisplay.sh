#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./display"
ALG="-lm"
METRIC="-coverage"			   # change

# extensions
EXT_IN=".tree"
EXT_OUT=".graph"

# input file info
LOC="../../input/"			   # change
GRAPH_NAME=(				   # change
	"adjnoun"
	"football"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_INPUT="${LOC}${name}${ALG}${METRIC}${EXT_IN}"
	FILE_OUTPUT="${LOC}${name}${ALG}${METRIC}${EXT_OUT}"
	$PROGRAM $FILE_INPUT > $FILE_OUTPUT
	# ./display ../../input/{name}-lm{METRIC}.tree > ../../input/{name}-lm{METRIC}.graph
done