#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./convert"

# extensions
EXT_IN=".pairs"
EXT_OUT=".bin"

# input file info
LOC="../../input/"			   # change
GRAPH_NAME=(				   # change
	"adjnoun"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_INPUT="${LOC}${name}${EXT_IN}"
	FILE_OUTPUT="${LOC}${name}${EXT_OUT}"
	$PROGRAM -i $FILE_INPUT -o $FILE_OUTPUT
	# ./convert -i ../../input/{name}.pairs -o ../../input/{name}.bin
done