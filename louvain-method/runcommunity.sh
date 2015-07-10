#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./community"
ALG="-lm"
METRIC="-coverage"			   # change
LEVEL="-1"

# extensions
EXT_IN=".bin"
EXT_OUT1=".info"
EXT_OUT2=".tree"

# input file info
LOC="../../input/"			   # change
GRAPH_NAME=(				   # change
	"adjnoun"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_INPUT="${LOC}${name}${EXT_IN}"
	FILE_OUTPUT1="${LOC}${name}${ALG}${METRIC}${EXT_OUT1}"
	FILE_OUTPUT2="${LOC}${name}${ALG}${METRIC}${EXT_OUT2}"
	${PROGRAM}${METRIC} $FILE_INPUT -l $LEVEL -v $FILE_OUTPUT1 > $FILE_OUTPUT2
	# ./community{METRIC} ../../input/{name}.bin -l -1 -v ../../input/{name}-lm{METRIC}.info > ../../input/{name}-lm{METRIC}.tree
done