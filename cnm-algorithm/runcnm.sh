#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./FastCommunityMH"
METRIC="modularity"		  	   # change

# flags - see description of alg for additional flags
FLAG_FILE="-f"
FLAG_LABEL="-l"

# extensions
EXT_IN=".pairs"

# input file info
PATH="../../input/"
GRAPH_NAME=(				   # change
	"soc-pokec-relationships"
)

for name in "${GRAPH_NAME[@]}"
do
	RUN_COMMAND="$PROGRAM $FLAG_FILE ${PATH}${name}${EXT_IN} $FLAG_LABEL $METRIC"
	echo $name
	$RUN_COMMAND
	# ./FastCommunityMH -f ../../input/{name}.pairs -l {METRIC}
done