#!/bin/bash
PROGRAM="./FastCommunityMH_w"
METRIC="modularity"

FLAG_FILE="-f"
FLAG_LABEL="-l"

EXT_IN=".wpairs"

PATH="../../input/"
GRAPH_NAME=(
	"soc-Slashdot0902"
)

for name in "${GRAPH_NAME[@]}"
do
	RUN_COMMAND="$PROGRAM $FLAG_FILE ${PATH}${name}${EXT_IN} $FLAG_LABEL $METRIC"
	echo $name
	$RUN_COMMAND
done