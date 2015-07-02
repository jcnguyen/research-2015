#!/bin/bash
PROGRAM="./FastCommunityMH"
METRIC="modularity"

FLAG_FILE="-f"
FLAG_LABEL="-l"

EXT_IN=".pairs"

PATH="../../input/"
GRAPH_NAME=(
	"com-amazon.ungraph"
	"com-dblp.ungraph"
	"email-Enron"
	"loc-gowalla_edges"
)

for name in "${GRAPH_NAME[@]}"
do
	RUN_COMMAND="$PROGRAM $FLAG_FILE ${PATH}${name}${EXT_IN} $FLAG_LABEL $METRIC"
	echo $name
	$RUN_COMMAND
done