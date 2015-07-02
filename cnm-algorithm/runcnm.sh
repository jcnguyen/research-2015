#!/bin/bash
PROGRAM="./FastCommunityMH"
METRIC="modularity"

FLAG_FILE="-f"
FLAG_LABEL="-l"

EXT_IN=".pairs"

PATH="../../input/"
GRAPH_NAME=(
	"amazon0302"
	"amazon0312"
	"amazon0505"
	"amazon0601"
	"as-22july06"
	"roadNet-CA"
	"roadNet-PA"
	"roadNet-TX"
	"soc-LiveJournal1"
	"soc-pokec-relationships"
	"web-NotreDame"
	"web-Stanford"
)

for name in "${GRAPH_NAME[@]}"
do
	RUN_COMMAND="$PROGRAM $FLAG_FILE ${PATH}${name}${EXT_IN} $FLAG_LABEL $METRIC"
	echo $name
	$RUN_COMMAND
done