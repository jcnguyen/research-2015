#!/bin/bash
PROGRAM="./FastCommunityMH"
METRIC="modularity"

FLAG_FILE="-f"
FLAG_LABEL="-l"

EXT_IN=".pairs"

PATH="../../input/"
GRAPH_NAME=(
	"adjnoun"
	"amazon0302"
	"amazon0312"
	"amazon0505"
	"amazon0601"
	"as-22july06"
	"com-amazon.ungraph"
	"com-dblp.ungraph"
	"dolphins"
	"email-Enron"
	"karate"
	"loc-gowalla_edges"
	"oregon1_010414"
	"oregon1_010421"
	"oregon1_010428"
	"oregon1_010505"
	"oregon1_010512"
	"oregon1_010519"
	"p2p-Gnutella04"
	"p2p-Gnutella05"
	"p2p-Gnutella06"
	"p2p-Gnutella08"
	"p2p-Gnutella25"
	"p2p-Gnutella30"
	"p2p-Gnutella31"
	"polbooks"
	"power"
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
	RUN_COMMAND="$PROGRAM $FLAG_FILE ${PATH}${name}${EXT_IN} $FLAG_LABEL $LABEL"
	echo $name
	$RUN_COMMAND
done