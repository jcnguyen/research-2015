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
	"soc-pokec-relationships"
	"web-Stanford"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_INPUT="${LOC}${name}${ALG}${METRIC}${EXT_IN}"
	FILE_OUTPUT="${LOC}${name}${ALG}${METRIC}${EXT_OUT}"
	$PROGRAM $FILE_INPUT > $FILE_OUTPUT
	# ./display ../../input/{name}-lm{METRIC}.tree > ../../input/{name}-lm{METRIC}.graph
done