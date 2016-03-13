#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="java ModularityOptimizer"
ALG="-lm"
METRIC="-silhouette" # change

# file info
LOC_IN="~/Documents/input/" # change
LOC_OUT="~/Documents/output/" # change
EXT_IN=".pairs" # change
EXT_OUT1=".graph"
EXT_OUT2=".info"
GRAPH_NAME=( # change
"power"
"p2p-Gnutella08"
"p2p-Gnutella06"
"p2p-Gnutella05"
"p2p-Gnutella04"
"oregon1_010414"
"oregon1_010428"
"oregon1_010505"
"oregon1_010512"
"oregon1_010519"
"oregon1_010421"
"as-22july06"
"p2p-Gnutella25"
"p2p-Gnutella30"
"p2p-Gnutella31"
"facebook"
"email-Enron"
"com-amazon.ungraph"
"com-dblp.ungraph"
"amazon0302"
"loc-gowalla_edges"
"roadNet-TX"
"web-Stanford"
"roadNet-PA"
"amazon0312"
"amazon0505"
"amazon0601"
"twitter"
"roadNet-CA"
"gplus"
"soc-pokec-relationships"
)

for name in "${GRAPH_NAME[@]}"
do
	echo ${name}
	FILE_IN="${LOC_IN}${name}${EXT_IN}"
	FILE_OUT1="${LOC_OUT}${name}${ALG}${METRIC}${EXT_OUT1}"
	FILE_OUT2="${LOC_OUT}${name}${ALG}${METRIC}${EXT_OUT2}"
	${PROGRAM} $FILE_IN $FILE_OUT1 3 1.0 1 1 10 0 1 > $FILE_OUT2
done