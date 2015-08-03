#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./FastCommunityMH"
METRIC="coverage"		  	   # change

# extensions
EXT_IN=".pairs"

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

	FILE_INPUT="${LOC}${name}${EXT_IN}"
	$PROGRAM -f $FILE_INPUT -l $METRIC "-c "
	# ./FastCommunityMH -f ../../input/{name}.pairs -l {METRIC}
done