#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
ALG="-cnm"
METRIC="-modularity"		   # change

# extensions
EXT_IN=".groups"
EXT_OUT=".graph"

# input file info
PATH="../../input/"
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
	"web-Stanford"
	"lesmis"
	"football"
	"celegansneural"
	"soc-sign-Slashdot081106"
	"soc-sign-Slashdot090216"
	"soc-sign-Slashdot090221"
)

for name in "${GRAPH_NAME[@]}"
do
	sort -n < "${PATH}${name}${ALG}${METRIC}${EXT_IN}" > "${PATH}${name}${ALG}${METRIC}${EXT_OUT}"
	# sort -n < ../../input/{name}-cnm{METRIC}.groups > ../../input/{name}-cnm{METRIC}.graph
done