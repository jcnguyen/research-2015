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
	"lesmis"
	"football"
	"celegansneural"
	"polblogs"
	"netscience"
	"hep-th"
	"astro-ph"
	"cond-mat"
	"cond-mat-2003"
	"cond-mat-2005"
	"soc-sign-Slashdot081106"
	"soc-sign-Slashdot090216"
	"soc-sign-Slashdot090221"
	"soc-Slashdot0902"
	"soc-sign-epinions"
)

for name in "${GRAPH_NAME[@]}"
do
	sort -n < "${PATH}${name}${ALG}${METRIC}${EXT_IN}" > "${PATH}${name}${ALG}${METRIC}${EXT_OUT}"
	# sort -n < ../../input/{name}-cnm{METRIC}.groups > ../../input/{name}-cnm{METRIC}.graph
done