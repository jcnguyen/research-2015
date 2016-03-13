#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="java ModularityOptimizer"
ALG="-lm"
METRIC="-silhouette" # change

# file info
LOC_IN="~/Documents/input/" # change
LOC_OUT="~/Documents/output/" # change
EXT_IN=".woairs" # change
EXT_OUT1=".graph"
EXT_OUT2=".info"
GRAPH_NAME=( # change
"celegansneural"
"netscience"
"hep-th"
"astro-ph"
"cond-mat"
"cond-mat-2003"
"cond-mat-2005"
"soc-sign-Slashdot081106"
"soc-sign-Slashdot090216"
"soc-sign-Slashdot090221"
)

for name in "${GRAPH_NAME[@]}"
do
	echo ${name}
	FILE_IN="${LOC_IN}${name}${EXT_IN}"
	FILE_OUT1="${LOC_OUT}${name}${ALG}${METRIC}${EXT_OUT1}"
	FILE_OUT2="${LOC_OUT}${name}${ALG}${METRIC}${EXT_OUT2}"
	${PROGRAM} $FILE_IN $FILE_OUT1 3 1.0 1 1 10 0 1 > $FILE_OUT2
done