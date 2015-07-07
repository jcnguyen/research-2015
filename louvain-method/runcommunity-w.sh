#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./community"
ALG="-lm"
METRIC="-coverage"			   # change
LEVEL="-1"

# flags - see description of alg for additional flags
FLAG_LEVEL="-l"
FLAG_VERBOSITY="-v"
FLAG_WEIGHT="-w"

# extensions
EXT_IN=".bin"
EXT_OUT1=".info"
EXT_OUT2=".tree"
EXT_WEIGHT=".weights"

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
	COMMUNITY="${PROGRAM}${METRIC} ${PATH}${name}${EXT_IN} $FLAG_LEVEL $LEVEL $FLAG_VERBOSITY ${PATH}${name}${ALG}${METRIC}${EXT_OUT1} $FLAG_WEIGHT ${PATH}${name}${EXT_WEIGHT}"
	echo $name
	$COMMUNITY > "${PATH}${name}${ALG}${METRIC}${EXT_OUT2}"
	# ./community{METRIC} ../../input/{name}.bin -l -1 -v ../../input/{name}-lm{METRIC}.info -w ../../input/{name}.weights > ../../input/{name}-lm{METRIC}.tree
	
done