#!/bin/bash
COMMAND="./community"
ALG="-lm"
METRIC="-modularity"
LEVEL="-1"

FLAG_LEVEL="-l"
FLAG_VERBOSITY="-v"
FLAG_WEIGHT="-w"

EXT_IN=".bin"
EXT_OUT1=".info"
EXT_OUT2=".tree"
EXT_WEIGHT=".weights"

PATH="../../input/"
GRAPH_NAME=(
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
	COMMUNITY="${COMMAND}${METRIC} ${PATH}${name}${EXT_IN} $FLAG_WEIGHT ${PATH}${name}${EXT_WEIGHT} $FLAG_LEVEL $LEVEL $FLAG_VERBOSITY ${PATH}${name}${ALG}${METRIC}${EXT_OUT1} > ${PATH}${name}${ALG}${METRIC}${EXT_OUT2}"
	echo $name
	$COMMUNITY
	
done