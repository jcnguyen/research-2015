#!/bin/bash
CONVERT="./convert"
COMMUNITY="./community"
FLAG_INPUT="-i"
FLAG_OUTPUT="-o"
FLAG_VERBOSITY="-v"
FLAG_WEIGHTED="-w"
FLAG_LEVEL=-"-l"
LEVEL="-1"
METRIC="-modularity"
EXT_CONV=".wpairs"
EXT_COMM=".bin"
EXT_HIER=".tree"
EXT_WEIGHT=".weights"
INPUT_PATH="../../input/"

NAME=(
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
	"soc-Slashdot0922"
	"soc-sign-epinions"
)

for name in "${NAME[@]}"
do
	CONVERT_COMMAND="$CONVERT $FLAG_INPUT ${INPUT_PATH}${name}${EXT_CONV} $FLAG_OUTPUT ${INPUT_PATH}${name}${EXT_COMM} $FLAG_WEIGHTED ${INPUT_PATH}${name}${EXT_WEIGHT} -d"
	COMMUNITY_COMMAND="${COMMUNITY}${METRIC} ${INPUT_PATH}${name}${EXT_COMM} $FLAG_LEVEL $LEVEL $FLAG_VERBOSITY $FLAG_WEIGHTED ${INPUT_PATH}${name}${EXT_WEIGHT} > ${INPUT_PATH}${name}${METRIC}${EXT_HIER}"
	echo $CONVERT_COMMAND
	$CONVERT_COMMAND
done

