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
LOC="../../input/"			   # change
GRAPH_NAME=(				   # change
	"lesmis"
	"football"
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
	echo $name
	FILE_INPUT="${LOC}${name}${EXT_IN}"
	FILE_OUTPUT1="${LOC}${name}${ALG}${METRIC}${EXT_OUT1}"
	FILE_OUTPUT2="${LOC}${name}${ALG}${METRIC}${EXT_OUT2}"
	FILE_WEIGHT="${LOC}${name}${EXT_WEIGHT}"
	${PROGRAM}${METRIC} $FILE_INPUT -l $LEVEL -v $FILE_OUTPUT1 -w $FILE_WEIGHT > $FILE_OUTPUT2
	# ./community{METRIC} ../../input/{name}.bin -l -1 -v ../../input/{name}-lm{METRIC}.info -w ../../input/{name}.weights > ../../input/{name}-lm{METRIC}.tree
	
done