#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./convert"

# flags - see description of alg for additional flags
FLAG_INPUT="-i"
FLAG_OUTPUT="-o"
FLAG_RENUM="-r"
FLAG_DISPLAY="-d"
FLAG_WEIGHT="-w"

# extensions
EXT_IN=".wpairs"
EXT_OUT=".bin"
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
	CONVERT="$PROGRAM $FLAG_INPUT ${PATH}${name}${EXT_IN} $FLAG_OUTPUT ${PATH}${name}${EXT_OUT} $FLAG_WEIGHT ${PATH}${name}${EXT_WEIGHT}"
	echo $name
	$CONVERT
	# ./convert -i ../../input/{name}.wpairs -o ../../input/{name}.bin -w ../../input/{name}.weights
done