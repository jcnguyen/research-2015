#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./convert"

# extensions
EXT_IN=".wpairs"
EXT_OUT=".bin"
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
	FILE_OUTPUT="${LOC}${name}${EXT_OUT}"
	FILE_WEIGHT="${LOC}${name}${EXT_WEIGHT}"
	$PROGRAM -i $FILE_INPUT -o $FILE_OUTPUT -w $FILE_WEIGHT
	# ./convert -i ../../input/{name}.wpairs -o ../../input/{name}.bin -w ../../input/{name}.weights
done