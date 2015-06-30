#!/bin/bash
COMMAND="./convert"

FLAG_INPUT="-i"
FLAG_OUTPUT="-o"
FLAG_RENUM="-r"
FLAG_DISPLAY="-d"
FLAG_WEIGHT="-w"

EXT_IN=".wpairs"
EXT_OUT=".bin"
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
	CONVERT="$COMMAND $FLAG_INPUT ${PATH}${name}${EXT_IN} $FLAG_OUTPUT ${PATH}${name}${EXT_OUT} $FLAG_WEIGHT ${PATH}${name}${EXT_WEIGHT}"
	echo $name
	$CONVERT
done