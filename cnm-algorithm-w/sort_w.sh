#!/bin/bash
ALG="-cnm"
METRIC="-modularity"

EXT_IN=".groups"
EXT_OUT=".graph"

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
	sort -n < "${PATH}${name}${ALG}${METRIC}${EXT_IN}" > "${PATH}${name}${ALG}${METRIC}${EXT_OUT}"
done