#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="java ModularityOptimizer"
ALG="-lm"
METRIC="-silhouette" # change

# output files
LOC_OUT="C:/Users/jcngu_000/Documents/lm-silhouette/" # change
EXT_OUT1=".graph"
EXT_OUT2=".info"

# input file
LOC_IN="C:/Users/jcngu_000/Documents/Graphs/" # change
GRAPH_NAME=( # change
"karate.pairs"
"dolphins.pairs"
"lesmis.wpairs"
"polbooks.pairs"
"adjnoun.pairs"
"football.wpairs"
"netscience.wpairs"
"power.pairs"
"hep-th.wpairs"
"oregon1_010414.pairs"
"oregon1_010421.pairs"
"oregon1_010428.pairs"
"oregon1_010505.pairs"
"oregon1_010512.pairs"
"oregon1_010519.pairs"
"astro-ph.wpairs"
"cond-mat.wpairs"
"as-22july06.pairs"
"cond-mat-2003.wpairs"
"email-Enron.pairs"
"cond-mat-2005.wpairs"
"loc-gowalla_edges.pairs"
"com-dblp.ungraph.pairs"
"com-amazon.ungraph.pairs"
"roadNet-PA.pairs"
"roadNet-TX.pairs"
"roadNet-CA.pairs"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_IN="${LOC_IN}${name}"
	FILE_OUT1="${LOC_OUT}${name}${ALG}${METRIC}${EXT_OUT1}"
	FILE_OUT2="${LOC_OUT}${name}${ALG}${METRIC}${EXT_OUT2}"
	${PROGRAM} $FILE_IN $FILE_OUT1 3 1.0 1 1 10 0 1 > $FILE_OUT2
	# java ModularityOptimizer C:/Users/jcngu_000/Documents/Graphs/netscience.wpairs C:/Users/jcngu_000/Documents/lm-silhouette/netscience-lm-silhouette.graph 3 1.0 1 1 10 0 1 > C:/Users/jcngu_000/Documents/lm-silhouette/netscience-lm-silhouette.info
done