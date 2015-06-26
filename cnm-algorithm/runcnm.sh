#!/bin/bash
PROGRAM_NAME="./FastCommunityMH"
F_FLAG="-f"
INPUT_PATH="../../input/"
INPUT_NAME=(
	"adjnoun.pairs"
	"dolphins.pairs"
	"karate.pairs"
	"polbooks.pairs"
)
L_FLAG="-l"
LABEL="modularity"

for name in "${INPUT_NAME[@]}"
do
	RUN_COMMAND="$PROGRAM_NAME $F_FLAG ${INPUT_PATH}${name} $L_FLAG $LABEL"
	$RUN_COMMAND
done