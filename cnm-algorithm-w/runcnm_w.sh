#!/bin/bash
PROGRAM_NAME="./FastCommunityMH_w"
F_FLAG="-f"
INPUT_PATH="../../input/"
INPUT_NAME=(
	"lesmis.wpairs"
	"football.wpairs"
	"celegansneural.wpairs"
	"polblogs.wpairs"
)
L_FLAG="-l"
LABEL="modularity"

for name in "${INPUT_NAME[@]}"
do
	RUN_COMMAND="$PROGRAM_NAME $F_FLAG ${INPUT_PATH}${name} $L_FLAG $LABEL"
	$RUN_COMMAND
done