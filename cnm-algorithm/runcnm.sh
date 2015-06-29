#!/bin/bash
PROGRAM_NAME="./FastCommunityMH"
F_FLAG="-f"
INPUT_PATH="../../input/"
INPUT_NAME=(
	"p2p-Gnutella08.pairs"
	"p2p-Gnutella06.pairs"
	"p2p-Gnutella05.pairs"
	"oregon1_010414.pairs"
	"oregon1_010421.pairs"
	"p2p-Gnutella04.pairs"
	"oregon1_010428.pairs"
	"oregon1_010505.pairs"
	"oregon1_010512.pairs"
	"oregon1_010519.pairs"
	"p2p-Gnutella25.pairs"
)
L_FLAG="-l"
LABEL="modularity"

for name in "${INPUT_NAME[@]}"
do
	RUN_COMMAND="$PROGRAM_NAME $F_FLAG ${INPUT_PATH}${name} $L_FLAG $LABEL"
	$RUN_COMMAND
done