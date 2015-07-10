#!/bin/bash
# verify that the areas with "# change" marker are correct

# program info
PROGRAM="./FastCommunityMH_w"
METRIC="modularity"		  	   # change

# extensions
EXT_IN=".wpairs"

# input file info
LOC="../../input/"			   # change
GRAPH_NAME=(				   # change
	"soc-Slashdot0902"
)

for name in "${GRAPH_NAME[@]}"
do
	echo $name
	FILE_INPUT="${LOC}${name}${EXT_IN}"
	$PROGRAM -f $FILE_INPUT -l $METRIC
	# ./FastCommunityMH_w -f ../../input/{name}.wpairs -l {METRIC}
done