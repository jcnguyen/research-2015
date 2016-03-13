#!/bin/bash

GRAPH_NAME=( # change
java ModularityOptimizer ~/Documents/input/cond-mat-2003.wpairs ~/Documents/output/cond-mat-2003 1 1.0 35.2 1 1 0 1 > ~/Documents/output/cond-mat-2003-lm-performance-maxM.info
java ModularityOptimizer ~/Documents/input/cond-mat-2005.wpairs ~/Documents/output/cond-mat-2005 1 1.0 46 1 1 0 1 > ~/Documents/output/cond-mat-2005-lm-performance-maxM.info
java ModularityOptimizer ~/Documents/input/cond-mat-2003.wpairs ~/Documents/output/cond-mat-2003 1 1.0 0.5 1 1 0 1 > ~/Documents/output/cond-mat-2003-lm-performance-rehoM.info
java ModularityOptimizer ~/Documents/input/cond-mat-2005.wpairs ~/Documents/output/cond-mat-2005 1 1.0 0.5 1 1 0 1 > ~/Documents/output/cond-mat-2005-lm-performance-rehoM.info
)

for name in "${GRAPH_NAME[@]}"
do
	echo ${name}
	${name}
done