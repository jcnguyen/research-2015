#!/bin/bash

GRAPH_NAME=( # change
"java ModularityOptimizer ~/Documents/input/p2p-Gnutella08.pairs ~/Documents/output/p2p-Gnutella08 1 1.0 1 1 1 0 1 > ~/Documents/output/p2p-Gnutella08-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/p2p-Gnutella06.pairs ~/Documents/output/p2p-Gnutella06 1 1.0 1 1 1 0 1 > ~/Documents/output/p2p-Gnutella06-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/p2p-Gnutella05.pairs ~/Documents/output/p2p-Gnutella05 1 1.0 1 1 1 0 1 > ~/Documents/output/p2p-Gnutella05-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/p2p-Gnutella04.pairs ~/Documents/output/p2p-Gnutella04 1 1.0 1 1 1 0 1 > ~/Documents/output/p2p-Gnutella04-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/oregon1_010414.pairs ~/Documents/output/oregon1_010414 1 1.0 1 1 1 0 1 > ~/Documents/output/oregon1_010414-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/oregon1_010428.pairs ~/Documents/output/oregon1_010428 1 1.0 1 1 1 0 1 > ~/Documents/output/oregon1_010428-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/oregon1_010505.pairs ~/Documents/output/oregon1_010505 1 1.0 1 1 1 0 1 > ~/Documents/output/oregon1_010505-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/oregon1_010512.pairs ~/Documents/output/oregon1_010512 1 1.0 1 1 1 0 1 > ~/Documents/output/oregon1_010512-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/oregon1_010519.pairs ~/Documents/output/oregon1_010519 1 1.0 1 1 1 0 1 > ~/Documents/output/oregon1_010519-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/oregon1_010421.pairs ~/Documents/output/oregon1_010421 1 1.0 1 1 1 0 1 > ~/Documents/output/oregon1_010421-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/as-22july06.pairs ~/Documents/output/as-22july06 1 1.0 1 1 1 0 1 > ~/Documents/output/as-22july06-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/p2p-Gnutella25.pairs ~/Documents/output/p2p-Gnutella25 1 1.0 1 1 1 0 1 > ~/Documents/output/p2p-Gnutella25-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/p2p-Gnutella30.pairs ~/Documents/output/p2p-Gnutella30 1 1.0 1 1 1 0 1 > ~/Documents/output/p2p-Gnutella30-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/p2p-Gnutella31.pairs ~/Documents/output/p2p-Gnutella31 1 1.0 1 1 1 0 1 > ~/Documents/output/p2p-Gnutella31-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/facebook.pairs ~/Documents/output/facebook 1 1.0 1 1 1 0 1 > ~/Documents/output/facebook-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/email-Enron.pairs ~/Documents/output/email-Enron 1 1.0 1 1 1 0 1 > ~/Documents/output/email-Enron-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/com-amazon.ungraph.pairs ~/Documents/output/com-amazon.ungraph 1 1.0 1 1 1 0 1 > ~/Documents/output/com-amazon.ungraph-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/com-dblp.ungraph.pairs ~/Documents/output/com-dblp.ungraph 1 1.0 1 1 1 0 1 > ~/Documents/output/com-dblp.ungraph-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/amazon0302.pairs ~/Documents/output/amazon0302 1 1.0 1 1 1 0 1 > ~/Documents/output/amazon0302-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/loc-gowalla_edges.pairs ~/Documents/output/loc-gowalla_edges 1 1.0 1 1 1 0 1 > ~/Documents/output/loc-gowalla_edges-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/roadNet-TX.pairs ~/Documents/output/roadNet-TX 1 1.0 1 1 1 0 1 > ~/Documents/output/roadNet-TX-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/web-Stanford.pairs ~/Documents/output/web-Stanford 1 1.0 1 1 1 0 1 > ~/Documents/output/web-Stanford-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/roadNet-PA.pairs ~/Documents/output/roadNet-PA 1 1.0 1 1 1 0 1 > ~/Documents/output/roadNet-PA-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/amazon0312.pairs ~/Documents/output/amazon0312 1 1.0 1 1 1 0 1 > ~/Documents/output/amazon0312-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/amazon0505.pairs ~/Documents/output/amazon0505 1 1.0 1 1 1 0 1 > ~/Documents/output/amazon0505-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/amazon0601.pairs ~/Documents/output/amazon0601 1 1.0 1 1 1 0 1 > ~/Documents/output/amazon0601-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/twitter.pairs ~/Documents/output/twitter 1 1.0 1 1 1 0 1 > ~/Documents/output/twitter-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/roadNet-CA.pairs ~/Documents/output/roadNet-CA 1 1.0 1 1 1 0 1 > ~/Documents/output/roadNet-CA-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/gplus.pairs ~/Documents/output/gplus 1 1.0 1 1 1 0 1 > ~/Documents/output/gplus-lm-performance-maxM.info"
"java ModularityOptimizer ~/Documents/input/soc-pokec-relationships.pairs ~/Documents/output/soc-pokec-relationships 1 1.0 1 1 1 0 1 > ~/Documents/output/soc-pokec-relationships-lm-performance-maxM.info"
)

for name in "${GRAPH_NAME[@]}"
do
	echo ${name}
	${name}
done