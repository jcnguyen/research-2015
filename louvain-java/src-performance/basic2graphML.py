import sys
import string
from sets import Set
#
# Converts from input graph and communities to graphML format
#

#
# The following needs to be modified for each run
#

numNodes = 34
numEdges = 78

# UNCOMMENT ONE OF THESE
# alg = "lm"
# metrics = Set(['modularity', 'silhouette', 'coverage', 'performance-maxM', 'performance-maxM-opt'])

alg = "cnm"
metrics = Set(['modularity', 'coverage'])

for metric in metrics:
	inputGraph = "/Users/christina_tong/Box Sync/ComplexNetworks-CREU/Graphs/karate.pairs"

	metric_no_suffix = metric.split('-')[0]
	inputComm = "/Users/christina_tong/Box Sync/ComplexNetworks-CREU/Data/" + alg + "-" + \
		metric_no_suffix + "/karate-" + alg + "-" + metric + ".graph"
	outfilePath = "./temp-gml/karate-" + alg + "-" + metric + ".graphml"

	#
	# The rest of the file should not need to be modified
	#
	infileGraph = open(inputGraph, "rU")
	infileCommunity = open(inputComm, "rU")
	outfile = open(outfilePath, "w")

	try:
	    communities = infileCommunity.readlines()
	    graph = infileGraph.readlines()
	finally:
	    infileGraph.close()
	    infileCommunity.close()

	xmlHeader = "<?xml version= \"1.0\" encoding= \"UTF-8\"?> "

	graphMLHeader = """
	<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"  
	    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"
	    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">
	"""

	graphHeader = "  <graph id=\"G\" edgedefault=\"undirected\">\n"

	# key for community identification
	keyHeader = """
	  <key id=\"community\" for=\"node\" attr.name=\"community\" attr.type=\"int\">
	    <default>0</default>
	  </key>
	"""

	# Reading in the community information into dict
	communityDict = {}
	for l in communities:
	    info = l.split()
	    communityDict[info[0]] = info[1] 

	# Processing files, writing output

	outfile.write(xmlHeader)
	outfile.write(graphMLHeader)
	outfile.write(keyHeader)

	outfile.write(graphHeader)

	edgeCount = 0
	for i in range(numNodes):
	    outfile.write("  <node id=\"n" + str(i) + "\">\n")
	    outfile.write("    <data key=\"community\">" + communityDict[str(i)] + "</data>\n")
	    outfile.write("  </node>\n")
	    
	for e in graph:
	    edge = e.split()
	    outfile.write("  <edge id=\"e"+str(edgeCount)+"\" source=\"n"+edge[0]+"\" target=\"n"+edge[1]+"\"/>\n")
	    edgeCount = edgeCount + 1

	outfile.write("</graph>\n")
	outfile.write("</graphml>\n")

	if edgeCount != numEdges:
	    print "number of lines in file != number of edges expected"

	# clean up
	outfile.close()
