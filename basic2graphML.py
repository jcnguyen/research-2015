import sys
import string

#
# Converts from input graph and communities to graphML format
#

#
# The following may need to be modified for each run
#

dir = "../Data/2015-06-25/"
graphName = "as-22july06"
numNodes = 22963
numEdges = 48437
algorithmName = "-LM-"
metricName = "modularity"

inputGraphSuffix = ".txt"
inputCommunitySuffix = ".graph"

#
# The rest of the file should not need to be modified
#

outputSuffix = ".graphML"


infileGraph = open(dir+graphName+inputGraphSuffix, "rU")
infileCommunity = open(dir+graphName+algorithmName+metricName+inputCommunitySuffix, "rU")
outfile = open(dir+graphName+algorithmName+metricName+outputSuffix, "w")

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
