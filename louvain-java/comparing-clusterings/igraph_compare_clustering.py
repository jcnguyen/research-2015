import igraph
import sys

############ READ CLI ARGS ############
if len(sys.argv) == 5:
	input_network_path = sys.argv[1]
	first_clustering_path = sys.argv[2]
	second_clustering_path = sys.argv[3]
	v_raw = sys.argv[4]
else:
	input_network_path = raw_input("Input network edgefile (.pairs, .wpairs): ")
	first_clustering_path = raw_input("Clustering 1 (.graph, .perfgraph): ")
	second_clustering_path = raw_input("Clustering 2 (.graph, .perfgraph): ")
	v_raw = raw_input("Verbose? (0 or 1): ")

if int(v_raw) == 1:
	verbosity = True
else:
	verbosity = False

############ INPUT NETWORK ############
print "Reading input network from file...",
# read in the network: no names, weights if present, and undirected
input_network = igraph.Graph.Read_Ncol(input_network_path, names=False, weights="if_present", directed=False)
if verbosity:
	print input_network
print "finished reading input network."

############ GET VERTEX CLUSTERINGS ############
print "Reading and forming clusterings...",

# read in the first clustering and store in membership list
# the file is in format 'vertex cluster', so we need to extract
# the cluster numbers and put them into the membership list
f = open(first_clustering_path, 'r')
first_clustering_mem_list = []
for line in f:
	split_string = line.split(" ")
	first_clustering_mem_list.append(int(split_string[1].replace("\n", "")))
f.close()

# do the same for second clustering file
f = open(second_clustering_path, 'r')
second_clustering_mem_list = []
for line in f:
	split_string = line.split(" ")
	second_clustering_mem_list.append(int(split_string[1].replace("\n", "")))
f.close()

# create corresponding Vertex Clusterings
first_clustering = igraph.VertexClustering(input_network, first_clustering_mem_list)
second_clustering = igraph.VertexClustering(input_network, second_clustering_mem_list)

print "done creating clusterings."
if verbosity:
	print first_clustering
	print second_clustering

############ COMPARE CLUSTERINGS ############
vi = igraph.compare_communities(first_clustering, second_clustering, method='vi', remove_none=False)
nmi = igraph.compare_communities(first_clustering, second_clustering, method='nmi', remove_none=False)
split_join = igraph.compare_communities(first_clustering, second_clustering, method='split-join', remove_none=False)
rand = igraph.compare_communities(first_clustering, second_clustering, method='rand', remove_none=False)
adj_rand = igraph.compare_communities(first_clustering, second_clustering, method='adjusted_rand', remove_none=False)

print "Meila (2003):", vi
print "Normalized mutual information, Danon et al (2005):", nmi
print "Split-join distance, van Dongen (2000):", split_join
print "Rand index, Rand (1971):", rand
print "Adjusted Rand index, Hubert and Arabie (1986)", adj_rand
