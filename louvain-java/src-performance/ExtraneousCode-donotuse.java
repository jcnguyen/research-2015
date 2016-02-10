/***********************************************************************
 * CODE FROM VOSCLUSTERINGTECHNIQUE.JAVA
 ***********************************************************************/

/**
 * Only called by Modularity Optimizer, to calculate the modularity of the final
 * community groupings.
 *
 * @return qualityFunction - calculating the metric (modularity in this case) for the graph
 */

public double calcModularityFunction()

{
double qualityFunction;
double[] clusterWeight;
int i, j, k;

qualityFunction = 0;
/*for each node in network, if neighbors are in the same cluster, adds their edge weight to total*/
for (i = 0; i < network.nNodes; i++)
{
j = clustering.cluster[i];  /*j is the cluster of node i*/
for (k = network.firstNeighborIndex[i]; k < network.firstNeighborIndex[i + 1]; k++) /*for every neighbor (e    dge?!) that i has??*/
/*if the cluster of the neighbor is the same as cluster of i, add the edge weight of neighbor to quality function*/
if (clustering.cluster[network.neighbor[k]] == j)
qualityFunction += network.edgeWeight[k];  // TODO: does edgeweight take in a vertex or an edge?!
}
qualityFunction += network.totalEdgeWeightSelfLinks;

/*each element of clusterWeight stores the total weight of the nodes in that cluster*/
clusterWeight = new double[clustering.nClusters];
for (i = 0; i < network.nNodes; i++) {
clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
}
/*subtract square of total weights of nodes in each cluser from the quality function*/
for (i = 0; i < clustering.nClusters; i++) {
qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution;
}

qualityFunction /= 2 * network.getTotalEdgeWeight() + network.totalEdgeWeightSelfLinks;

return qualityFunction;
}

/**
 *
 * @param maxNInterations - max iterations you want to run Louvain
 * @return run the algorithm
 */
public boolean runIteratedLouvainAlgorithm(int maxNIterations, int modularityFunction)
{
return runIteratedLouvainAlgorithm(maxNIterations, new Random(), modularityFunction);
}

/**
 *
 * @param maxNInterations - max iterations you want to run Louvain
 * @param random - random number generator
 * @return run the algorithm
 */
public boolean runIteratedLouvainAlgorithm(int maxNIterations, Random random, int modularityFunction)
{
boolean update;
int i;

i = 0;
do
{
update = runLouvainAlgorithm(random, modularityFunction);
i++;
}
while ((i < maxNIterations) && update);
return ((i > 1) || update);
}

public boolean runLouvainAlgorithmWithMultilevelRefinement()
{
return runLouvainAlgorithmWithMultilevelRefinement(new Random());
}

public boolean runLouvainAlgorithmWithMultilevelRefinement(Random random)
{
boolean update, update2;
VOSClusteringTechnique VOSClusteringTechnique;

if (network.nNodes == 1)
return false;

update = runLocalMovingAlgorithm(random);

if (clustering.nClusters < network.nNodes)
{
VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

update2 = VOSClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);

if (update2)
{
update = true;

clustering.mergeClusters(VOSClusteringTechnique.clustering);

runLocalMovingAlgorithm(random);
}
}

return update;
}

public boolean runIteratedLouvainAlgorithmWithMultilevelRefinement(int maxNIterations)
{
return runIteratedLouvainAlgorithmWithMultilevelRefinement(maxNIterations, new Random());
}

public boolean runIteratedLouvainAlgorithmWithMultilevelRefinement(int maxNIterations, Random random)
{
boolean update;
int i;

i = 0;
do
{
update = runLouvainAlgorithmWithMultilevelRefinement(random);
i++;
}
while ((i < maxNIterations) && update);
return ((i > 1) || update);
}


public boolean runSmartLocalMovingAlgorithm()
{
return runSmartLocalMovingAlgorithm(new Random());
}

public boolean runSmartLocalMovingAlgorithm(Random random)
{
boolean update;
int i, j, k;
int[] nNodesPerClusterReducedNetwork;
int[][] nodePerCluster;
Network[] subnetwork;
VOSClusteringTechnique VOSClusteringTechnique;

if (network.nNodes == 1)
return false;

update = runLocalMovingAlgorithm(random);

if (clustering.nClusters < network.nNodes)
{
subnetwork = network.createSubnetworks(clustering);

nodePerCluster = clustering.getNodesPerCluster();

clustering.nClusters = 0;
nNodesPerClusterReducedNetwork = new int[subnetwork.length];
for (i = 0; i < subnetwork.length; i++)
{
VOSClusteringTechnique = new VOSClusteringTechnique(subnetwork[i], resolution);

VOSClusteringTechnique.runLocalMovingAlgorithm(random);

for (j = 0; j < subnetwork[i].nNodes; j++)
clustering.cluster[nodePerCluster[i][j]] = clustering.nClusters + VOSClusteringTechnique.clustering.cluster[j];
clustering.nClusters += VOSClusteringTechnique.clustering.nClusters;
nNodesPerClusterReducedNetwork[i] = VOSClusteringTechnique.clustering.nClusters;
}

VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

i = 0;
for (j = 0; j < nNodesPerClusterReducedNetwork.length; j++)
for (k = 0; k < nNodesPerClusterReducedNetwork[j]; k++)
{
VOSClusteringTechnique.clustering.cluster[i] = j;
i++;
}
VOSClusteringTechnique.clustering.nClusters = nNodesPerClusterReducedNetwork.length;

update |= VOSClusteringTechnique.runSmartLocalMovingAlgorithm(random);

clustering.mergeClusters(VOSClusteringTechnique.clustering);
}

return update;
}


public int removeCluster(int cluster)
{
double maxQualityFunction, qualityFunction;
double[] clusterWeight, totalEdgeWeightPerCluster;
int i, j;

clusterWeight = new double[clustering.nClusters];
totalEdgeWeightPerCluster = new double[clustering.nClusters];
for (i = 0; i < network.nNodes; i++)
{
clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
if (clustering.cluster[i] == cluster)
for (j = network.firstNeighborIndex[i]; j < network.firstNeighborIndex[i + 1]; j++)
totalEdgeWeightPerCluster[clustering.cluster[network.neighbor[j]]] += network.edgeWeight[j];
}

i = -1;
maxQualityFunction = 0;
for (j = 0; j < clustering.nClusters; j++)
if ((j != cluster) && (clusterWeight[j] > 0))
{
qualityFunction = totalEdgeWeightPerCluster[j] / clusterWeight[j];
if (qualityFunction > maxQualityFunction)
{
i = j;
maxQualityFunction = qualityFunction;
}
}

if (i >= 0)
{
for (j = 0; j < network.nNodes; j++)
if (clustering.cluster[j] == cluster)
clustering.cluster[j] = i;
if (cluster == clustering.nClusters - 1)
clustering.nClusters = Arrays2.calcMaximum(clustering.cluster) + 1;
}

return i;
}

public void removeSmallClusters(int minNNodesPerCluster)
{
int i, j, k;
int[] nNodesPerCluster;
VOSClusteringTechnique VOSClusteringTechnique;

VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

nNodesPerCluster = clustering.getNNodesPerCluster();

do
{
i = -1;
j = minNNodesPerCluster;
for (k = 0; k < VOSClusteringTechnique.clustering.nClusters; k++)
if ((nNodesPerCluster[k] > 0) && (nNodesPerCluster[k] < j))
{
i = k;
j = nNodesPerCluster[k];
}

if (i >= 0)
{
j = VOSClusteringTechnique.removeCluster(i);
if (j >= 0)
nNodesPerCluster[j] += nNodesPerCluster[i];
nNodesPerCluster[i] = 0;
}
}
while (i >= 0);

clustering.mergeClusters(VOSClusteringTechnique.clustering);
}


/*
 * Compute the performance of the initial partition of G,
 * in which every node is its own cluster. This is a special
 * case of calculating performance. That is,
 * f(G) = the sum of weights of all self loops,
 * because every node is its own clustering and thus there are
 * no internal edges within communities other than self loops.
 *
 * g(G) = nPossibleEdges - |E|, because at this point
 * every nonexistent edge (disregarding loops) between
 * 2 vertices is between two separate communities, so they
 * count in the g(G) score.
 *
 * g_w is the difference of weight that would be counted if no inter-cluster edges were present,
 * minus the weight that is assigned to the actual inter-cluster edges:
 * g_w = meaningful maximum M * number of nonexistent edges - sum of weights of all nonexistent edges
 *
 * v_scalingParam is the scaling parameter [0,1] that rates the importance of the
 * weight of intercluster edges (with respect to the weight of the
 * intra-cluster edges)
 *
 * Finally, we divide g(G) by the total number of possible edges
 * between all vertices multiplied by the bound on nonexistent edge weight (W).
 *
 * @return perf - performance score
 */
private double initPerformanceFunction() {

// Look at diagonals in the adjacency matrix, calculating f (sum of self loops)
for (int v=0; v<network.nNodes; v++) {
f += adjMatrix[v][v];
}

// inter-community sparsity score
g = M*(network.nPossibleEdges - network.nEdges);

// final performance calculations
double numerator = f + g + v_scalingParam*g_w;
return numerator / denominator;
}

/*
 * Computes the potential performance change from removing
 * vertex from its current community and putting it into its own community.
 *
 * Removing node from comm only changes its relations with
 * vertices in comm. To all other vertices, it is still just
 * some vertex in a different community.
 *
 * @param vertex - the node to be removed
 * @param comm - community containing node
 * @return performance change from removing node from comm
 */
private double removePerfChange(int vertex) {
assert(vertex>=0 && vertex<network.nNodes);
double numeratorChange = 0;
int comm = clustering.cluster[vertex];

// look at all possible edges, existent or nonexistent, with one end at node
for (int u=0; u<network.nNodes; u++) {

// ignore self loops
// if node in our target community
if (u != vertex && clustering.cluster[u] == comm) {

// edge (u, vertex) exists and was previously included in the
// intracluster density score, so its contribution must be removed
if (adjMatrix[vertex][u] != INF) {
numeratorChange -= adjMatrix[vertex][u];
}

// the edge doesn't exist, so increase the intercluster sparsity score
else {
numeratorChange += M;
}
}
}

return numeratorChange/denominator;
}



/*
 * Compute the performance of the initial partition of G,
 * in which every node is its own cluster. This is a special
 * case of calculating performance. That is,
 * f(G) = 0,
 * because every node is its own clustering and thus there are
 * no internal edges within communities.
 *
 * g(G) = nPossibleEdges - |E|, because at this point
 * every nonexistent edge (disregarding loops) between
 * 2 vertices is between two separate communities, so they
 * count in the g(G) score.
 *
 * g_w is the difference of weight that would be counted if no inter-cluster edges were present,
 * minus the weight that is assigned to the actual inter-cluster edges:
 * g_w = meaningful maximum M * number of possible inter-cluster edges
 *       - sum of weights of existent inter-cluster edges
 *
 * v_scalingParam is the scaling parameter [0,1] that rates the importance of the
 * weight of intercluster edges (with respect to the weight of the
 * intra-cluster edges)
 *
 * Finally, we divide g(G) by the total number of possible edges
 * between all vertices multiplied by the bound on nonexistent edge weight (M).
 *
 * @return perf - performance score
 */
private double initPerformanceFunction() {
double  f=0,        // intra-cluster density
g=0;        // inter-cluster sparsity

// inter-community sparsity score
g = M*(network.nPossibleEdges() - network.nEdges);

// final performance calculations
double numerator = f + g + v_scalingParam*g_w;
return numerator / denominator;
}