/**
 * VOSClusteringTechnique
 *
 * Contains the algorithms to maximize metric values. Also includes the metric
 * functions.
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/23/14
 *
 * @author Complex Network Research
 *         http://research.pomona.edu/complexnetworks/
 * @version 11/18/15
 *
 * The Floyd-Warshall algorithm for APSP problem was found at geeksforgeeks.org.
 */

import java.util.Random;
import java.util.Arrays;

public class VOSClusteringTechnique {
    private static final boolean TEST = false;
    private static final boolean TEST2 = false;

    // metric functions
    private static final int MODULARITY_STANDARD = 1;
    private static final int MODULARITY_ALTERNATIVE = 2;
    private static final int SILHOUETTE_INDEX = 3;
    private static final double INF = Double.POSITIVE_INFINITY;

    protected Network network;
    protected Clustering clustering;
    protected double resolution;

    /**
     * Constructor where every vertex is its own cluster.
     *
     * @param network     object representation of graph
     * @param resolution  granularity level at which communities are detected;
     *                    1.0 for standard modularity-based community detection
     */
    public VOSClusteringTechnique(Network network, double resolution) {
        this.network = network;

        // on initialization, every node is in its own community
        clustering = new Clustering(network.nNodes);
        clustering.initSingletonClusters();

        this.resolution = resolution;
    }

    /**
     * Standard constructor.
     *
     * @param network     object representation of graph
     * @param clustering  object representation of the communities in graph
     * @param resolution  granularity level at which communities are detected
     */
    public VOSClusteringTechnique(
        Network network, Clustering clustering, double resolution) {

        this.network = network;
        this.clustering = clustering;
        this.resolution = resolution;
    }

    /**
     * @return the object representation of a graph
     */
    public Network getNetwork() {
        return network;
    }

    /**
     * 
     * @return the object representation of the communities in graph
     */
    public Clustering getClustering() {
        return clustering;
    }

    /**
     * 
     * @return the granularity level at which communities are detected
     */
    public double getResolution() {
        return resolution;
    }

    /**
     * Set the network.
     *
     * @param network  object representation of graph
     */
    public void setNetwork(Network network) {
        this.network = network;
    }

    /**
     * Set the clustering.
     *
     * @param clustering  object representation of the communities in graph 
     */
    public void setClustering(Clustering clustering) {
        this.clustering = clustering;
    }
    
    /**
     * Set the resolution.
     *
     * @param resolution  granularity level at which communities are detected 
     */
    public void setResolution(double resolution) {
        this.resolution = resolution;
    }


    /**
     * Run the Louvain algorithm using a defined metric function.
     *
     * @param metricFunction  the metric function
     * @return true if nodes are moved to a different community; false otherwise 
     */
    public boolean runLouvainAlgorithm(int metricFunction) {
        return runLouvainAlgorithm(new Random(), metricFunction);
    }

    /**
     * Run the Louvain algorithm using a defined metric function.
     *
     * @param random          random num generator
     * @param metricFunction  the metric function
     * @return true if nodes are moved to a different community; false otherwise 
     */
    public boolean runLouvainAlgorithm(Random random, int metricFunction) {
        boolean update, update2;
        
        VOSClusteringTechnique VOSClusteringTechnique;

        if (network.nNodes == 1)
            return false;

        // run an iteration of the Louvain algorithm based on the metric
        if (metricFunction == MODULARITY_STANDARD || 
            metricFunction == MODULARITY_ALTERNATIVE) {
            update = runLocalMovingAlgorithmWithModularity(random);
        } else if (metricFunction == SILHOUETTE_INDEX) {
            update = runLocalMovingAlgorithmWithSilhouetteIndex(random);
        } else {
            update = runLocalMovingAlgorithmWithModularity(random);
        }

        if (clustering.nClusters < network.nNodes) {
            VOSClusteringTechnique = new VOSClusteringTechnique(
                network.createReducedNetwork(clustering), resolution);

            update2 = VOSClusteringTechnique.runLouvainAlgorithm(
                random, metricFunction);
            if (update2) {
                update = true;

                clustering.mergeClusters(VOSClusteringTechnique.clustering);
            }
        }

        return update;
    }

    public boolean runIteratedLouvainAlgorithm(
        int maxNIterations, int metricFunction) {
        return runIteratedLouvainAlgorithm(
            maxNIterations, new Random(), metricFunction);
    }

    public boolean runIteratedLouvainAlgorithm(
        int maxNIterations, Random random, int metricFunction) {
        boolean update;
        int i;

        i = 0;
        do {
            update = runLouvainAlgorithm(random, metricFunction);
            i++;
        }
        while ((i < maxNIterations) && update);
        return ((i > 1) || update);
    }

    public boolean runLouvainAlgorithmWithMultilevelRefinement() {
        return runLouvainAlgorithmWithMultilevelRefinement(new Random());
    }

    public boolean runLouvainAlgorithmWithMultilevelRefinement(Random random) {
        boolean update, update2;
        VOSClusteringTechnique VOSClusteringTechnique;

        if (network.nNodes == 1)
            return false;

        update = runLocalMovingAlgorithmWithModularity(random);

        if (clustering.nClusters < network.nNodes) {
            VOSClusteringTechnique = new VOSClusteringTechnique(
                network.createReducedNetwork(clustering), resolution);

            update2 = VOSClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);

            if (update2) {
                update = true;

                clustering.mergeClusters(VOSClusteringTechnique.clustering);

                runLocalMovingAlgorithmWithModularity(random);
            }
        }

        return update;
    }

    public boolean runIteratedLouvainAlgorithmWithMultilevelRefinement(
        int maxNIterations) {
        return runIteratedLouvainAlgorithmWithMultilevelRefinement(
            maxNIterations, new Random());
    }

    public boolean runIteratedLouvainAlgorithmWithMultilevelRefinement(
        int maxNIterations, Random random) {
        boolean update;
        int i;

        i = 0;
        do {
            update = runLouvainAlgorithmWithMultilevelRefinement(random);
            i++;
        }
        while ((i < maxNIterations) && update);
        return ((i > 1) || update);
    }

    public boolean runSmartLocalMovingAlgorithm() {
        return runSmartLocalMovingAlgorithm(new Random());
    }

    public boolean runSmartLocalMovingAlgorithm(Random random) {
        boolean update;
        int i, j, k;
        int[] nNodesPerClusterReducedNetwork;
        int[][] nodePerCluster;
        Network[] subnetwork;
        VOSClusteringTechnique VOSClusteringTechnique;

        if (network.nNodes == 1)
            return false;

        update = runLocalMovingAlgorithmWithModularity(random);

        if (clustering.nClusters < network.nNodes) {
            subnetwork = network.createSubnetworks(clustering);

            nodePerCluster = clustering.getNodesPerCluster();

            clustering.nClusters = 0;
            nNodesPerClusterReducedNetwork = new int[subnetwork.length];
            for (i = 0; i < subnetwork.length; i++) {
                VOSClusteringTechnique = new VOSClusteringTechnique(
                    subnetwork[i], resolution);

                VOSClusteringTechnique.runLocalMovingAlgorithmWithModularity(
                    random);

                for (j = 0; j < subnetwork[i].nNodes; j++)
                    clustering.cluster[nodePerCluster[i][j]] = clustering.nClusters + VOSClusteringTechnique.clustering.cluster[j];
                clustering.nClusters += VOSClusteringTechnique.clustering.nClusters;
                nNodesPerClusterReducedNetwork[i] = VOSClusteringTechnique.clustering.nClusters;
            }

            VOSClusteringTechnique = new VOSClusteringTechnique(
                network.createReducedNetwork(clustering), resolution);

            i = 0;
            for (j = 0; j < nNodesPerClusterReducedNetwork.length; j++)
                for (k = 0; k < nNodesPerClusterReducedNetwork[j]; k++) {
                    VOSClusteringTechnique.clustering.cluster[i] = j;
                    i++;
                }
            VOSClusteringTechnique.clustering.nClusters = nNodesPerClusterReducedNetwork.length;

            update |= VOSClusteringTechnique.runSmartLocalMovingAlgorithm(random);

            clustering.mergeClusters(VOSClusteringTechnique.clustering);
        }

        return update;
    }

    public boolean runIteratedSmartLocalMovingAlgorithm(int nIterations) {
        return runIteratedSmartLocalMovingAlgorithm(nIterations, new Random());
    }

    public boolean runIteratedSmartLocalMovingAlgorithm(
        int nIterations, Random random) {

        boolean update;
        int i;

        update = false;
        for (i = 0; i < nIterations; i++)
            update |= runSmartLocalMovingAlgorithm(random);
        return update;
    }

    public int removeCluster(int cluster) {
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, totalEdgeWeightPerCluster;
        int i, j;

        clusterWeight = new double[clustering.nClusters];
        totalEdgeWeightPerCluster = new double[clustering.nClusters];
        for (i = 0; i < network.nNodes; i++) {
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
            if (clustering.cluster[i] == cluster)
                for (j = network.firstNeighborIndex[i]; j < network.firstNeighborIndex[i + 1]; j++)
                    totalEdgeWeightPerCluster[clustering.cluster[network.neighbor[j]]] += network.edgeWeight[j];
        }

        i = -1;
        maxQualityFunction = 0;
        for (j = 0; j < clustering.nClusters; j++)
            if ((j != cluster) && (clusterWeight[j] > 0)) {
                qualityFunction = totalEdgeWeightPerCluster[j] / clusterWeight[j];
                if (qualityFunction > maxQualityFunction) {
                    i = j;
                    maxQualityFunction = qualityFunction;
                }
            }

        if (i >= 0) {
            for (j = 0; j < network.nNodes; j++)
                if (clustering.cluster[j] == cluster)
                    clustering.cluster[j] = i;
            if (cluster == clustering.nClusters - 1)
                clustering.nClusters = Arrays2.calcMaximum(clustering.cluster) + 1;
        }

        return i;
    }

    public void removeSmallClusters(int minNNodesPerCluster) {
        int i, j, k;
        int[] nNodesPerCluster;
        VOSClusteringTechnique VOSClusteringTechnique;

        VOSClusteringTechnique = new VOSClusteringTechnique(
            network.createReducedNetwork(clustering), resolution);

        nNodesPerCluster = clustering.getNNodesPerCluster();

        do {
            i = -1;
            j = minNNodesPerCluster;
            for (k = 0; k < VOSClusteringTechnique.clustering.nClusters; k++)
                if ((nNodesPerCluster[k] > 0) && (nNodesPerCluster[k] < j)) {
                    i = k;
                    j = nNodesPerCluster[k];
                }

            if (i >= 0) {
                j = VOSClusteringTechnique.removeCluster(i);
                if (j >= 0)
                    nNodesPerCluster[j] += nNodesPerCluster[i];
                nNodesPerCluster[i] = 0;
            }
        }
        while (i >= 0);

        clustering.mergeClusters(VOSClusteringTechnique.clustering);
    }

    /***********************************************************************
     * MODULARITY
     ***********************************************************************/

    /**
     * Calculate the modularity score of a graph.
     *
     * Note that this method is only called by ModularityOptimizer to 
     * calculate the modularity of the final clustering/community graph.
     * 
     * @return the modularity score of a graph
     */
    public double calcModularityFunction() {
        double qualityFunction;
        double[] clusterWeight;
        int i, j, k;

        // calculate the total edge weight of the network
        qualityFunction = 0;
        for (i = 0; i < network.nNodes; i++) {
            j = clustering.cluster[i];
            for (k = network.firstNeighborIndex[i]; k < network.firstNeighborIndex[i + 1]; k++)
                if (clustering.cluster[network.neighbor[k]] == j) 
                    qualityFunction += network.edgeWeight[k];
        }
        qualityFunction += network.totalEdgeWeightSelfLinks;

        // store the total node weight of each cluster
        clusterWeight = new double[clustering.nClusters];
        for (i = 0; i < network.nNodes; i++) {
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
        }

        // calculate the modularity score
        for (i = 0; i < clustering.nClusters; i++) {
            qualityFunction -= clusterWeight[i] * clusterWeight[i] * resolution;
        }
        qualityFunction /= (2 * network.getTotalEdgeWeight() + network.totalEdgeWeightSelfLinks);

        return qualityFunction;
    }

    /**
     * Find the partition of nodes that maximizes the modularity score.
     *
     * Note that while this method only returns a boolean, its work is saved in 
     * the instance variable clustering.
     *
     * @param random  a random num generator 
     * @return true if nodes are moved to a different community 
     *         false otherwise 
     */
    public boolean runLocalMovingAlgorithmWithModularity(Random random) {
        boolean update;
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, edgeWeightPerCluster;
        int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
        int[] neighboringCluster, newCluster, nNodesPerCluster, nodePermutation, unusedCluster;

        update = false;

        if (network.nNodes == 1)
            return update;


        // store the total degree of and the number of nodes in each cluster
        clusterWeight = new double[network.nNodes];
        nNodesPerCluster = new int[network.nNodes]; 
        for (i = 0; i < network.nNodes; i++) {
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
            nNodesPerCluster[clustering.cluster[i]]++;
        }

        // store the empty clusters
        nUnusedClusters = 0;
        unusedCluster = new int[network.nNodes];
        for (i = 0; i < network.nNodes; i++) {
            if (nNodesPerCluster[i] == 0) {
                unusedCluster[nUnusedClusters] = i;
                nUnusedClusters++;
            }
        }

        nodePermutation = Arrays2.generateRandomPermutation(network.nNodes, random);
        edgeWeightPerCluster = new double[network.nNodes];
        neighboringCluster = new int[network.nNodes - 1];
        nStableNodes = 0;
        i = 0;

        // compute the clustering that achieves the highest modularity score
        do {
            j = nodePermutation[i]; // start with a random vertex j

            // get the neighboring clusters of vertex j and their edge weights
            nNeighboringClusters = 0;
            for (k = network.firstNeighborIndex[j]; k < network.firstNeighborIndex[j + 1]; k++) {
                l = clustering.cluster[network.neighbor[k]];
                if (edgeWeightPerCluster[l] == 0) { 
                    neighboringCluster[nNeighboringClusters] = l;
                    nNeighboringClusters++;
                }
                edgeWeightPerCluster[l] += network.edgeWeight[k];
            }

            // remove vertex j from its cluster and update structures 
            clusterWeight[clustering.cluster[j]] -= network.nodeWeight[j]; 
            nNodesPerCluster[clustering.cluster[j]]--;
            if (nNodesPerCluster[clustering.cluster[j]] == 0) {
                unusedCluster[nUnusedClusters] = clustering.cluster[j];
                nUnusedClusters++;
            }

            // prepare to simulate adding vertex j to each neighboring cluster
            bestCluster = -1;
            maxQualityFunction = 0;

            // determine the cluster that vertex j should be added to
            for (k = 0; k < nNeighboringClusters; k++) {

                // move vertex j into cluster l and calculate modularity score
                l = neighboringCluster[k];
                qualityFunction = edgeWeightPerCluster[l] - (
                    network.nodeWeight[j] * clusterWeight[l] * resolution);
                
                // adding vertex j to neighboring cluster l gives a 
                // better modularity score
                // TODO tie-breaking step differs from SI version
                if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) && (l < bestCluster))) {
                    bestCluster = l;
                    maxQualityFunction = qualityFunction;
                }
                edgeWeightPerCluster[l] = 0;

            }
            // the original cluster of vertex j gives a better modularity score
            if (maxQualityFunction == 0) {
                bestCluster = unusedCluster[nUnusedClusters - 1];
                nUnusedClusters--;
            }

            // add vertex j to the cluster that gives the best modularity score
            // and update structures accordingly
            clusterWeight[bestCluster] += network.nodeWeight[j];
            nNodesPerCluster[bestCluster]++;
            if (bestCluster == clustering.cluster[j]) { 
                nStableNodes++;
            } else { 
                clustering.cluster[j] = bestCluster;
                nStableNodes = 1;
                update = true;
            }

            i = (i < network.nNodes - 1) ? (i + 1) : 0;

        } 
        while (nStableNodes < network.nNodes);

        // update structures to reflect the new clustering
        newCluster = new int[network.nNodes];
        clustering.nClusters = 0;
        for (i = 0; i < network.nNodes; i++) {
            if (nNodesPerCluster[i] > 0) {
                newCluster[i] = clustering.nClusters;
                clustering.nClusters++;
            }
        }
        for (i = 0; i < network.nNodes; i++) {
            clustering.cluster[i] = newCluster[clustering.cluster[i]];
        }

        return update;
    }

    /***********************************************************************
     * SILHOUETTE INDEX
     ***********************************************************************/

    /**
     * Calculate the distance of the shortest path between all vertices.
     *
     * @param graph   a matrix representing the input network
     * @param nNodes  the number of vertices/nodes in the network
     * @return a matrix of shortest path distances
     */
    public double[][] floydWarshall(double[][] graph, int nNodes) {
        double[][] distances;
        int i, j, k;

        // initialize the solution matrix to the input graph matrix
        distances = new double[nNodes][nNodes];
        for (i = 0; i < nNodes; i++) {
            for (j = 0; j < nNodes; j++) {
                if (i == j) { // diagonal
                    distances[i][j] = 0;
                } else if (graph[i][j] == INF) {
                    distances[i][j] = 1000000000;
                } else {
                    distances[i][j] = graph[i][j];
                }
            }
        }

        // calculate the shortest distance between vertices i and j, 
        // given an intermediate vertex k
        for (k = 0; k < nNodes; k++) { 
            for (i = 0; i < nNodes; i++) {
                for (j = 0; j < nNodes; j++) { 
                    // if vertex k is on the shortest path from
                    // i to j, then update the value of distances[i][j]
                    if (distances[i][k] + distances[k][j] < distances[i][j]) {
                        distances[i][j] = distances[i][k] + distances[k][j];
                    }
                }
            }
        }

        return distances;
    }

    /**
     * Calculate the silhouette index score of a graph.
     * 
     * @return the silhouette index score of a graph
     */
    public double calcSilhouetteFunction() {
        return calcSilhouetteFunction(
            floydWarshall(network.getMatrix(), network.nNodes));
    }

    /**
     * Calculate the silhouette index score of a graph.
     * 
     * @param shortestPath  the distances of the shortest path between 
     *                      all vertices
     * @return the silhouette index score of a graph
     */
    public double calcSilhouetteFunction(double[][] shortestPath) {
        int c, c2, i, j, k;
        int cumulativeNNodes, nClusters, nNodesInClusterC, nNodesInClusterC2, nEmptyClusters;
        int[] nNodesPerCluster;
        int[][] nodesPerCluster;
        double averageInnerDistance, distance, minAverageOuterDistance, silhouetteIndexGlobal;
        double[] averageOuterDistances, silhouetteIndexPerCluster, silhouetteWidths;

        nClusters = clustering.nClusters;
        nNodesPerCluster = clustering.getNNodesPerCluster();
        nodesPerCluster = clustering.getNodesPerCluster();

        silhouetteIndexPerCluster = new double[nClusters];

        if (TEST) {
            System.out.println("________________________________________________________");
        }

        // iterate through the clusters
        cumulativeNNodes = 0;
        nEmptyClusters = 0;
        for (c = 0; c < nClusters; c++) {
            nNodesInClusterC = nNodesPerCluster[c];

            if (TEST) {
                System.out.println();
                System.out.println("CLUSTER " + c);
                System.out.println("-----------------------------------------");
                System.out.println("nNodesInClusterC: " + nNodesInClusterC);
            }

            // calculate the silhouette index of the entire cluster c
            if (TEST) {
                System.out.println("1. calculating silhouette of vertex j in cluster " + c);
            }
            if (nNodesInClusterC == 0) {
                silhouetteIndexPerCluster[c] = 0;
                nEmptyClusters++;
            } else {
                silhouetteWidths = new double[nNodesInClusterC];
                for(j = 0; j < nNodesInClusterC; j++) {
                    if (TEST) {
                        System.out.println("    node " + (j + cumulativeNNodes));
                    }

                    // calculate the average of the shortest distances between 
                    // vertex j and every other vertex k in the same cluster
                    if (TEST) {
                        System.out.println("        -calculating average inner distance");
                    }
                    if (nNodesInClusterC == 1) { //singleton
                        averageInnerDistance = 1.0;
                    }  else {
                        distance = 0.0;
                        for (i = 0; i < nNodesInClusterC; i++) {
                            k = nodesPerCluster[c][i];
                            if ((j + cumulativeNNodes) != k) {
                                distance += shortestPath[j + cumulativeNNodes][k];
                            }
                        }
                        averageInnerDistance = distance / (nNodesInClusterC - 1);
                    } 
                    if (TEST) {
                        System.out.println("            averageInnerDistance: " + averageInnerDistance);
                    }

                    // calculate the minimum average distance between vertex j 
                    // and the vertices in each cluster c2
                    if (TEST) {
                        System.out.println("        -calculating minimum outer distance");
                    }
                    averageOuterDistances = new double[nClusters];
                    for (c2 = 0; c2 < nClusters; c2++) {
                        nNodesInClusterC2 = nNodesPerCluster[c2];
                        if (c2 == c || nNodesInClusterC2 == 0) {
                            averageOuterDistances[c2] = INF;
                        } else {
                            distance = 0;
                            for (i = 0; i < nNodesInClusterC2; i++) {
                                k = nodesPerCluster[c2][i];
                                distance += shortestPath[j + cumulativeNNodes][k];

                            }
                            averageOuterDistances[c2] = distance / (nNodesInClusterC2);

                        }
                    }

                    minAverageOuterDistance = Arrays2.calcMinimum2(averageOuterDistances);
                    if (TEST) {

                        System.out.println();
                        System.out.println("            minAverageOuterDistance: " + minAverageOuterDistance);
                    }

                    // calculate the silhouette width value
                    silhouetteWidths[j] = (
                        (minAverageOuterDistance - averageInnerDistance) / 
                        (Math.max(averageInnerDistance, minAverageOuterDistance))); 

                    if (TEST) {
                        System.out.println("            silhouetteWidths[j]: " + silhouetteWidths[j]);
                    }
                }
                silhouetteIndexPerCluster[c] = Arrays2.calcSum(silhouetteWidths) / nNodesInClusterC;
            }

            if (TEST) {
                System.out.println("2. calculate silhouette of entire cluster " + c);
                System.out.println("    silhouetteIndexPerCluster[c]: " + silhouetteIndexPerCluster[c]);
            }

            cumulativeNNodes += nNodesInClusterC;
        }

        // get the silhouette of entire graph
        silhouetteIndexGlobal = Arrays2.calcSum(silhouetteIndexPerCluster) / (nClusters - nEmptyClusters);

        if (TEST) {
            System.out.println();
            System.out.println("silhouetteIndexGlobal: " + silhouetteIndexGlobal);
        }
        return silhouetteIndexGlobal;
    }

    /**
     * Find the partition of nodes that maximizes the silhouette index score.
     *
     * Note that while this method only returns a boolean, its work is saved in 
     * the instance variable clustering.
     *
     * @param random  a random num generator 
     * @return true if nodes are moved to a different community; false otherwise 
     */
    public boolean runLocalMovingAlgorithmWithSilhouetteIndex(Random random) {
        boolean update;
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, edgeWeightPerCluster;
        double[][] shortestPath;
        int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters, originalCluster;
        int[] neighboringCluster, newCluster, nNodesPerCluster, nodePermutation, unusedCluster;
 
        update = false;

        if (network.nNodes == 1)
            return update;


        // store the total degree of and the number of nodes in each cluster
        clusterWeight = new double[network.nNodes];
        nNodesPerCluster = new int[network.nNodes];
        for (i = 0; i < network.nNodes; i++) {
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
            nNodesPerCluster[clustering.cluster[i]]++;
        }

        // store the empty clusters
        nUnusedClusters = 0;
        unusedCluster = new int[network.nNodes];
        for (i = 0; i < network.nNodes; i++) {
            if (nNodesPerCluster[i] == 0) {
                unusedCluster[nUnusedClusters] = i;
                nUnusedClusters++;
            }
        }

        nodePermutation = Arrays2.generateRandomPermutation(
            network.nNodes, random);
        edgeWeightPerCluster = new double[network.nNodes];
        neighboringCluster = new int[network.nNodes - 1];
        nStableNodes = 0;
        i = 0;
        shortestPath = floydWarshall(network.getMatrix(), network.nNodes);

        // compute the clustering that maximizes the silhouette index value
        int nDoWhileIterations = 0; 
        do {
            j = nodePermutation[i]; // start with a random vertex j
            originalCluster = clustering.cluster[j];
           
            if (TEST2) {
                System.out.println();
                System.out.println("------------------------------------------");
                System.out.println("ITERATION " + nDoWhileIterations + "; nStableNodes is " + nStableNodes + "; current vertex is " + j + " in cluster " + originalCluster);
                nDoWhileIterations++;
            }

            // get the neighboring clusters of vertex j and their edge weights
            nNeighboringClusters = 0;
            for (k = network.firstNeighborIndex[j]; k < network.firstNeighborIndex[j + 1]; k++) {
                l = clustering.cluster[network.neighbor[k]];            
                if (edgeWeightPerCluster[l] == 0) {
                    neighboringCluster[nNeighboringClusters] = l;
                    nNeighboringClusters++;
                }
                edgeWeightPerCluster[l] += network.edgeWeight[k];
            }

            // remove vertex j from its cluster and update structures
            clusterWeight[clustering.cluster[j]] -= network.nodeWeight[j]; 
            nNodesPerCluster[clustering.cluster[j]]--;
            if (nNodesPerCluster[clustering.cluster[j]] == 0) {
                unusedCluster[nUnusedClusters] = clustering.cluster[j];
                nUnusedClusters++;
            }

            // prepare to simulate adding j to each neighboring cluster
            bestCluster = -1;
            maxQualityFunction = calcSilhouetteFunction(shortestPath);
            if (TEST2) {
                System.out.println("    -initial SI: " + maxQualityFunction);
            }

            // determine the cluster that vertex j should be added to
            for (k = 0; k < nNeighboringClusters; k++) {

                // move vertex j into cluster l and calculate SI score
                l = neighboringCluster[k]; 
                clustering.setCluster(j, l);
                qualityFunction = calcSilhouetteFunction(shortestPath);

                if (TEST2) {
                    System.out.println("    -adding vertex " + j + " to cluster " + l + ": newSI = " + qualityFunction);
                }

                // adding vertex j to neighboring cluster l gives a 
                // better SI score
                // TODO tie-breaking step here differs from modularity version
                if ( (qualityFunction > maxQualityFunction) || ( (qualityFunction == maxQualityFunction) && ((l < bestCluster) || (bestCluster == -1)) ) ) {
                    if (TEST2) {
                        System.out.println("    -bestCluster updated: from " + bestCluster + " to " + l);
                        System.out.println("    -maxQualityFunction updated from " + maxQualityFunction + " to " + qualityFunction);
                    }
                    bestCluster = l;
                    maxQualityFunction = qualityFunction;
                }
                edgeWeightPerCluster[l] = 0;
            }
            // the original cluster of vertex j gives a better SI score
            if (bestCluster == -1) {
                if (TEST2) {
                    System.out.println("    -bestCluster is original: from " + bestCluster + " to " + unusedCluster[nUnusedClusters - 1] + "maxQualityFunction not updated: " + maxQualityFunction);
                }

                bestCluster = unusedCluster[nUnusedClusters - 1];
                nUnusedClusters--;
            } 

            // add vertex j to the cluster that gives the best SI score and
            // update structures accordingly
            clusterWeight[bestCluster] += network.nodeWeight[j];
            nNodesPerCluster[bestCluster]++;
            if (bestCluster == originalCluster) {
                nStableNodes++;
                if (TEST2) {
                    System.out.println("    -bestCluster is original cluster " + bestCluster);
                    System.out.println("    -nStableNodes increments 1: " + nStableNodes);
                }
            } else { 
                clustering.cluster[j] = bestCluster;
                nStableNodes = 1;
                update = true;
                if (TEST2) {
                    System.out.println("    -bestCluster is neighboring cluster " + bestCluster);
                    System.out.println("    -nStableNodes resets to 1: " + nStableNodes);
                }
            }
            clustering.setCluster(j, bestCluster);

            i = (i < network.nNodes - 1) ? (i + 1) : 0;
        }
        while (nStableNodes < network.nNodes);

        // update structures to reflect the new clustering
        newCluster = new int[network.nNodes];
        clustering.nClusters = 0;
        for (i = 0; i < network.nNodes; i++) {
            if (nNodesPerCluster[i] > 0) {
                newCluster[i] = clustering.nClusters;
                clustering.nClusters++;
            }
        }
        for (i = 0; i < network.nNodes; i++) {
            clustering.cluster[i] = newCluster[clustering.cluster[i]];
        }

        return update;
    }

    /***********************************************************************
     * COVERAGE
     ***********************************************************************/

    public double calcCoverageFunction() {
        double qualityFunction;
        double[] clusterWeight;
        int i, j, k;

        // calculate the total edge weight of the network
        qualityFunction = 0;
        for (i = 0; i < network.nNodes; i++) {
            j = clustering.cluster[i];
            for (k = network.firstNeighborIndex[i]; k < network.firstNeighborIndex[i + 1]; k++)
                if (clustering.cluster[network.neighbor[k]] == j) 
                    qualityFunction += network.edgeWeight[k];
        }
        qualityFunction += network.totalEdgeWeightSelfLinks;

        // store the total node weight of each cluster
        clusterWeight = new double[clustering.nClusters];
        for (i = 0; i < network.nNodes; i++) {
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
        }
   
        qualityFunction /= 2 * network.getTotalEdgeWeight() + network.totalEdgeWeightSelfLinks;

        return qualityFunction;
    }

}