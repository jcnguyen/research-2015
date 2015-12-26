/**
 * VOSClusteringTechnique
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/23/14
 */

/*
 * Floyd-Warshall algorithm for APSP problem
 * contributed by Aakash Hasija at geeksforgeeks.org
 */

import java.util.Random;
import java.util.Arrays;

public class VOSClusteringTechnique {
    private static final boolean TEST = false;
    private static final boolean TEST2 = true;
    private static final boolean TEST3 = false;
    private static final double INF = Double.POSITIVE_INFINITY;

    // the network and its communities
    protected Network network;
    protected Clustering clustering;
    protected double resolution;

    /**
     * This constructor is called where every vertex is in its own cluster.
     *
     * @param network - object representation of graph
     * @param resolution - granularity level at which communities are detected, 1.0 for standard modularity-based community detection
     */
    public VOSClusteringTechnique(Network network, double resolution) {
        this.network = network;

        // on initialization, every node is in its own community
        clustering = new Clustering(network.nNodes);
        clustering.initSingletonClusters();

        this.resolution = resolution;
    }

    /**
     * 
     * @param network - object representation of graph
     * @param clustering - object representation of the communities in graph
     * @param resolution - granularity level at which communities are detected
     */
    public VOSClusteringTechnique(Network network, Clustering clustering, double resolution) {
        this.network = network;
        this.clustering = clustering;
        this.resolution = resolution;
    }

    /**
     * 
     * @return network - object representation of a graph
     */
    public Network getNetwork() {
        return network;
    }

    /**
     * 
     * @return clustering - object representation of the communities in graph
     */
    public Clustering getClustering() {
        return clustering;
    }

    /**
     * 
     * @return resolution - granularity level at which communities are detected
     */
    public double getResolution() {
        return resolution;
    }

    /**
     * 
     * @param network - object representation of graph
     */
    public void setNetwork(Network network) {
        this.network = network;
    }

    /**
     * 
     * @param clustering - object representation of the communities in graph 
     */
    public void setClustering(Clustering clustering) {
        this.clustering = clustering;
    }
    
    /**
     * 
     * @param resolution - granularity level at which communities are detected 
     */
    public void setResolution(double resolution) {
        this.resolution = resolution;
    }

    /**
     * 
     * @return runs the local moving algorithm, with new random num generator as input
     */
    public boolean runLocalMovingAlgorithm() {
        return runLocalMovingAlgorithm(new Random());
    }

    /**
     * 
     * Finds what partition of nodes maximizes the modularity, puts nodes
     * in those clusters.
     *
     * Note that while this method only returns a boolean, its work is saved in 
     * the instance variable clustering.
     *
     * @param random - a random num generator 
     * @return whether or not we updated what nodes are in what communties 
     */
    public boolean runLocalMovingAlgorithm(Random random) {
        boolean update;
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, edgeWeightPerCluster;
        int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
        int[] neighboringCluster, newCluster, nNodesPerCluster, nodePermutation, unusedCluster;

        if (network.nNodes == 1)
            return false;

        update = false;

        clusterWeight = new double[network.nNodes]; /*elements contain the total node weights of everything in that cluster*/
        nNodesPerCluster = new int[network.nNodes]; /*elements contain num nodes in that cluster*/
        for (i = 0; i < network.nNodes; i++) {
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
            nNodesPerCluster[clustering.cluster[i]]++;
        }

        /* an unusedCluster is an array that keeps track of clusters with no nodes in them */
        nUnusedClusters = 0;
        unusedCluster = new int[network.nNodes];
        for (i = 0; i < network.nNodes; i++)
            if (nNodesPerCluster[i] == 0) {
                unusedCluster[nUnusedClusters] = i;
                nUnusedClusters++;
            }

        // array of nodes in a random order
        nodePermutation = Arrays2.generateRandomPermutation(network.nNodes, random);

        edgeWeightPerCluster = new double[network.nNodes];
        neighboringCluster = new int[network.nNodes - 1];
        nStableNodes = 0;
        i = 0;

        /*stay in loop until a local optimal modularity value is moved, 
        moving any node would not increase it*/
        do {
            j = nodePermutation[i]; /*start with some node*/

            nNeighboringClusters = 0;

            /* for each neighboring node of j, find the cluster of that node, 
            find the number of neighboring clusters, find total edge weight for each neighbor cluster */
            for (k = network.firstNeighborIndex[j]; k < network.firstNeighborIndex[j + 1]; k++) {
                l = clustering.cluster[network.neighbor[k]];
                if (edgeWeightPerCluster[l] == 0) {
                    neighboringCluster[nNeighboringClusters] = l;
                    nNeighboringClusters++;
                }
                edgeWeightPerCluster[l] += network.edgeWeight[k];
            }

            /*remove j from its cluster. update number nodes in cluster and
            if the cluster became empty, add it to unusedClusters */
            clusterWeight[clustering.cluster[j]] -= network.nodeWeight[j]; 
            nNodesPerCluster[clustering.cluster[j]]--;
            if (nNodesPerCluster[clustering.cluster[j]] == 0) {
                unusedCluster[nUnusedClusters] = clustering.cluster[j];
                nUnusedClusters++;
            }

            // prepare to simulate adding j to each neighboring cluster
            bestCluster = -1;
            maxQualityFunction = 0;

            /*simulate adding j to each neighboring cluster, 
            see if any of these give a better quality function, 
            find one that gives the best quality function that j should go in*/
            for (k = 0; k < nNeighboringClusters; k++) {

                // target cluster
                l = neighboringCluster[k];

                // calculate modularity change
                qualityFunction = edgeWeightPerCluster[l] - network.nodeWeight[j] * clusterWeight[l] * resolution;
                
                // if we've found a better cluster, or if we have found an equivalent cluster
                // that ranks higher in our tie breaking algorithm (alphanumeric ordering)
                if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) && (l < bestCluster))) {

                    bestCluster = l;
                    maxQualityFunction = qualityFunction;
                }

                // TODO: why?
                edgeWeightPerCluster[l] = 0;
            }

            // TODO: if no improvement, why do we mess with unused clusters?
            if (maxQualityFunction == 0) {
                bestCluster = unusedCluster[nUnusedClusters - 1];
                nUnusedClusters--;
            }

            /*add j into the best cluster it should be in*/
            clusterWeight[bestCluster] += network.nodeWeight[j];
            nNodesPerCluster[bestCluster]++;

            /*if it ends up in original cluster, update stable nodes*/
            if (bestCluster == clustering.cluster[j]) 
                nStableNodes++;
            else { /*j was moved to a new cluster that is better*/
                clustering.cluster[j] = bestCluster;
                nStableNodes = 1;
                update = true;
            }

            // iterate through all the nodes in our randomly determined order
            // if we reach the end, start over again at the beginning (node 0)
            i = (i < network.nNodes - 1) ? (i + 1) : 0;

        } 
        while (nStableNodes < network.nNodes);

        /*update nunmber of clusters that exist now, and
        what nodes are in what cluster*/
        newCluster = new int[network.nNodes];
        clustering.nClusters = 0;
        for (i = 0; i < network.nNodes; i++)
            if (nNodesPerCluster[i] > 0) {
                newCluster[i] = clustering.nClusters;
                clustering.nClusters++;
            }
        for (i = 0; i < network.nNodes; i++) {
            clustering.cluster[i] = newCluster[clustering.cluster[i]];
        }

        return update;
    }

/*******************************************************************
    IN PROGRESS
********************************************************************/

    public boolean runLocalMovingAlgorithm2(Random random, int modularityFunction) {
        boolean update;
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, edgeWeightPerCluster;
        double[][] shortestPath;
        int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
        int[] neighboringCluster, newCluster, nNodesPerCluster, nodePermutation, unusedCluster;

        double maxSI, originalSI;

        update = false;
 
        if (network.nNodes == 1)
            return update;

        // store the total degree of each cluster and 
        // the number of nodes in each cluster
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

        // randomizes the node ordering
        nodePermutation = Arrays2.generateRandomPermutation(network.nNodes, random);

        // initialize structures
        edgeWeightPerCluster = new double[network.nNodes];
        neighboringCluster = new int[network.nNodes - 1];
        nStableNodes = 0;
        i = 0;

        // compute the shortest distances between every pair of nodes
        shortestPath = floydWarshall(network.getMatrix(), network.nNodes);

        // compute the clustering that achieves the optimal metric value
        boolean bestClusterIsOriginal = false;
        int nDoWhileIterations = 0;
        do {

            j = nodePermutation[i]; // start with a random vertex j
           
            if (TEST2) {
                System.out.println();
                System.out.println("------------------------------------------");
                System.out.println("ITERATION " + nDoWhileIterations + "; nStableNodes is " + nStableNodes + "; current vertex is " + j);
                nDoWhileIterations++;
            }

            // calculate the SI given this current clustering
            maxSI = calcSilhouetteFunction(shortestPath);
            originalSI = maxSI;
            if (TEST2) {
                System.out.println("    -initial SI: " + maxSI);
            }

            // get the neighboring clusters of vertex j and the edge weights of those clusters
            nNeighboringClusters = 0;
            for (k = network.firstNeighborIndex[j]; k < network.firstNeighborIndex[j + 1]; k++) {
                l = clustering.cluster[network.neighbor[k]]; // the neighboring cluster // TODO: this is global clustering ID                
                if (edgeWeightPerCluster[l] == 0) { // l is not yet added to neighboringCluster
                    neighboringCluster[nNeighboringClusters] = l;
                    nNeighboringClusters++;
                }
                edgeWeightPerCluster[l] += network.edgeWeight[k];
            }

            // remove vertex j from its cluster and update structures accordingly
            clusterWeight[clustering.cluster[j]] -= network.nodeWeight[j]; 
            nNodesPerCluster[clustering.cluster[j]]--;
            if (nNodesPerCluster[clustering.cluster[j]] == 0) {
                unusedCluster[nUnusedClusters] = clustering.cluster[j];
                nUnusedClusters++;
            }

            bestCluster = -1;

            // determine the cluster that vertex j should be added to
            // the best cluster is the cluster that gives the highest SI value
            // options are the original cluster and the neighboring clusters of j
            for (k = 0; k < nNeighboringClusters; k++) {

                // move vertex j into cluster l and calculate SI of new clustering
                l = neighboringCluster[k]; 
                clustering.setCluster(j, l);
                qualityFunction = calcSilhouetteFunction(shortestPath);

                if (TEST2) {
                    System.out.println("    -adding vertex " + j + " to cluster " + l);
                    System.out.println("        -qualityFunction: " + qualityFunction);
                }

                // TODO METRIC STUFF HERE - modularize when you've implemented SI
                // if (modularityFunction <= 2) { // modularity
                //     qualityFunction = edgeWeightPerCluster[l] - network.nodeWeight[j] * clusterWeight[l] * resolution; 
                // } else if (modularityFunction == 3) { // silhouette index
                //     qualityFunction = edgeWeightPerCluster[l] - network.nodeWeight[j] * clusterWeight[l] * resolution;
                // } else { // default
                //     qualityFunction = edgeWeightPerCluster[l] - network.nodeWeight[j] * clusterWeight[l] * resolution; 
                // }

                if ((qualityFunction > maxSI) || ((qualityFunction == maxSI) && ( (l < bestCluster) || (bestCluster == -1) ))) { // best cluster is neighboring cluster l
                    if (TEST2) {
                        System.out.println("        -bestCluster updated: from " + bestCluster + " to " + l);
                        System.out.println("        -maxSI updated from " + maxSI + " to " + qualityFunction);
                    }
                    bestCluster = l;
                    maxSI = qualityFunction;
                }
                edgeWeightPerCluster[l] = 0;
            }
            Double aMaxSI = (Double) maxSI;
            if (bestCluster == -1 || aMaxSI.isNaN()) {  // best cluster is original
                if (TEST2) {
                    System.out.println("        -bestCluster is original: from " + bestCluster + " to " + unusedCluster[nUnusedClusters - 1]);
                    System.out.println("        -maxSI not updated: " + maxSI);
                }

                bestCluster = unusedCluster[nUnusedClusters - 1]; // TODO should this just be originalCluster?? - yes, but we'll just use the given code to be consistent
                nUnusedClusters--;
                bestClusterIsOriginal = true;
            } 

            // add vertex j to the cluster that gives the best SI and
            // update structures accordingly
            if (TEST2) {
                System.out.println("    -updating nStableNodes");
            }
            clusterWeight[bestCluster] += network.nodeWeight[j];
            nNodesPerCluster[bestCluster]++;
            if (bestClusterIsOriginal) { // vertex j is in original cluster
                nStableNodes++;
                if (TEST2) {
                    System.out.println("        -bestCluster is original: nStableNodes++: " + nStableNodes);
                }
            } else { // vertex j is in one of its neighboring clusters
                clustering.cluster[j] = bestCluster;
                nStableNodes = 1;
                update = true;
                if (TEST2) {
                    System.out.println("        -bestCluster is neighbor: nStableNodes = 1: " + nStableNodes);
                }
            }
            clustering.setCluster(j, bestCluster);

            i = (i < network.nNodes - 1) ? (i + 1) : 0;
        }
        while (nStableNodes < network.nNodes);

        // update structures to reflect the new clustering
        newCluster = new int[network.nNodes];
        clustering.nClusters = 0;
        for (i = 0; i < network.nNodes; i++)
            if (nNodesPerCluster[i] > 0) {
                newCluster[i] = clustering.nClusters;
                clustering.nClusters++;
            }
        for (i = 0; i < network.nNodes; i++) {
            clustering.cluster[i] = newCluster[clustering.cluster[i]];
        }

        return update;
    }

/*****************************************************************************
    END IN PROGRESS
******************************************************************************/

    /**
     *
     * @return running Louvain Algorithm with new random num generator as param
     */
    public boolean runLouvainAlgorithm(int modularityFunction) {
        return runLouvainAlgorithm(new Random(), modularityFunction);
    }

    /**
     * 
     * @param random - a random num generator 
     * @return whether or not we updated what nodes are in what communties 
     */
    public boolean runLouvainAlgorithm(Random random, int modularityFunction) {
        boolean update, update2;
        
        VOSClusteringTechnique VOSClusteringTechnique;

        /*no update if only one node*/
        if (network.nNodes == 1)
            return false;

        /*see if moving any nodes will increase modularity*/
        // TODO METRIC STUFF HERE
        update = runLocalMovingAlgorithm2(random, modularityFunction);
        // update = runLocalMovingAlgorithm(random);


        if (clustering.nClusters < network.nNodes) {
            VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

            /*run louvain again to see if another update*/
            update2 = VOSClusteringTechnique.runLouvainAlgorithm(random, modularityFunction);

            if (update2) {
                update = true;

                // updates this instance's clustering such that we can run this iterative loop
                clustering.mergeClusters(VOSClusteringTechnique.clustering);
            }
        }

        return update;
    }

    /**
     * 
     * @param maxNInterations - max iterations you want to run Louvain
     * @return run the algorithm 
     */
    public boolean runIteratedLouvainAlgorithm(int maxNIterations, int modularityFunction) {
        return runIteratedLouvainAlgorithm(maxNIterations, new Random(), modularityFunction);
    }

    /**
     * 
     * @param maxNInterations - max iterations you want to run Louvain
     * @param random - random number generator 
     * @return run the algorithm 
     */
    public boolean runIteratedLouvainAlgorithm(int maxNIterations, Random random, int modularityFunction) {
        boolean update;
        int i;

        i = 0;
        do {
            update = runLouvainAlgorithm(random, modularityFunction);
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

        update = runLocalMovingAlgorithm(random);

        if (clustering.nClusters < network.nNodes) {
            VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

            update2 = VOSClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);

            if (update2) {
                update = true;

                clustering.mergeClusters(VOSClusteringTechnique.clustering);

                runLocalMovingAlgorithm(random);
            }
        }

        return update;
    }

    public boolean runIteratedLouvainAlgorithmWithMultilevelRefinement(int maxNIterations) {
        return runIteratedLouvainAlgorithmWithMultilevelRefinement(maxNIterations, new Random());
    }

    public boolean runIteratedLouvainAlgorithmWithMultilevelRefinement(int maxNIterations, Random random) {
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

        update = runLocalMovingAlgorithm(random);

        if (clustering.nClusters < network.nNodes) {
            subnetwork = network.createSubnetworks(clustering);

            nodePerCluster = clustering.getNodesPerCluster();

            clustering.nClusters = 0;
            nNodesPerClusterReducedNetwork = new int[subnetwork.length];
            for (i = 0; i < subnetwork.length; i++) {
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

    public boolean runIteratedSmartLocalMovingAlgorithm(int nIterations, Random random) {
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

        VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

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
     * Calculates the modularity score of the final community graph.
     *
     * Note: only called by ModularityOptimizer to calculate the modularity of 
     * the final clustering/community graph.
     * 
     * @return the modularity score of the community graph
     */
    public double calcModularityFunction() {
        double qualityFunction;
        double[] clusterWeight;
        int i, j, k;

        /*for each node in network, if neighbors are in the same cluster, adds their edge weight to total*/
        qualityFunction = 0;
        for (i = 0; i < network.nNodes; i++) {
            j = clustering.cluster[i];  /*j is the cluster of node i*/
            for (k = network.firstNeighborIndex[i]; k < network.firstNeighborIndex[i + 1]; k++)
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

    /***********************************************************************
     * SILHOUETTE INDEX
     ***********************************************************************/

    /**
     * Calculates the silhouette index score of the final community graph.
     *
     * Note: only called by ModularityOptimizer to calculate the modularity of 
     * the final clustering/community graph.
     * 
     * @return the silhouette index score of the community graph
     */
    public double calcSilhouetteFunction() {
        return calcSilhouetteFunction(floydWarshall(network.getMatrix(), network.nNodes));
    }

    // c = current cluster
    // c2 = other cluster
    // j = node in cluster c
    // i = iterator for other nodes (either in cluster c or cluster c2) 
    // k = node in other cluster at i
    public double calcSilhouetteFunction(double[][] shortestPath) {
        int c, c2, i, j, k;
        int cumulativeNNodes, nClusters, nNodesInClusterC, nNodesInClusterC2, nEmptyClusters;
        int[] nNodesPerCluster;
        int[][] nodesPerCluster;

        double averageInnerDistance, distance, minAverageOuterDistance;
        double[] averageOuterDistances;

        double silhouetteIndexGlobal;
        double[] silhouetteIndexPerCluster;
        double[] silhouetteWidths;

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
                        System.out.print("            averageOuterDistances: ");
                        for (int jj = 0; jj < averageOuterDistances.length; jj++) {
                            System.out.print(averageOuterDistances[jj] + ", ");
                        }
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
     * Calculates the distance of shortest path between all vertices.
     *
     * @param graph   a matrix representing the input network
     * @param nNodes  the number of vertices/nodes in the network
     * @return a matrix of shortest path distances
     */
    public double[][] floydWarshall(double[][] graph, int nNodes) {
        double[][] distances; // the solution matrix
        int i, j, k; // iterators

        distances = new double[nNodes][nNodes];
 
        // initialize the solution matrix to the input graph matrix
        for (i = 0; i < nNodes; i++) {
            for (j = 0; j < nNodes; j++) {
                if (i == j) { // diagonal
                    distances[i][j] = 0;
                } else {
                    distances[i][j] = graph[i][j];
                }
            }
        }

        // calculate the shortest distance between vertices i and j, 
        // given an intermediate vertex k
        for (k = 0; k < nNodes; k++) { // intermediate vertex
            for (i = 0; i < nNodes; i++) { // source vertex
                for (j = 0; j < nNodes; j++) { // destination vertex, given source vertex
                    // if vertex k is on the shortest path from
                    // i to j, then update the value of distances[i][j]
                    if (distances[i][k] + distances[k][j] < distances[i][j])
                        distances[i][j] = distances[i][k] + distances[k][j];
                }
            }
        }
        
        if (TEST3) {
            for (int jj = 0; jj < nNodes; jj++) {
                for (int kk = 0; kk < nNodes; kk++) {
                    System.out.print(distances[jj][kk] + ", ");
                }
                System.out.println();
            }
        }

        return distances;
    }

}
