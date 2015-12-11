/**
 * VOSClusteringTechnique
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/23/14
 */

/*
 * APSP contributed by Aakash Hasija
 * at geeksforgeeks.org
 */

import java.util.Random;
import java.util.Arrays;

public class VOSClusteringTechnique {
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
    public VOSClusteringTechnique(Network network, double resolution)
    {
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
    public VOSClusteringTechnique(Network network, Clustering clustering, double resolution)
    {
        this.network = network;
        this.clustering = clustering;
        this.resolution = resolution;
    }

    /**
     * 
     * @return network - object representation of a graph
     */
    public Network getNetwork()
    {
        return network;
    }

    /**
     * 
     * @return clustering - object representation of the communities in graph
     */
    public Clustering getClustering()
    {
        return clustering;
    }

    /**
     * 
     * @return resolution - granularity level at which communities are detected
     */
    public double getResolution()
    {
        return resolution;
    }

    /**
     * 
     * @param network - object representation of graph
     */
    public void setNetwork(Network network)
    {
        this.network = network;
    }

    /**
     * 
     * @param clustering - object representation of the communities in graph 
     */
    public void setClustering(Clustering clustering)
    {
        this.clustering = clustering;
    }
    
    /**
     * 
     * @param resolution - granularity level at which communities are detected 
     */
    public void setResolution(double resolution)
    {
        this.resolution = resolution;
    }

    /**
     * 
     * @return runs the local moving algorithm, with new random num generator as input
     */
    public boolean runLocalMovingAlgorithm()
    {
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
    public boolean runLocalMovingAlgorithm(Random random)
    {
        boolean update;
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, edgeWeightPerCluster;
        int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
        int[] neighboringCluster, newCluster, nNodesPerCluster, nodePermutation, unusedCluster;

        /*don't need to run alg if only 1 node*/ 
        if (network.nNodes == 1)
            return false;

        update = false;

        clusterWeight = new double[network.nNodes]; /*elements contain the total node weights of everything in that cluster*/
        nNodesPerCluster = new int[network.nNodes]; /*elements contain num nodes in that cluster*/
        for (i = 0; i < network.nNodes; i++)
        {
            clusterWeight[clustering.cluster[i]] += network.nodeWeight[i];
            nNodesPerCluster[clustering.cluster[i]]++;
        }

        /*keep track of the num clusters with no nodes in them,
         as well as what clusters these are*/
        nUnusedClusters = 0;
        unusedCluster = new int[network.nNodes];
        for (i = 0; i < network.nNodes; i++)
            if (nNodesPerCluster[i] == 0)
            {
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
        do
        {
            j = nodePermutation[i]; /*start with some node*/

            nNeighboringClusters = 0;

            /* for each neighboring node of j, find the cluster of that node, 
            find the number of neighboring clusters, find total edge weight for each neighbor cluster */
            for (k = network.firstNeighborIndex[j]; k < network.firstNeighborIndex[j + 1]; k++)
            {
                l = clustering.cluster[network.neighbor[k]];
                if (edgeWeightPerCluster[l] == 0)
                {
                    neighboringCluster[nNeighboringClusters] = l;
                    nNeighboringClusters++;
                }
                edgeWeightPerCluster[l] += network.edgeWeight[k];
            }

            /*remove j from its cluster. update number nodes in cluster and
            if the cluster became empty, add it to unusedClusters */
            clusterWeight[clustering.cluster[j]] -= network.nodeWeight[j]; 
            nNodesPerCluster[clustering.cluster[j]]--;
            if (nNodesPerCluster[clustering.cluster[j]] == 0)
            {
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
            if (maxQualityFunction == 0)
            {
                bestCluster = unusedCluster[nUnusedClusters - 1];
                nUnusedClusters--;
            }

            /*add j into the best cluster it should be in*/
            clusterWeight[bestCluster] += network.nodeWeight[j];
            nNodesPerCluster[bestCluster]++;

            /*if it ends up in original cluster, update stable nodes*/
            if (bestCluster == clustering.cluster[j]) 
                nStableNodes++;
            else /*j was moved to a new cluster that is better*/
            {
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
            if (nNodesPerCluster[i] > 0)
            {
                newCluster[i] = clustering.nClusters;
                clustering.nClusters++;
            }
        for (i = 0; i < network.nNodes; i++) {
            clustering.cluster[i] = newCluster[clustering.cluster[i]];
        }

        return update;
    }

    /**
     *
     * @return running Louvain Algorithm with new random num generator as param
     */
    public boolean runLouvainAlgorithm()
    {
        return runLouvainAlgorithm(new Random());
    }

    /**
     * 
     * @param random - a random num generator 
     * @return whether or not we updated what nodes are in what communties 
     */
    public boolean runLouvainAlgorithm(Random random)
    {
        boolean update, update2;
        
        VOSClusteringTechnique VOSClusteringTechnique;

        /*no update if only one node*/
        if (network.nNodes == 1)
            return false;

        /*see if moving any nodes will increase modularity*/
        update = runLocalMovingAlgorithm(random);

        if (clustering.nClusters < network.nNodes)
        {
            VOSClusteringTechnique = new VOSClusteringTechnique(network.createReducedNetwork(clustering), resolution);

            /*run louvain again to see if another update*/
            update2 = VOSClusteringTechnique.runLouvainAlgorithm(random);

            if (update2)
            {
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
    public boolean runIteratedLouvainAlgorithm(int maxNIterations)
    {
        return runIteratedLouvainAlgorithm(maxNIterations, new Random());
    }

    /**
     * 
     * @param maxNInterations - max iterations you want to run Louvain
     * @param random - random number generator 
     * @return run the algorithm 
     */
    public boolean runIteratedLouvainAlgorithm(int maxNIterations, Random random)
    {
        boolean update;
        int i;

        i = 0;
        do
        {
            update = runLouvainAlgorithm(random);
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

    public boolean runIteratedSmartLocalMovingAlgorithm(int nIterations)
    {
        return runIteratedSmartLocalMovingAlgorithm(nIterations, new Random());
    }

    public boolean runIteratedSmartLocalMovingAlgorithm(int nIterations, Random random)
    {
        boolean update;
        int i;

        update = false;
        for (i = 0; i < nIterations; i++)
            update |= runSmartLocalMovingAlgorithm(random);
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

    /***********************************************************************
     * MODULARITY
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

    /***********************************************************************
     * SILHOUETTE INDEX
     ***********************************************************************/

    /**
     * 
     * @return qualityFunction - calculating the metric (modularity in this case) for the graph
     */
    // TODO
    // public double calcSilhouetteFunction()
    public double calcQualityFunction()
    {
        int i, j, k, c, n; // iterators
        int nodeK; // other node
        int nNodes = network.nNodes;
        int numCommunities = clustering.nClusters;
        int[] numNodesInCluster = clustering.getNNodesPerCluster();
        int[][] nodesInCluster = clustering.getNodesPerCluster();

        double[][] adjMatrix = network.getMatrix();
        double[][] shortestPath = floydWarshall(adjMatrix, nNodes);

        double averageInnerDistance;
        double currentDistance = 0;
        int sizeOfCommunityN;
        int sizeOfCommunity;
        double minOfAverageOutsideDistance;
        double maxOfDistances;
        double communitySilhouette;
        double globalSilhouetteIndex = 0;

        double[] allCommunitySilhouettes = new double[numCommunities];

        // TODO delete
        for (int jj = 0; jj < numNodesInCluster.length; jj++) {
            System.out.println("index: " + jj + " " + numNodesInCluster[jj] + " ");
        }

        //for every community
        for (c = 0; c < numCommunities; c++) {
            
            sizeOfCommunity = numNodesInCluster[c];

            // stores the silhouette for each vertex in this community
            double[] silhouetteWidth = new double[sizeOfCommunity];

            //for every node in community c
            for(i = 0; i < sizeOfCommunity; i++) {

                averageInnerDistance = 0.0;

                //find avg of shortest distance between v_i and every other vertex in same community
                // System.out.println("sizeOfCommunity: " + sizeOfCommunity); // TODO delete
                if (sizeOfCommunity == 1) { //singleton
                    averageInnerDistance = 1.0;
                }  else {
                    for (k = 0; k < sizeOfCommunity; k++) { // all node in community
                        nodeK = nodesInCluster[c][k];
                        if (i != nodeK) {
                            averageInnerDistance += shortestPath[i][nodeK];
                        }
                    }
                    averageInnerDistance = averageInnerDistance/(sizeOfCommunity-1);
                } 

                // System.out.println("avgInnerDistance: " + silhouetteWidth[i]); // TODO delete

                // finds distance between vertex i and every vertex in different community
                double[] distanceOptions = new double[sizeOfCommunity];
                distanceOptions[c] = 0;
                for (n = 0; n < numCommunities; n++) {
                    currentDistance = 0;

                    // for every other community
                    if (n != c) {
                        sizeOfCommunityN = numNodesInCluster[n];

                        for (k = 0; k < sizeOfCommunityN; k++) {
                            nodeK = nodesInCluster[n][k];
                            currentDistance += shortestPath[i][nodeK];
                        }
                        currentDistance = currentDistance/(sizeOfCommunityN);
                        distanceOptions[n] = currentDistance;
                    }
                }

                // get the min of all the average distances
                Arrays.sort(distanceOptions);
                minOfAverageOutsideDistance = distanceOptions[0];

                // get the max between the average distance within and
                // min average distance
                if (averageInnerDistance > minOfAverageOutsideDistance) {
                    maxOfDistances = averageInnerDistance;
                } else {
                    maxOfDistances = minOfAverageOutsideDistance;
                }

                silhouetteWidth[i] = ((minOfAverageOutsideDistance - averageInnerDistance)/maxOfDistances); 
                // System.out.println("SI at " + i + ": " + silhouetteWidth[i]); TOOD delete
            }

            // calculate the silhouette width for the entire community
            communitySilhouette = 0;
            for (k=0;k<sizeOfCommunity;k++) {
                communitySilhouette += silhouetteWidth[k];
            }

            allCommunitySilhouettes[c] = communitySilhouette/sizeOfCommunity;
        }

        // get the silhouette of entire graph
        for (k=0;k<numCommunities;k++) {
            globalSilhouetteIndex += allCommunitySilhouettes[k];
        }

        return globalSilhouetteIndex/numCommunities;
    }


    public double[][] floydWarshall(double[][] graph, int nNodes) {
        double[][] dist;
        int i, j, k;

        dist = new double[nNodes][nNodes];
 
        /* Initialize the solution matrix same as input graph matrix.
           Or we can say the initial values of shortest distances
           are based on shortest paths considering no intermediate
           vertex. */
        for (i = 0; i < nNodes; i++) {
            for (j = 0; j < nNodes; j++) {
                if (i == j) {
                    dist[i][j] = 0;
                } else {
                    dist[i][j] = graph[i][j];
                }
            }
        }
 
        /* Add all vertices one by one to the set of intermediate
           vertices.
          ---> Before start of a iteration, we have shortest
               distances between all pairs of vertices such that
               the shortest distances consider only the vertices in
               set {0, 1, 2, .. k-1} as intermediate vertices.
          ----> After the end of a iteration, vertex no. k is added
                to the set of intermediate vertices and the set
                becomes {0, 1, 2, .. k} */
        for (k = 0; k < nNodes; k++)
        {
            // Pick all vertices as source one by one
            for (i = 0; i < nNodes; i++)
            {
                // Pick all vertices as destination for the
                // above picked source
                for (j = 0; j < nNodes; j++)
                {
                    // If vertex k is on the shortest path from
                    // i to j, then update the value of dist[i][j]
                    if (dist[i][k] + dist[k][j] < dist[i][j])
                        dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
 
        return dist;
    }

    /***********************************************************************
     * PERFORMANCE
     ***********************************************************************/

    /**
    * Calculates total modularity of the graph
    *
    */
    public double COPYcalcModularityFunction() {
        double qualityFunction;
        double[] clusterWeight;
        int i, j, k;

        qualityFunction = 0;
        /*for each node in network, if neighbors are in the same cluster, adds their edge weight to total*/
        for (i = 0; i < network.nNodes; i++)
        {
            j = clustering.cluster[i];  /*j is the cluster of node i*/
            for (k = network.firstNeighborIndex[i]; k < network.firstNeighborIndex[i + 1]; k++) /*for every neighbor (edge?!) that i has??*/
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

    /* 
    * Compute the performance of the initial partition of G,
    * in which every node is its own cluster. This is a special 
    * case of calculating performance. That is, f(G) = 0 
    * because every node is one clustering and thus there are 
    * no internal edges within communities. 
    *
    * g(G) = num_possible_edges - |E|, because at this point
    * every nonexistent edge (disregarding loops) between 
    * 2 vertices is between two separate communities, so they
    * count in the g(G) score.
    * 
    * If this a weighted graph, we multiply the number of nonexistent 
    * external edges by the average weight of the nonexistent external edges,
    * which we calculate from the other weighted edges.
    *
    * Finally, we divide g(G) by the total number of possible edges
    * between all vertices (num_possible_edges).
    *
    * @return perf - performance score
    */

    private double initPerformanceFunction() {

        // TODO make num_possible_edges an instance variable
        // see if protected status of nNodes and nEdges will be a problem
        // TODO comments about how we're calculating num_possible_edges
        // where does it calculate original modularity?
        int num_possible_edges = network.nNodes*(network.nNodes + 1)/2 - network.nNodes;
        double g = num_possible_edges - network.nEdges;

        return g / (double)num_possible_edges;

    }

    /* 
    * @return - the performance score of the current graph partition
    */
    private calcPerformanceFunction() {
        double perf = 0;

      /* MY PSEUDOCODE:

      for every vertex u in the graph 
        for every vertex v > u (strict inequality to prevent duplicates)
            if u and v are in the same community
              if the edge (u,v) exists,
                f = f + w(u,v)
              else do nothing
            else 
              if the edge between them does not exist, 
                g = g + w (where w is 1 for unweighted graph, 
                          otherwise it's the average weight for nonexistent edges)
              else do nothing

      add f and g
      divide by num_possible_edges
      */

        return perf;
    }

    /*
* @param node - the node to be removed
* @param comm - community containing node
* @param dnodecomm - number of links from node to comm
* 
* @return performance change from removing node from comm
*/
private void remove_perf_change(int node, int comm, double dnodecomm) {
  assert(node>=0 && node<size);
  int change = 0;
  
  // look at all possible edges, existent or nonexistent, with one end at node
  for(int i=0; i<size; i++) {

      // Removing node from comm only changes its relations with
      // vertices in comm. To all other vertices, it is still just
      // some vertex in a different community.
      if (n2c[i] == comm) {

          // since node and vertices v in comm are no longer in the same
          // community, we increase g(G) for nonexistent edges (node,v),
          if (A[node][i] == NO_EDGE) {
              change += 1;    // TODO: change for weighted graphs
          } 

          // removing node from its community deducts from f(G)
          // those edges that went from (node, v), where v is a vertex in comm 
          else {
              change -= A[node][i];
          }

      }
  }

  return (float)change/(float)num_possible_edges;
}

/* 
* computes the potential performance gain
*
* TODO: remove some of these inputs?
* @param node - node to place
* @param comm - community to put node into
* @param dnodecomm - number of links from node to comm
* @param w_degree - node degree
* @return - change in performance 
*/
private double performance_gain(int node, int comm, double dnodecomm, double w_degree) {
  assert(node>=0 && node<size);
  int change = 0;
  
  // look at all possible edges, existent or nonexistent, with one end at node
  for(int i=0; i<size; i++) {

      // Adding a node into comm only changes its relations with
      // vertices in comm. To all other vertices, it is still just
      // some vertex in a different community.
      if (n2c[i] == comm) {

          // subtract from g(G) the nonexistent edges (node, v),
          // where v is in comm. We don't count nonexistent edges
          // within the community toward our performance score
          if (A[node][i] == NO_EDGE) {
              change -= 1;    // TODO: change for weighted graphs
          } 

          // adding node to comm adds to f(G)
          // for existent edges (node, v) where v is a vertex in comm
          else {
              change += A[node][i];
          }
      }
  }

  return (float)change/(float)num_possible_edges;
}



}
