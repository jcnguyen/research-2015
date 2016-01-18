/**
 * VOSClusteringTechnique
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/23/14
 */

import java.util.Random;
import java.util.Arrays;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.text.SimpleDateFormat;

public class VOSClusteringTechnique {
    // to denote nonexistent edges in the adjacency matrix
    private static final double INF = Double.POSITIVE_INFINITY;

    // keep track of what recursive level we're on
    public static int level = 0;

    // where to print to
    private static String fileName;

    // the network and its communities
    protected Network network;
    protected Clustering clustering;
    protected double resolution;

    /* FOR PERFORMANCE */

    // adjacency matrix representation of current network
    protected double[][] adjMatrix;

    // meaningful maximum of edge weights
    private static double M;

    /** 
    * v_scalingParam is the scaling parameter that rates the importance of the
    * weight of intercluster edges (with respect to the weight of the
    * intra-cluster edges)
    * 0 <= v_scalingParam <= 1 
    */
    private static double v_scalingParam;

    /*
    * Number of possible edges, including self loops, multiplied by M
    * Denominator in performance calculation;
    * constant for any given network regardless of local shifts
    */
    private double denominator; 

    /*
    * Keep track of the clusterings at each level, in order to later
    * calculate the level that performance should have stopped at.
    */
    private static ArrayList<Clustering> levelClusterings; 
    private static boolean singleNode;
    static {
        levelClusterings = new ArrayList<Clustering>();
    }

    /**
     * This constructor is called when every vertex is in its own cluster.
     *
     * @param network - object representation of graph
     * @param resolution - granularity level at which communities are detected, 1.0 for standard modularity-based community detection
     */
    public VOSClusteringTechnique(double v_scalingParam, double maxM, Network network, double resolution, String fileName) {

        // on initialization, every node is in its own community
        clustering = new Clustering(network.nNodes);
        clustering.initSingletonClusters();

        this.resolution = resolution;
        this.network = network;
        this.fileName = fileName;
        initPerfVariables(v_scalingParam, maxM);
    }

    /**
     * 
     * @param network - object representation of graph
     * @param clustering - object representation of the communities in graph
     * @param resolution - granularity level at which communities are detected
     */
    public VOSClusteringTechnique(double v_scalingParam, double maxM, Network network, Clustering clustering, 
                                    double resolution, String fileName) {
        this.network = network;
        this.clustering = clustering;
        this.resolution = resolution;
        this.fileName = fileName;
        initPerfVariables(v_scalingParam, maxM);
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
    public boolean runLocalMovingAlgorithm(Random random) {
        boolean update;
        double maxQualityFunction, qualityFunction;
        double[] clusterWeight, edgeWeightPerCluster;
        int bestCluster, i, j, k, l, nNeighboringClusters, nStableNodes, nUnusedClusters;
        int[] neighboringCluster, newCluster, nNodesPerCluster, nodePermutation, unusedCluster;

        /*don't need to run alg if only 1 node*/ 
        if (network.nNodes == 1) {
            return false;
        }

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

        /*stay in loop until a local optimal performance value is found, 
        and moving any node would not increase it*/
        do {
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

                // calculate performance change
                qualityFunction = addPerfChange(j, l);

                // if we've found a better cluster, or if we have found an equivalent cluster
                // that ranks higher in our tie breaking algorithm (alphanumeric ordering)
                if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) 
                                                            && (l < bestCluster))) {
                    bestCluster = l;
                    maxQualityFunction = qualityFunction;
                }

                edgeWeightPerCluster[l] = 0;
            }


            // if we move ourselves into another community, our previous community becomes unused
            if (maxQualityFunction == 0) {
                bestCluster = unusedCluster[nUnusedClusters - 1];
                nUnusedClusters--;
            }

            /*add j into the best cluster it should be in*/
            clusterWeight[bestCluster] += network.nodeWeight[j];
            nNodesPerCluster[bestCluster]++;

            /*if it ends up in original cluster, update stable nodes*/
            if (bestCluster == clustering.cluster[j]) {
                nStableNodes++;
            } else {
                
                /*j was moved to a new cluster that is better*/
                clustering.cluster[j] = bestCluster;
                nStableNodes = 1;
                update = true;
            }

            // iterate through all the nodes in our randomly determined order
            // if we reach the end, start over again at the beginning (node 0)
            i = (i < network.nNodes - 1) ? (i + 1) : 0;

        } while (nStableNodes < network.nNodes);

        /*update number of clusters that exist now, and
        what nodes are in what cluster*/
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


    /**
     * @return running Louvain Algorithm with new random num generator as param
     */
    public boolean runLouvainAlgorithm()
    {
        return runLouvainAlgorithm(new Random());
    }

    /**
     * Runs the original Louvain algorithm
     *
     * @param random - a random num generator 
     * @return whether or not we updated what nodes are in what communties 
     */
    public boolean runLouvainAlgorithm(Random random) {
        boolean update, update2;
        VOSClusteringTechnique new_VOSClusteringTechnique;
        
        /* Print info */
        System.out.println("\n\tPass " + level);
        level++;
        System.out.printf("\tnetwork size: %d nodes, %d edges\n", network.getNNodes(), network.getNEdges());
        System.out.println("\tstart computation: " + new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date()));

        /*no update if only one node*/
        if (network.nNodes == 1) {
            System.out.println("\tend computation (1 node): " + new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date()));

            // remember we hit a single node, for printing to .tree and .perfgraph
            singleNode = true;
            return false;
        }

        /* Phase 1: see if moving any nodes will increase modularity */
        update = runLocalMovingAlgorithm(random);

        System.out.println("\tend computation: " + new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date()));
        
        /* Remember the current clustering */
        levelClusterings.add((Clustering)clustering.clone());

        /* if we ended up moving any nodes into different communities, 
        begin recursive procedure. same as asking: if (update) {} */
        if (clustering.nClusters < network.nNodes) {

            // Phase 2: recursively call the next pass of the algorithm
            new_VOSClusteringTechnique = new VOSClusteringTechnique(v_scalingParam, M, 
                                            network.createReducedNetwork(clustering), resolution, fileName);

            /*run Louvain again to see if another update*/
            update2 = new_VOSClusteringTechnique.runLouvainAlgorithm(random);

            if (update2) {
                update = true;

                // Merge results from the next level deep recursive call into this call
                clustering.mergeClusters(new_VOSClusteringTechnique.clustering);

            }
        }

        /* keep track of which iteration we're on */
        if (level > 1) { level--; }

        /* On our lowest level iteration (base case) */
        if (level == 1) {

            /* Highest level clustering is currently saved. We couldn't add it before because of
            the recursion's ordering. Add this clustering to the end of the list. */
            if (singleNode) {
                int[] oneNodeAr = {0};
                levelClusterings.add(new Clustering(oneNodeAr));
            }

            /* Print the communities at each level to a .tree file */
            print_current_communities();

            /* After we've merged all recursive calls and built up the record of clusterings 
            at each level, we calculate the  level at which the algorithm should've stopped 
            to maximize performance. Print the result to a .perfgraph file. */
            bestOverallClustering();
        }

        return update;
    }

    /***********************************************************************
     * PERFORMANCE
     ***********************************************************************/

    /* 
    * Overloaded method - calculating performance function for the current clustering
    * @return - the performance score of the current graph partition
    */
    public double calcPerformanceFunction() {
        return calcPerformanceFunction(this.clustering);
    }

    /* 
    * @param network - network to calculate on
    * @param clust - network's clustering
    * @return - the performance score of a
    */
    public double calcPerformanceFunction(Clustering clust) {
        
        int u, v;        // vertices
        double  f=0,        // intra-cluster density
                g=0;        // inter-cluster sparsity

        // used in calculating g_w
        int numIntercommEdges = 0;
        double sumWeightIntercommEdges = 0;

        // for every vertex 
        for (u = 0; u < network.nNodes; u++) {

            // for all vertices lexicographically equal to or after u (exclude self loops and duplicates)
            for (v = u+1; v < network.nNodes; v++) {

                // u and v are in the same community
                if (clust.cluster[u] == clust.cluster[v]) {

                    // intra-community edge exists, so increase intra-community density score
                    if (adjMatrix[u][v] != INF) {
                        f += adjMatrix[u][v]; 
                    }

                // u and v are in different communities
                } else if (clust.cluster[u] != clust.cluster[v]) {

                    // inter-community edge doesn't exist, so increase inter-community sparsity score
                    if (adjMatrix[u][v] == INF) {
                        g += M;
                    } else {        

                        // edge exists: do calculations for g_w: inter-community edges
                        numIntercommEdges++; 
                        sumWeightIntercommEdges += adjMatrix[u][v];
                    }
                } 

                /* 
                * ELSE: Where u and v are in the same community, but there is no edge (u,v), OR 
                * where u and v are in different communities, but there is an edge (u,v),
                * this is not a "correctly interpreted" edge. So, we don't count it toward our
                * performance score.
                */ 
            }
        }        

        /**
        * g_w is the difference of weight that would be counted if no inter-cluster edges were present,
        * minus the weight that is assigned to the actual inter-cluster edges:
        */
        double g_w = M*numIntercommEdges - sumWeightIntercommEdges;
        double numerator = f + g + v_scalingParam*g_w;
        return numerator / denominator;
    }


    /* 
    * Computes the potential performance change adding vertex,
    * currently a community by itself, into a target comm.
    *
    * Adding a node into tcomm changes its relations with
    * vertices in tcomm. To all other vertices, it is still just
    * some vertex in a different community.
    *
    * We don't look at self loops because self loops always 
    * count toward the f score (correctly interpreted intracommunity 
    * edges), so there is no change due to self loops.
    *
    * @param vertex - node to place
    * @param tcomm - community to put node into
    * @return - change in performance 
    */
    private double addPerfChange(int vertex, int tcomm) {
        assert(vertex>=0 && vertex<network.nNodes);
        double numeratorChange = 0;

        // look at all possible edges, existent or nonexistent, with one end at node
        for (int u=0; u<network.nNodes; u++) {

            // ignore self loops
            // if node in our target community
            if (u != vertex && clustering.cluster[u] == tcomm) {

                // the edge exists
                if (adjMatrix[vertex][u] != INF) {

                    // increase intracluster density score
                    numeratorChange += adjMatrix[vertex][u];

                    // part of g_w: subtract 1, scaled, from the count of inter-community edges
                    numeratorChange -= v_scalingParam*M;

                    // part of g_w: subtract the weight, scaled, from 
                    // the sum total of weight of the intra-cluster edges,
                    // which is itself subtracted from the numerator. so, add.
                    numeratorChange += v_scalingParam*adjMatrix[vertex][u];
    
                } 

                // subtract from g(G) the nonexistent edges (vertex, u),
                // where u is in tcomm. We previously counted this edge in g, but
                // we don't count nonexistent edges within the community toward performance,
                // so we must remove its weight
                else {
                    numeratorChange -= M;
                }
            }
        }

        return numeratorChange/denominator;
    }


    /**
    * Initialize performance variables that are constant for a given network
    */
    private void initPerfVariables(double v_scalingParam, double maxM) {

        // initialize scaling parameter
        this.v_scalingParam = v_scalingParam;

        // get adjacency matrix representation of current network
        adjMatrix = network.getMatrix();

        /* Initialize the meaningful maximum M of edge weights for this pass
        of the Louvain algorithm.

        If we're in the initial call of the algorithm, we use the M passed
        in from the command line. If we're in later calls of the algorithm,
        there may be some edges (ignoring self loops) that are greater than M, 
        because the weights of edges (u,v) in the condensed community are the 
        sum of inter-community edges between those communities represented 
        by u and v. */

        if (level == 0) {
            this.M = maxM;

            /* if any weights are higher than M, set them equal to M*/
            network.adjustEdgeWeights(M); 
        } else {

            /* if we need to set a higher M, setM will return that higher M. 
            otherwise, this statement is equivalent to this.M = maxM */
            this.M = setMtoMax(M);
        }

        // calculate denominator (number of possible edges * meaningful maximum M)
        denominator = ((double)network.nPossibleEdges()) * M;
    }

    /*
    * Given the original network and the iterated clusterings at each level,
    * finds the clustering that gives the best performance. Prints it to a
    * .perfgraph file.
    *
    */
    public void bestOverallClustering() {
        System.out.println("\nCalculating the best overall clustering for all passes...");
        int nLevels = levelClusterings.size();
        double[] levelPerformances = new double[nLevels];
        Clustering nextLowerClustering;

        /* the first clustering is the lowest level, so there's nothing to merge,
        but we still calculate performance */
        levelPerformances[0] = calcPerformanceFunction(levelClusterings.get(0));

        /* at every subsequent level, calculate the merged network and performance */
        for (int i = 1; i < nLevels; i++) {

            /* copy of clustering at next lowest (unfurled) level; need copy so we 
            don't overwrite that clustering */
            nextLowerClustering = (Clustering)levelClusterings.get(i-1).clone();

            /* merge to next lowest level, and save the result */
            nextLowerClustering.mergeClusters(levelClusterings.get(i));
            levelClusterings.set(i, nextLowerClustering);

            /* calculate performance for merged and saved clustering */
            levelPerformances[i] = calcPerformanceFunction(levelClusterings.get(i));
        }

        /* find maximum performance and corresponding clustering */
        int indexOfBestClustering = Arrays2.calcMaximumIndex(levelPerformances);
        Clustering bestClustering = levelClusterings.get(indexOfBestClustering);

        System.out.format("Best clustering is Pass %d with performance: %.4f%n", 
                    indexOfBestClustering, levelPerformances[indexOfBestClustering]);

        /* print .perfgraph file. does not append */
        String perfgraph_fileName = fileName + ".perfgraph";
        try {

            System.out.println("Writing to " + perfgraph_fileName + "...");
            BufferedWriter bufferedWriter;
            bufferedWriter = new BufferedWriter(new FileWriter(perfgraph_fileName));

            /* print level and perf calculation */
            bufferedWriter.write("Level " + indexOfBestClustering + " performance: " 
                                + levelPerformances[indexOfBestClustering]);

            /* print clustering */
            for (int i = 0; i < bestClustering.getNNodes(); i++) {
                bufferedWriter.newLine();
                bufferedWriter.write(i + " " + Integer.toString(bestClustering.getCluster(i)));
            } 

            bufferedWriter.close();
            System.out.println("Finished writing to " + perfgraph_fileName + ".");


        } catch (IOException e) {
            System.out.println("Error printing to .perfgraph file: " + e.getMessage());
        }
    }


    /**
    * Prints the adjacency matrix
    */
    private void printAdjMatrix() {

        for(int u=0; u < adjMatrix.length; u++) {
            for(int v=0; v < adjMatrix.length; v++) {
                System.out.print(adjMatrix[u][v]);
            }
            System.out.println();
        }
    }

    /**
    * Prints to .tree file the communities at each level
    * Does not append.
    */
    private void print_current_communities() {
        Clustering curClustering;
        int nNodes;
        int nLevels = levelClusterings.size();
        String tree_fileName = fileName + ".tree";  

        try {

            System.out.println("\nWriting to " + tree_fileName + "...");
            BufferedWriter bufferedWriter;
            bufferedWriter = new BufferedWriter(new FileWriter(tree_fileName));

            /* for every level */
            for (int l = 0; l < nLevels; l++) {

                /* get the current clustering and its size */
                curClustering = levelClusterings.get(l);
                nNodes = curClustering.getNNodes();

                /* write every node in that clustering to the file */
                for (int i = 0; i < nNodes; i++) {
                    bufferedWriter.write(i + " " + Integer.toString(curClustering.getCluster(i)));
                    bufferedWriter.newLine();
                } 
            }

            bufferedWriter.close();
            System.out.println("Finished writing to " + tree_fileName + ".");

        } catch (IOException e) {
            System.out.println("Error printing to .tree file: " + e.getMessage());
        }

    }

    /*
    * Method called by initPerfVariables().
    * If any edge weights (not self loops) are higher than M,
    * set M to the highest edge weight.
    *
    * @pre - adjacency matrix already generated
    * @post - M updated, adjMatrix unchanged
    */
    private double setMtoMax(double origM) {
        double newM = origM;

        // All edges except self-loops, no repeats
        for (int u=0; u<network.nNodes; u++) {
            for (int v=u+1; v<network.nNodes; v++) {

                // keep track of max edge weight
                if (adjMatrix[u][v] > newM) {
                    newM = adjMatrix[u][v];
                }
            }
        }

        if (newM > origM) {
            System.out.println("\tUpdated M from " + origM + " to " + newM);
        }

        return newM;
    }

}
