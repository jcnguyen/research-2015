/**
 * Clustering.java
 *
 * A clustering is the set of all the clusters of a network. 
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1 11/17/14
 */

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Arrays;

public class Clustering implements Cloneable, Serializable {

    private static final long serialVersionUID = 1;

    protected int nNodes;    // total number of nodes
    protected int nClusters; // number of clusters
    protected int[] cluster; // the cluster that each node i belongs to

    /**
     * Loads a clustering off of a file.
     * 
     * @param  fileName                the file that contains the clustering
     * @throws ClassNotFoundException  occurs if a class is unable to load a 
     * @throws IOException             occurs if there's an input/output error 
     * @return a clustering from a file
     **/
    public static Clustering load(String fileName
        ) throws ClassNotFoundException, IOException {
        Clustering clustering;
        ObjectInputStream objectInputStream;

        objectInputStream = new ObjectInputStream(
            new FileInputStream(fileName));

        clustering = (Clustering)objectInputStream.readObject();

        objectInputStream.close();

        return clustering;
    }

    /** 
     * Constructs a clustering using a defined number of nodes.
     * 
     * @param nNodes  the number of nodes in the network
     **/
    public Clustering(int nNodes) {
        this.nNodes = nNodes;
        cluster = new int[nNodes];
        nClusters = 1;
    }

    /** 
     * Constructs a clustering using a defined cluster.
     * 
     * @param cluster  the cluster used to construct this clustering
     **/
    public Clustering(int[] cluster) {
        nNodes = cluster.length;
        this.cluster = (int[])cluster.clone(); // shallow copy
        nClusters = Arrays2.calcMaximum(cluster) + 1;
    }

    /**
     * Creates a shallow copy of this clustering instance. 
     * That is, if the clone clustering is modified, then the modification is 
     * reflected in this clustering. But if this clustering is modified, then
     * the modification is not reflected in the clone. 
     * 
     * @throws CloneNotSupportedException  occurs if the Cloneable interface 
     *                                     is not implemented by this instance
     * @return a shallow copy of this clustering
     **/
    public Object clone() {
        Clustering clonedClustering;

        try {
            clonedClustering = (Clustering)super.clone();
            clonedClustering.cluster = getClusters();
            return clonedClustering;
        }
        catch (CloneNotSupportedException e) {
            return null;
        }
    }

    /**
     * Saves this clustering into a file.
     *
     * @param  fileName     the file where this clustering is saved
     * @throws IOException  occurs if there's an input or output error
     */
    public void save(String fileName) throws IOException {
        ObjectOutputStream objectOutputStream;

        objectOutputStream = new ObjectOutputStream(
            new FileOutputStream(fileName));

        objectOutputStream.writeObject(this);

        objectOutputStream.close();
    }

    /**
     * @return the total number of nodes in this clustering
     */
    public int getNNodes() {
        return nNodes;
    }

    /**
     * @return the number of clusters in this clustering
     */
    public int getNClusters() {
        return nClusters;
    }

    /**
     * @return an array that stores the cluster of each node
     */
    public int[] getClusters() {
        return (int[])cluster.clone();
    }

    /**
     * Gets the cluster that the node is in.
     *
     * @param  node  the node
     * @return the cluster in which node belongs
     **/
    public int getCluster(int node) {
        return cluster[node];
    }

    /**
     * @return an array that stores the number of nodes of each cluster 
     *         (value at index i is the number of nodes in cluster i)
     **/
    public int[] getNNodesPerCluster() {
        int i;
        int[] nNodesPerCluster;

        nNodesPerCluster = new int[nClusters];
        for (i = 0; i < nNodes; i++) {
            nNodesPerCluster[cluster[i]]++;
        }
        return nNodesPerCluster;
    }

    /**
     * @return an array that stores the nodes of each cluster 
     *         (value at index i is the array storing the nodes in cluster i)
     **/
    public int[][] getNodesPerCluster() {
        int i;
        int[] nNodesPerCluster;
        int[][] nodePerCluster;

        nodePerCluster = new int[nClusters][];
        nNodesPerCluster = getNNodesPerCluster();

        // initialize the inner arrays
        for (i = 0; i < nClusters; i++) { 
            nodePerCluster[i] = new int[nNodesPerCluster[i]];
            nNodesPerCluster[i] = 0;
        }
        
        // fill out the inner arrays
        for (i = 0; i < nNodes; i++) { 
            nodePerCluster[cluster[i]][nNodesPerCluster[cluster[i]]] = i;
            nNodesPerCluster[cluster[i]]++;
        }

        return nodePerCluster;
    }

    /**
     * Move the node into a different cluster.
     *
     * @param node     the node to move
     * @param cluster  the cluster where the node will be moved
     **/
    public void setCluster(int node, int cluster) {
        this.cluster[node] = cluster;
        nClusters = Math.max(nClusters, cluster + 1);
    }

    public void setCluster2(int node, int cluster) {
        this.cluster[node] = cluster;
        nClusters--;
    }

    /**
     * Set each node as its own cluster (all clusters are singletons)
     **/
    public void initSingletonClusters() {
        int i;

        for (i = 0; i < nNodes; i++)
            cluster[i] = i;
        nClusters = nNodes;
    }

    /**
     * Order the clusters in increasing order of number of nodes.
     **/
    public void orderClustersByNNodes() {

        /**
         * Class to compare the number of nodes in two clusters. 
         **/
        class ClusterNNodes implements Comparable<ClusterNNodes> {
            public int cluster;
            public int nNodes;

            /**
             * Keep track of the number of nodes in a cluster.
             *
             * @param cluster  the cluster
             * @param nNodes   the number of nodes in this cluster
             **/
            public ClusterNNodes(int cluster, int nNodes) {
                this.cluster = cluster;
                this.nNodes = nNodes;
            }

            /**
             * Compare the number of nodes of this cluster to another cluster
             *
             * @param  clusterNNodes  the cluster to compare to against
             * @return 1 if compared cluster has more nodes than this cluster
             *        -1 if compared cluster has less nodes than this cluster
             *         0 if the number of nodes in both clusters are equal
             **/
            public int compareTo(ClusterNNodes clusterNNodes) {
                return (clusterNNodes.nNodes > nNodes) ? 1 : (
                    (clusterNNodes.nNodes < nNodes) ? -1 : 0);
            }
        }

        ClusterNNodes[] clusterNNodes;
        int i;
        int[] newCluster, nNodesPerCluster;

        // initialize the array to reflect the number of nodes of each cluster
        nNodesPerCluster = getNNodesPerCluster();
        clusterNNodes = new ClusterNNodes[nClusters];
        for (i = 0; i < nClusters; i++)
            clusterNNodes[i] = new ClusterNNodes(i, nNodesPerCluster[i]);

        // sort the clusters based on the number of nodes
        Arrays.sort(clusterNNodes);
        newCluster = new int[nClusters];
        i = 0;
        do {
            newCluster[clusterNNodes[i].cluster] = i;
            i++;
        }
        while ((i < nClusters) && (clusterNNodes[i].nNodes > 0));

        // update the clusters in this clustering
        nClusters = i;
        for (i = 0; i < nNodes; i++)
            cluster[i] = newCluster[cluster[i]];

    }

    /**
     * Update this clustering according to a new, passed-in clustering.
     *
     * @param clustering  the clustering to copy
     **/
    public void mergeClusters(Clustering clustering) {
        int i;

        for (i = 0; i < nNodes; i++)
            cluster[i] = clustering.cluster[cluster[i]];
        nClusters = clustering.nClusters;
    }
}



/**
 A clustering contains all the clusters in the graph at that time.

       0------3---4 
      / \     |   |
     1   2    5---6
   t = 0: each node is in its own cluster
     cluster = |0|1|2|3|4|5|6| // values = cluster ID
                0 1 2 3 4 5 6  // indices = node ID

   t = 1: cluster0 = {0, 1, 2}, cluster1 = {3, 4, 5, 6}
     cluster = |0|0|0|1|1|1|1|
                0 1 2 3 4 5 6    
 **/