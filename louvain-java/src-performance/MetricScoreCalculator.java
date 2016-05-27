import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Console;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

public class MetricScoreCalculator {
    public static void main(String[] args) throws IOException {

        System.out.println("-------------------------------------");
        System.out.println("Metric Score Calculator");
        System.out.println("-------------------------------------");

        // read in the arguments
        String inputFileName = args[0];
        String clusteringFileName = args[1];
        double mVal = Double.parseDouble(args[2]);

        // read files
        Network network = readInputFile(inputFileName);
        Clustering clustering = readClusteringFile(clusteringFileName);

        // calculate metric scores
        VOSClusteringTechnique VOSClusteringTechnique = new VOSClusteringTechnique(1, mVal, network, clustering, 1, "output");
		double performance = VOSClusteringTechnique.calcPerformanceFunction();

        System.out.println("performance: " + performance);

    }

    private static Network readInputFile(String fileName) throws IOException {

        BufferedReader bufferedReader;
        double[] edgeWeight1, edgeWeight2, nodeWeight;
        int i, j, nEdges, nLines, nNodes, start, destination;
        int[] firstNeighborIndex, neighbor, nNeighbors, node1, node2;
        Network network;
        String[] splittedLine;

        // get the number of lines in the file
        bufferedReader = new BufferedReader(new FileReader(fileName));
        nLines = 0;
        while (bufferedReader.readLine() != null)
            nLines++;
        bufferedReader.close();

        // read in the start and destination nodes of each line and 
        // compute the total number of nodes
        bufferedReader = new BufferedReader(new FileReader(fileName));
        node1 = new int[nLines]; // first column 
        node2 = new int[nLines]; // second column
        edgeWeight1 = new double[nLines];
        i = -1;
        for (j = 0; j < nLines; j++) {

            // read in start and destination vertices
            splittedLine = bufferedReader.readLine().split(" ");
            start = Integer.parseInt(splittedLine[0]);
            destination = Integer.parseInt(splittedLine[1]);

            // enforce that the start node is lower than the destination node
            if (start < destination) {
                node1[j] = start;
                node2[j] = destination;
            } else if (start > destination) {
                node1[j] = destination;
                node2[j] = start;
            }

            // update the counter of number of nodes
            // i keeps track of the max node
            if (node1[j] > i)
                i = node1[j];
            if (node2[j] > i)
                i = node2[j];

            // store edge weights, if any
            edgeWeight1[j] = (splittedLine.length > 2) ? Double.parseDouble(
                splittedLine[2]) : 1;
        }
        nNodes = i + 1;
        bufferedReader.close();

        // compute the number of neighbors each node has
        nNeighbors = new int[nNodes];
        for (i = 0; i < nLines; i++) {
            nNeighbors[node1[i]]++;
            nNeighbors[node2[i]]++;
        }

        // finds number of half-edges
        firstNeighborIndex = new int[nNodes + 1];
        nEdges = 0;
        for (i = 0; i < nNodes; i++) {
            firstNeighborIndex[i] = nEdges;
            nEdges += nNeighbors[i];
        }
        firstNeighborIndex[nNodes] = nEdges;

        // computes the neighbor array and weight of each half-edge
        neighbor = new int[nEdges];
        edgeWeight2 = new double[nEdges];
        Arrays.fill(nNeighbors, 0);
        for (i = 0; i < nLines; i++) {
            j = firstNeighborIndex[node1[i]] + nNeighbors[node1[i]];
            neighbor[j] = node2[i];
            edgeWeight2[j] = edgeWeight1[i];
            nNeighbors[node1[i]]++;

            j = firstNeighborIndex[node2[i]] + nNeighbors[node2[i]];
            neighbor[j] = node1[i];
            edgeWeight2[j] = edgeWeight1[i];
            nNeighbors[node2[i]]++;
        }

        // construct the network based on the metric function
        network = new Network(
                nNodes, firstNeighborIndex, neighbor, edgeWeight2);
        
        return network;

    }

    private static Clustering readClusteringFile(
    	String fileName) throws IOException {

        BufferedReader bufferedReader;

        // get the number of lines in the file
        bufferedReader = new BufferedReader(new FileReader(fileName));
        int nLines = 0;
        while (bufferedReader.readLine() != null)
            nLines++;
        bufferedReader.close();

        // get the clustering
        bufferedReader = new BufferedReader(new FileReader(fileName));
        int[] cluster = new int[nLines];
        int node;
        int community;
        String[] splittedLine;
        for (int i = 0; i < nLines; i++) {
        	splittedLine = bufferedReader.readLine().split(" ");
            node = Integer.parseInt(splittedLine[0]);
            community = Integer.parseInt(splittedLine[1]);
            cluster[node] = community;
        }
        bufferedReader.close();

    	Clustering clustering = new Clustering(cluster);

    	return clustering;

    }

}
