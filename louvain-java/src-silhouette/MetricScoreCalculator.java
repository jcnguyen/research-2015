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

        String inputFileName
        int metricFunction;
        double resolution;
        Network network;
        VOSClusteringTechnique VOSClusteringTechnique;

        System.out.println("-------------------------------------");
        System.out.println("Metric Score Calculator");
        System.out.println("-------------------------------------");

        // read in the arguments
        inputFileName = args[0];
        metricFunction = Integer.parseInt(args[1]);
        resolution = Double.parseDouble(args[2]);

        // read input file
        network = readInputFile(inputFileName, metricFunction);

        // calculate metric scores
        VOSClusteringTechnique = new VOSClusteringTechnique(network, resolution);
// TODO
		VOSClusteringTechnique.calcModularityFunction();
        VOSClusteringTechnique.calcSilhouetteFunction();

    }

    /** 
     * Construct a network based on the metric function and an input file that 
     * contains the list of edges.
     *
     * @param  fileName        the input file
     * @param  metricFunction  the metric function
     * @throws IOException     occurs if there's an input or output error
     * @return a network based on the input file and the chosen metric function
     **/
    private static Network readInputFile(
        String fileName, int metricFunction) throws IOException {

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
        if (metricFunction != MODULARITY_ALTERNATIVE) {
            network = new Network(
                nNodes, firstNeighborIndex, neighbor, edgeWeight2);
        } else {
            nodeWeight = new double[nNodes];
            Arrays.fill(nodeWeight, 1);
            network = new Network(
                nNodes, nodeWeight, firstNeighborIndex, neighbor, edgeWeight2);
        }

        return network;

    }

}
