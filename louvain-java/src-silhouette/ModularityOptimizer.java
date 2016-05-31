/**
 * ModularityOptimizer
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.0, 08/31/15
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Console;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

public class ModularityOptimizer {

    private static final int MODULARITY_STANDARD = 1;
    private static final int MODULARITY_ALTERNATIVE = 2;
    private static final int SILHOUETTE_INDEX = 3;

    /** 
     * The main sequence. 
     * Runs the complete algorithm to maximize the metric score.
     *
     * @param  args         expects 9 arguments
     * @throws IOException  occurs if there's an input or output error
     **/
    public static void main(String[] args) throws IOException {
        boolean printOutput, update;
        Clustering clustering;
        Console console;
        double metricValue, maxMetricValue, resolution, resolution2;
        int algorithm, i, j, metricFunction, nIterations, nRandomStarts;
        long beginTime, endTime, randomSeed;
        Network network;
        Random random;
        String inputFileName, outputFileName;
        VOSClusteringTechnique VOSClusteringTechnique;

        System.out.println("-------------------------------------");
        System.out.println("Modularity Optimizer");
        System.out.println("Version 1.3.0");
        System.out.println("by Ludo Waltman and Nees Jan van Eck");
        System.out.println("-------------------------------------");

        // read in the arguments
        if (args.length == 9) {
            inputFileName = args[0];
            outputFileName = args[1];
            metricFunction = Integer.parseInt(args[2]);
            resolution = Double.parseDouble(args[3]);
            algorithm = Integer.parseInt(args[4]);
            nRandomStarts = Integer.parseInt(args[5]);
            nIterations = Integer.parseInt(args[6]);
            randomSeed = Long.parseLong(args[7]);
            printOutput = (Integer.parseInt(args[8]) > 0);
        } else {
            console = System.console();
            inputFileName = console.readLine("Input file name: ");
            outputFileName = console.readLine("Output file name: ");
            metricFunction = Integer.parseInt(console.readLine("Metric function (1 = modularity standard; 2 = modularity alternative; 3 = silhouette index; 4 = performance): ")); 
            resolution = Double.parseDouble(console.readLine("Resolution parameter (e.g. 1.0): "));
            algorithm = Integer.parseInt(console.readLine("Algorithm (1 = Louvain; 2 = Louvain with multilevel refinement; 3 = smart local moving): "));
            nRandomStarts = Integer.parseInt(console.readLine("Number of random starts (e.g., 10): "));
            nIterations = Integer.parseInt(console.readLine("Number of iterations (e.g., 10): "));
            randomSeed = Long.parseLong(console.readLine("Random seed (e.g., 0): "));
            printOutput = (Integer.parseInt(console.readLine("Print output (0 = no; 1 = yes): ")) > 0);
            System.out.println();
        }

        // read input file
        if (printOutput) System.out.println("Reading input file...");
        network = readInputFile(inputFileName, metricFunction);
        if (printOutput) System.out.println("Finish reading input file.");

        // print network characteristics
        if (printOutput) {
            System.out.println();
            System.out.format("Number of nodes: %d%n", network.getNNodes());
            System.out.format("Number of edges: %d%n", network.getNEdges());
            System.out.println();
            System.out.print("Running " + 
                ((algorithm == 1) ? "Louvain algorithm" : 
                ((algorithm == 2) ? "Louvain algorithm with multilevel refinement" : 
                "smart local moving algorithm")));
            System.out.println(" with metric " + 
                ((metricFunction == MODULARITY_ALTERNATIVE) ? "alternative modularity" : 
                ((metricFunction == SILHOUETTE_INDEX) ? "silhouette index" : 
                "standard modularity")));
            System.out.println();
        }

        // calculate resolution based on metric function
        resolution2 = ((metricFunction != MODULARITY_ALTERNATIVE) ? 
            (resolution / (2 * network.getTotalEdgeWeight() + 
                network.totalEdgeWeightSelfLinks)) : 
            resolution);

        beginTime = System.currentTimeMillis();
        clustering = null;
        maxMetricValue = Double.NEGATIVE_INFINITY;
        random = new Random(randomSeed);

        for (i = 0; i < nRandomStarts; i++) {
            if (printOutput && (nRandomStarts > 1))
                System.out.format("\tRandom start: %d%n", i + 1);

            VOSClusteringTechnique = new VOSClusteringTechnique(network, resolution2);

            j = 0;
            update = true;
            do {
                if (printOutput && (nIterations > 1))
                    System.out.format("\tIteration %d", (j + 1));

                // run the algorithm
                if (algorithm == 1)
                    update = VOSClusteringTechnique.runLouvainAlgorithm(random, metricFunction);
                else if (algorithm == 2)
                    update = VOSClusteringTechnique.runLouvainAlgorithmWithMultilevelRefinement(random);
                else if (algorithm == 3)
                    VOSClusteringTechnique.runSmartLocalMovingAlgorithm(random);
                j++;

                // calculate the metric score of the final clustering
                if (metricFunction == MODULARITY_STANDARD || metricFunction == MODULARITY_ALTERNATIVE) {
                    metricValue = VOSClusteringTechnique.calcModularityFunction();
                } else if (metricFunction == SILHOUETTE_INDEX) {
                    metricValue = VOSClusteringTechnique.calcSilhouetteFunction();
                } else { 
                    metricValue = VOSClusteringTechnique.calcModularityFunction();
                }

                if (printOutput && (nIterations > 1))
                    System.out.format("\tMetric value: %.4f%n", metricValue);
            }
            while ((j < nIterations) && update);

            if (metricValue > maxMetricValue) {
                clustering = VOSClusteringTechnique.getClustering();
                maxMetricValue = metricValue;
            }

            if (printOutput && (nRandomStarts > 1)) {
                if (nIterations == 1)
                    System.out.format("\tMetric value: %.4f%n", metricValue);
                System.out.println();
            }
        }
        endTime = System.currentTimeMillis();

        // print end-of-algorithm information
        if (printOutput) {
            if (nRandomStarts == 1) {
                if (nIterations > 1)
                    System.out.println();
                System.out.format("Metric value: %.4f%n", maxMetricValue);
            } else {
                System.out.format("Maximum metric value in %d random starts: %.4f%n", nRandomStarts, maxMetricValue);
            }

            System.out.format("Number of communities: %d%n", clustering.getNClusters());
            System.out.format("Elapsed time: %d seconds%n", Math.round((endTime - beginTime) / 1000.0));
            System.out.println();
        }

        // write to output file
        if (printOutput) System.out.println("Writing output file...");
        writeOutputFile(outputFileName, clustering);
        if (printOutput) System.out.println("Finish writing to output file.");

        // because writeOutputFile reorders the clustering, the SI score t be recalculated to reflect the reordered clusters
        if (metricFunction == SILHOUETTE_INDEX) {
            Network nt = readInputFile(inputFileName, metricFunction);
            Clustering cl = readClusteringFile(outputFileName);
            double res = 1/(2 * nt.getTotalEdgeWeight() + nt.totalEdgeWeightSelfLinks);
            VOSClusteringTechnique vos = new VOSClusteringTechnique(nt, cl, res);

            // calculate metric scores
            double silhouette = vos.calcSilhouetteFunction();
            System.out.println();
            System.out.format("SI value after reordering clusters: %.4f%n", silhouette);
        }
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

    /** 
     * Write the clusterings out into a file.
     *
     * @param  fileName     the file to write out to
     * @param  clustering   the clusters to write out
     * @throws IOException  occurs if there's an input or output error
     **/
    private static void writeOutputFile(
        String fileName, Clustering clustering) throws IOException {

        BufferedWriter bufferedWriter;
        int i, nNodes;

        nNodes = clustering.getNNodes();
        clustering.orderClustersByNNodes();
        bufferedWriter = new BufferedWriter(new FileWriter(fileName));
        for (i = 0; i < nNodes; i++) {
            bufferedWriter.write(
                i + " " + Integer.toString(clustering.getCluster(i)));
            bufferedWriter.newLine();
        }
        bufferedWriter.close();
    }

    /** 
     * Construct a clustering based on a file.
     *
     * @param  fileName        the clustering file
     * @throws IOException     occurs if there's an input or output error
     * @return a clustering based on the cluster file
     **/
    private static Clustering readClusteringFile(String fileName) throws IOException {

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
