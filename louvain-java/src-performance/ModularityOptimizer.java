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

    /** 
     * The main sequence. TODO
     *
     * @param  args         expects 9 arguments
     * @throws IOException  occurs if there's an input or output error
     **/
    public static void main(String[] args) throws IOException {
        boolean printOutput, update;
        Clustering clustering;
        Console console;
        double performance, maxPerformance, resolution, resolution2;
        int meaningfulMaxM, i, j, v_scalingParam, nIterations, nRandomStarts;
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
            v_scalingParam = Integer.parseInt(args[2]);
            resolution = Double.parseDouble(args[3]);
            meaningfulMaxM = Integer.parseInt(args[4]);
            nRandomStarts = Integer.parseInt(args[5]);
            nIterations = Integer.parseInt(args[6]);
            randomSeed = Long.parseLong(args[7]);
            printOutput = (Integer.parseInt(args[8]) > 0);
        } else {
            console = System.console();
            inputFileName = console.readLine("Input file path/name: ");
            outputFileName = console.readLine("Output file path/name (without file extension): ");
            v_scalingParam = Integer.parseInt(console.readLine("Scaling parameter V (0 to 1; rates the" + 
                                                            " importance of the weight of inter-cluster edges," + 
                                                            " with respect to the weight of intra-cluster edges): ")); 
            resolution = Double.parseDouble(console.readLine("Resolution parameter (e.g., 1.0): "));
            meaningfulMaxM = Integer.parseInt(console.readLine("Meaningful maximum M of edge weights: "));
            nRandomStarts = Integer.parseInt(console.readLine("Number of random starts (e.g., 10): "));
            nIterations = Integer.parseInt(console.readLine("Number of iterations (e.g., 10): "));
            randomSeed = Long.parseLong(console.readLine("Random seed (e.g., 0): "));
            printOutput = (Integer.parseInt(console.readLine("Print output (0 = no; 1 = yes): ")) > 0);
            System.out.println();
        }

        // read input file
        if (printOutput) System.out.println("Reading input file...");
        network = readInputFile(inputFileName);
        if (printOutput) System.out.println("Finish reading input file.");

        // print network characteristics
        if (printOutput) {
            System.out.println();
            System.out.format("Number of nodes: %d%n", network.getNNodes());
            System.out.format("Number of edges: %d%n", network.getNEdges());
            System.out.println("Running Louvain algorithm - performance");
            System.out.println("Scaling parameter V: " + v_scalingParam);
            System.out.println("Meaningful maximum M: " + meaningfulMaxM);
            System.out.println();
        }


        /** 
        * Run the Louvain algorithm nRandomStarts times,
        * beginning with a random node configuration every time.
        * After all the trials, return the trial that resulted in the
        * best performance score. 
        */
        beginTime = System.currentTimeMillis();
        clustering = null;
        maxPerformance = Double.NEGATIVE_INFINITY;
        random = new Random(randomSeed);

        for (i = 0; i < nRandomStarts; i++) {
            if (printOutput && (nRandomStarts > 1)) {
                System.out.format("\tRandom start: %d%n", i + 1);
            }

            VOSClusteringTechnique = new VOSClusteringTechnique(v_scalingParam, meaningfulMaxM, 
                                                                network, resolution, outputFileName);

            j = 0;
            update = true;

            // Print performance of unaltered graph once
            if (i == 0) {
                System.out.println("Performance of unaltered graph: " + VOSClusteringTechnique.calcPerformanceFunction());
            }

            // run for nIterations times, or until there is no increase in performance
            do {
                if (printOutput && (nIterations > 0)) {
                    System.out.format("\nIteration: %d%n", j + 1);
                    update = VOSClusteringTechnique.runLouvainAlgorithm(random);
                }

                j++;

                performance = VOSClusteringTechnique.calcPerformanceFunction();

                if (printOutput && (nIterations > 0))
                    System.out.format("\nIteration %d performance: %.4f%n", j, performance);

            } while ((j < nIterations) && update);

            // we have found our best clustering
            if (performance > maxPerformance) {
                clustering = VOSClusteringTechnique.getClustering();
                maxPerformance = performance;
            }

            /* Print all clusterings to .tree file */
            VOSClusteringTechnique.printClusteringsToTree();

            /* Find the best clustering and print it to .perfgraph */
            VOSClusteringTechnique.bestOverallClustering();
        }

        endTime = System.currentTimeMillis();

        // print end-of-algorithm information
        if (printOutput) {
            System.out.format("\nMaximum performance in %d random starts: %.4f%n", 
                    nRandomStarts, maxPerformance);
            System.out.format("Number of communities: %d%n", 
                clustering.getNClusters());
            System.out.format("Elapsed time: %d seconds%n", 
                Math.round((endTime - beginTime) / 1000.0));
            System.out.println();
        }

        // write to output file
        if (printOutput) System.out.println("Writing to .graph...");
        writeOutputFile(outputFileName, clustering);
        if (printOutput) System.out.println("Finish writing to .graph.");
    }

    /** 
     * Constructs a network based on the modularity function and 
     * an input file that contains the list of edges.
     *
     * @param  fileName            the input file
     * @throws IOException         occurs if there's an input or output error
     * @return a network based on the input file and the chosen modularity 
     *         function
     **/
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

        // read in the start and end nodes of each line and 
        // compute the total number of nodes
        bufferedReader = new BufferedReader(new FileReader(fileName));
        node1 = new int[nLines];    // first column 
        node2 = new int[nLines];    // second column
        edgeWeight1 = new double[nLines];
        i = -1;
        for (j = 0; j < nLines; j++) {

            // read in start and destination vertices
            splittedLine = bufferedReader.readLine().split(" ");
            start = Integer.parseInt(splittedLine[0]);
            destination = Integer.parseInt(splittedLine[1]);

            // enforce the starting node is lower than the dest node
            if (start < destination) {
                node1[j] = start;
                node2[j] = destination;
            } else if (start > destination) {
                node1[j] = destination;
                node2[j] = start;
            }

            // update the counter of number of nodes
            // i keeps track of the max 
            if (node1[j] > i)
                i = node1[j];
            if (node2[j] > i)
                i = node2[j];

            // if there's a 3rd argument in the line, we know we have edge weights
            // store weights in edgeWeight1
            edgeWeight1[j] = (splittedLine.length > 2) ? Double.parseDouble(
                splittedLine[2]) : 1;
        }
        nNodes = i + 1;
        bufferedReader.close();

        /** DONE READING INPUT - 
        now calculate useful information about the graph **/

        // compute the number of neighbors each node has
        nNeighbors = new int[nNodes];
        for (i = 0; i < nLines; i++) {
            nNeighbors[node1[i]]++;
            nNeighbors[node2[i]]++;
        }

        // count the number of edges
        // set firstNeighborIndex[i]-firstNeighborValue[i-1] = nNeighbors[i-1]

        // finds number of half-edges
        // TODO what's the purpose of firstNeighborIndex?
        firstNeighborIndex = new int[nNodes + 1];
        nEdges = 0;
        for (i = 0; i < nNodes; i++) {
            firstNeighborIndex[i] = nEdges;
            nEdges += nNeighbors[i];
        }
        firstNeighborIndex[nNodes] = nEdges;

        // TODO delete: prints out firstNeighborIndex
        // System.out.println("index firstNeighborIndex");
        // for (int jj = 0; jj < nNodes + 1; jj++) {
        //     System.out.println(jj + " " + firstNeighborIndex[jj]);
        // }

        // computes the neighbor array and weight of each half-edge
        // i'th spot in neighbor is the destination vertex d of the edge between v_c and d,
        // where v_c is the vertex corresponding to the "chunk" of neighbor that contains i
        // i'th spot in edgeWeight2 is the weight of the edge indicated in 
        neighbor = new int[nEdges];
        edgeWeight2 = new double[nEdges];

        // initialize nNeighbors to 0
        Arrays.fill(nNeighbors, 0);

        // for every edge, in the order given in the input file
        // initialize 
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

        // construct the network 
        
            network = new Network(
                nNodes, firstNeighborIndex, neighbor, edgeWeight2);

        return network;
    }

    /** 
     * Write the clusterings out into a file.
     *
     * @param  fileName     the file to write out to
     * @param  clustering   the clusters to write out
     * @throws IOException  occurs if there's an input or output error
     **/
    private static void writeOutputFile(String fileName, Clustering clustering
        ) throws IOException {

        BufferedWriter bufferedWriter;
        int i, nNodes;

        // get the information
        nNodes = clustering.getNNodes();
        clustering.orderClustersByNNodes();

        // writing to the .graph: final communities output
        fileName = fileName + ".graph";
        bufferedWriter = new BufferedWriter(new FileWriter(fileName));

        for (i = 0; i < nNodes; i++) {
            bufferedWriter.write(i + " " + Integer.toString(clustering.getCluster(i))); 
            bufferedWriter.newLine();
        }

        bufferedWriter.close();
    }


}
