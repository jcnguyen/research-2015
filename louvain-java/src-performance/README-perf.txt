
Java code for the Louvain method for performance. Original code can be found here:
http://www.ludowaltman.nl/slm/

For use in our complex network research: 
http://research.pomona.edu/complexnetworks/

Download
----------
The code can be found here: https://github.com/jcnguyen/research-2015
This documentation is specifically for running the Java code for the Louvain 
method using performance, which is found in “/research-2015/louvain-java/src-performance.”

The test suite can be found in the Box account under "Graphs". 

Running
----------
1. Download the code and test suite (see "Download").
2. Open the command line.
3. Change the directory to the src directory of the Louvain Method Java code. 
   $ cd ../research-2015/louvain-java/src-performance
   where ".." is the path to the research folder.
4. To compile,
   $ javac *.java
5. To run, 
   $ java ModularityOptimizer input output v resolution M n_random_starts n_iterations random_seed print

The command line arguments are as follows:
input            Input file name. Extension should be .pairs for unweighted 
                   graphs and .wpairs for weighted graphs.
output           Output file name. No extension - the algorithm will append the file endings.
v	       	  V is the scaling parameter that rates the importance of the
    		    weight of intercluster edges (with respect to the weight of 
		    the intra-cluster edges). Can be between 0 and 1, inclusive.
		    This value is unique to performance.
resolution       Value of the resolution parameter
M        	  Meaningful maximum of edge weights. This value is unique to 			    	    performance. Setting M is as less than the input graph’s
		    maximum edge weight will result in the algorithm setting all edge 			    weights greater than M equal to M.
n_random_starts  Number of random starts
n_iterations     Number of iterations per random start
random_seed      Seed of the random number generator
print            Whether or not to print output to the console 
                   0 = no
                   1 = yes

Test Run
----------
We run the test case using the Louvain method with modularity with the 
following arguments:
    input            ../test/test.pairs
    output           ../output-perf/test
    v	              1 = weight them equally
    resolution       1.0
    M       	      1 for unweighted graphs.
		        Varies for weighted graphs;
		        For suggested M’s for our test suite, 
			see the column M(perf) in https://docs.google.com/spreadsheets/d/11j7ESk_izWIaDsI2iMIFlhAaivLqxRDWDn5OOGQTkO0/edit?usp=sharing. 
			The “Max” column sets M to the max edge weight;
			the REHO column sets M to the max edge weight 
			excluding extreme high outliers. The percentage in parentheses
			is the % of data considered an extreme high outlier.
    n_random_starts  1
    n_iterations     1
    random_seed      0
    print            1 = yes

Input the following in the command line:
$ javac *.java
$ java ModularityOptimizer ../test/test.pairs ../output-perf/test 1 1.0 M 1 1 0 1

Expected output:
1. Output files in “/research-2015/louvain-java/output”:
	1. test.graph: final clustering output
	2. test.tree: the intermediate clusterings.
	3. test.info: to save a record of the command line output (see below), add &> test.info at the end.

2. Command line:
    -------------------------------------
    Modularity Optimizer
    Version 1.3.0
    by Ludo Waltman and Nees Jan van Eck
    -------------------------------------
    Reading input file...
    Finish reading input file.

    Number of nodes: 7
    Number of edges: 10
    Running Louvain algorithm - performance
    Scaling parameter V: 1
    Meaningful maximum M: 1

    Performance of unaltered graph: 0.3333333333333333

    Iteration: 1

	Pass 0
	network size: 4 nodes, 4 edges
	start computation: 01/18/2016 14:09:15
	end computation: 01/18/2016 14:09:15

	Pass 1
	network size: 2 nodes, 1 edges
	start computation: 01/18/2016 14:09:15
	end computation: 01/18/2016 14:09:15

	Pass 2
	network size: 1 nodes, 0 edges
	start computation: 01/18/2016 14:09:15
	end computation (1 node): 01/18/2016 14:09:15

    Iteration 1 performance: 0.6667

    Writing to ../output/test.tree...
    Finished writing to ../output/test.tree.

    Calculating the best overall clustering for all passes...
    Best clustering is Pass 0 with performance: 0.8333
    Writing to ../output/test.perfgraph...
    Finished writing to ../output/test.perfgraph.

    Maximum performance in 1 random starts: 0.6667
    Number of communities: 1
    Elapsed time: 0 seconds

    Writing to .graph...
    Finish writing to .graph.
