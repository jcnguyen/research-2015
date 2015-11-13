// C Program for Floyd Warshall Algorithm
#include <stdlib.h>
 
/* Define Infinite as a large enough value. This value will be used
  for vertices not connected to each other */
#define INF 9999
 
// // Solves the all-pairs shortest path problem using Floyd Warshall algorithm
// void floydWarshall1 (int graph[][nb_nodes]) {

// 	printf ("PRINTING INPUT MATRIX FOR FW \n");
// 	for (int i = 0; i < nb_nodes; i++)
// 	{
// 		for (int j = 0; j < nb_nodes; j++)
// 		{
// 			if (graph[i][j] == INF)
// 				printf("%7s", "INF");
// 			else
// 				printf ("%7d", graph[i][j]);
// 		}
// 		printf("\n");
// 	}
  
// 	/* dist[][] will be the output matrix that will finally have the shortest 
// 	  distances between every pair of vertices */
// 	int dist[nb_nodes][nb_nodes], i, j, k;
 
// 	/* Initialize the solution matrix same as input graph matrix. Or 
// 	   we can say the initial values of shortest distances are based
// 	   on shortest paths considering no intermediate vertex. */
// 	printf ("INITIALIZING DISTANCE MATRIX\n");
// 	for (i = 0; i < nb_nodes; i++) {
// 		for (j = 0; j < nb_nodes; j++) {
// 			if (graph[i][j] == 0 & i != j) {
// 				dist[i][j] == INF;
// 				printf("%7s", "INF");
// 			} else {
// 				dist[i][j] = graph[i][j];
// 				printf ("%7d", dist[i][j]);
				
// 			}
// 		}
// 		printf ("\n");
// 	}
 
//  	printf ("----------------------------------\n");
// 	 Add all vertices one by one to the set of intermediate vertices.
// 	  ---> Before start of a iteration, we have shortest distances between all
// 	  pairs of vertices such that the shortest distances consider only the
// 	  vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
// 	  ----> After the end of a iteration, vertex no. k is added to the set of
// 	  intermediate vertices and the set becomes {0, 1, 2, .. k} 
// 	for (k = 0; k < nb_nodes; k++) {
// 		// Pick all vertices as source one by one
// 		for (i = 0; i < nb_nodes; i++) {
// 			// Pick all vertices as destination for the
// 			// above picked source
// 			for (j = 0; j < nb_nodes; j++) {
// 				// If vertex k is on the shortest path from
// 				// i to j, then update the value of dist[i][j]
// 				// printf("(%7d,%7d,%7d)", k, i, j);
// 				if (dist[i][k] + dist[k][j] < dist[i][j]) {
// 					dist[i][j] = dist[i][k] + dist[k][j];
// 					printf ("%7d\n", dist[i][j]);
// 				}
// 			}
// 		}
// 	}
 
// 	// Print the shortest distance matrix
// 	// printSolution(dist);

// 	printf ("PRINTING DISTANCE MATRIX FOR FW \n");
// 	for (int i = 0; i < nb_nodes; i++)
// 	{
// 		for (int j = 0; j < nb_nodes; j++)
// 		{
// 			if (dist[i][j] == INF)
// 				printf("%7s", "INF");
// 			else
// 				printf ("%7d", dist[i][j]);
// 		}
// 		printf("\n");
// 	}
// }

int** floydWarshall (unsigned int nb_nodes, int** graph) {

	// printf ("PRINTING INPUT MATRIX FOR FW \n");
	// for (int i = 0; i < nb_nodes; i++)
	// {
	// 	for (int j = 0; j < nb_nodes; j++)
	// 	{
	// 		if (graph[i][j] == INF)
	// 			printf("%7s", "INF");
	// 		else
	// 			printf ("%7d", graph[i][j]);
	// 	}
	// 	printf("\n");
	// }
  
	/* dist[][] will be the output matrix that will finally have the shortest 
	  distances between every pair of vertices */
	int i, j, k;

	// initializing adjMatrix
	int** dist = new int*[nb_nodes];
	for (int i = 0; i < nb_nodes; ++i) {
		dist[i] = new int[nb_nodes];
	}
 
	/* Initialize the solution matrix same as input graph matrix. Or 
	   we can say the initial values of shortest distances are based
	   on shortest paths considering no intermediate vertex. */
	// printf ("INITIALIZING DISTANCE MATRIX\n");
	for (i = 0; i < nb_nodes; i++) {
		for (j = 0; j < nb_nodes; j++) {
			if (graph[i][j] == 0 & i != j) {
				dist[i][j] == INF;
				// printf("%7s", "INF");
			} else {
				dist[i][j] = graph[i][j];
				// printf ("%7d", dist[i][j]);
				
			}
		}
		printf ("\n");
	}
 
 	// printf ("----------------------------------\n");
	/* Add all vertices one by one to the set of intermediate vertices.
	  ---> Before start of a iteration, we have shortest distances between all
	  pairs of vertices such that the shortest distances consider only the
	  vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
	  ----> After the end of a iteration, vertex no. k is added to the set of
	  intermediate vertices and the set becomes {0, 1, 2, .. k} */
	for (k = 0; k < nb_nodes; k++) {
		// Pick all vertices as source one by one
		for (i = 0; i < nb_nodes; i++) {
			// Pick all vertices as destination for the
			// above picked source
			for (j = 0; j < nb_nodes; j++) {
				// If vertex k is on the shortest path from
				// i to j, then update the value of dist[i][j]
				// printf("(%7d,%7d,%7d)", k, i, j);
				if (dist[i][k] + dist[k][j] < dist[i][j]) {
					dist[i][j] = dist[i][k] + dist[k][j];
					// printf ("%7d\n", dist[i][j]);
				}
			}
		}
	}
 
	// Print the shortest distance matrix
	printf ("PRINTING DISTANCE MATRIX FOR FW \n");
	for (int i = 0; i < nb_nodes; i++)
	{
		for (int j = 0; j < nb_nodes; j++)
		{
			if (dist[i][j] == INF)
				printf("%7s", "INF");
			else
				printf ("%7d", dist[i][j]);
		}
		printf("\n");
	}
 	return dist;
}