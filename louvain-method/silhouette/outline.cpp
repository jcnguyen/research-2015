allpairsshortestpath() = some structure called allPairsShortest
shortestDistance(vi. vk) // accesses allPairsShortest structure to get shortest path from vi to vk


numNodes = number of nodes in Community;
silhouetteWidths[numNodes]

// calculating the silhouette width for every vertex in the community
for (i=0;i<numNodes;i++){
	averageInnerDistance = 0;

	// finding the average of the shortest distances between vertex i and
	// every other vertex in the community
	for(k=0;k<numNodes;k++) {
		averageDistance += shortestDistance(vertex i, vertex k);
	}
	averageDistance=(averageDistance/numNodes);


	minimumOfAverageOutsideDistance = 0;
	distanceOptions[numCommunities-1];
	// finding the distances between vertex i, and every other vertex in
	// every other neighboring community
	for(n=0;n<numCommunties;n++){
		currentDistance = 0;
		numNodesInOtherCommunity = number of nodes in community n;
		for(k=0;k<numNodesInOtherCommunity; k++) {
			currentDistance +=shortestDistance(vertex i, vertex k);
		}
		currentDistance = currentDistance/numNodesInOtherCommunity;
		distanceOptions[n] = currentDistance;
	}

	minimumOfAverageOutsideDistance = min(distanceOptions);
	maxOfDistances = max(averageInnerDistance, minimumOfAverageOutsideDistance)
	silhouetteWidthOfi = (minimumOfAverageOutsideDistance-
							averageInnerDistance)/maxOfDistances)

	silhouetteWidths[i] = silhouetteWidthOfi
}

// calculating the silhouette width for the entire community
communitySilhouette = 0;
for (k=0;k<numNodes;k++) {
	communitySilhoette += silhoetteWidths[k]
}
communitySilhouette = communitySilhouette/numNodes

return communitySilhouette