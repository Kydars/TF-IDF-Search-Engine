// Centrality Measures API implementation
// COMP2521 Assignment 2
// 
// An API implementation which calculatess the closeness and betweenness 
// centralities of each vertices in a given graph.
//
// Author:
// Kyle Wu (z5363490@unsw.edu.au)
//
// Written: 19/03/2022

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "Graph.h"

// 	Function Prototypes
double countPaths(NodeData *nodeData, Vertex src, Vertex dest);
double countAppearances(NodeData *nodeData, Vertex i, Vertex src, Vertex dest);

//	Check CentralityMeasures.h for the description of this function.
double *closenessCentrality(Graph g) {
	int vertices = GraphNumVertices(g);
	double *new = calloc(vertices, sizeof(double));
	assert(new != NULL);
	for (int i = 0; i < vertices; i++) {
		NodeData *nodeData = dijkstra(g, i);
		double pathSum = 0;
		double nodeNum = 0;
		//	Calculate the shortest-path distance and count the number of 
		//	reachable nodes in the graph.
		for (int j = 0; j < vertices; j++) {
			if (nodeData[j].dist != INFINITY) {
				pathSum += nodeData[j].dist;
			}
			if (nodeData[j].pred != NULL) {
				nodeNum++;
			}
		}
		if (!pathSum) {
			new[i] = 0.0;
		} else {
			//	Use the Wasserman and Faust formula to calculate the closeness
			//	of vertex i.
			new[i] = ((nodeNum * nodeNum) / ((vertices - 1) * pathSum));
		}
		freeNodeData(nodeData, vertices);
	}
	return new;
}

//	Check CentralityMeasures.h for the description of this function.
double *betweennessCentrality(Graph g) {
	int vertices = GraphNumVertices(g);
	double *new = calloc(vertices, sizeof(double));
	assert(new != NULL);
	for (int i = 0; i < vertices; i++) {
		for (int src = 0; src < vertices; src++) {
			if (i == src) {
				continue;
			}
			NodeData *nodeData = dijkstra(g, src);
			//	Calculate the amount of appearances of i divided by the total 
			//	amount of paths, and store this value in the new, using the src
			//	value as the index.
			for (int dest = 0; dest < vertices; dest++) {
				if (i == dest || src == dest) {
					continue;
				}
				double totalPaths = countPaths(nodeData, src, dest);
				double appearances = countAppearances(nodeData, i, src, dest);
				if ((appearances) && (totalPaths)) {
					new[i] += appearances / totalPaths;
				}
			}
			freeNodeData(nodeData, vertices);
		}
	}
	return new;
}

//	Recursively calculate the number of shortest paths from the src vertex to 
//	the dest vertex.
double countPaths(NodeData *nodeData, Vertex src, Vertex dest) {
	if (dest == src) {
		return 1;
	}
	double count = 0;
	PredNode *curr = nodeData[dest].pred;
	while (curr != NULL) {
		count += countPaths(nodeData, src, curr->v);
		curr = curr->next;
	}
	return count;
}

//	Recursively count the number of appearances of i when traversing the 
//	shortest paths from the src vertex to the dest vertex.
double countAppearances(NodeData *nodeData, Vertex i, Vertex src, Vertex dest) {
	if (i == dest) {
		i = -1;
	}
	if (src == dest) {
		return (i == -1);
	}
	double count = 0;
	PredNode *curr = nodeData[dest].pred;
	while (curr != NULL) {
		count += countAppearances(nodeData, i, src, curr->v);
		curr = curr->next;
	}
	return count;
}

