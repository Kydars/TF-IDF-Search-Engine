// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2 
// 
// An API implementation of Lance-Williams Algorithm.

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX
#define SINGLE_GAMMA 	-0.5
#define	COMPLETE_GAMMA 	 0.5

//	Function Prototypes
double **createDistArray(Graph g);
double **createArray(int vertices);
void updateShortestDistance(double **distArray, AdjList list, Vertex v);
DNode *makeDNode(Vertex v);
void shiftDendrogram(Dendrogram *dendA, DNode *cluster12, Vertex cluster1, 
	Vertex cluster2, int dendASize);
void shiftDistArray(double **distArray, double **shiftedDistArray, 
	Vertex cluster1, Vertex cluster2, int dendASize, int method);
void shiftSingleClusteredValues(double **distArray, double **shiftedDistArray, 
	Vertex i, Vertex cluster1, Vertex cluster2, int dendASize);
void shiftCompleteClusteredValues(double **distArray, double **shiftedDistArray, 
	Vertex i, Vertex cluster1, Vertex cluster2, int dendASize);
void shiftDistArrayValue(double **distArray, double **shiftedDistArray, 
	Vertex i, Vertex j, Vertex cluster1, Vertex cluster2);
void copyDistArray(double **distArray, double **shiftedDistArray,
	int dendASize);
void freeDistArray(double **shiftedDistArray, int size);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
	int vertices = GraphNumVertices(g);
	double **distArray = createDistArray(g);
	Dendrogram *dendA = malloc(vertices * sizeof(DNode));
	for (int i = 0; i < vertices; i++) {
		dendA[i] = makeDNode(i);
	}
	int dendASize = vertices;
	//	Apply Lance Willians algorithm continuously until the complete 
	//	Dendrogram is formed. 
	for (int iter = 0; iter < vertices; iter++) {
		//	Locate the minimum distance in the distArray and store its indices.
		double min = INFINITY;
		int cluster1 = 0, cluster2 = 0;
		for (int i = 0; i < dendASize; i++) {
			for (int j = 0; j < dendASize; j++) {
				if (min > distArray[i][j]) {
					min = distArray[i][j];
					cluster1 = i;
					cluster2 = j;
				}
			}
		}
		//	Swap the values of the clusters, such that cluster1 < cluster2.
		if (cluster1 > cluster2) {
			int temp = cluster1;
			cluster1 = cluster2;
			cluster2 = temp;
		}
		//	Merge the two clusters into a single cluster.
		DNode *cluster12 = makeDNode(iter);
		cluster12->left = dendA[cluster1];
		cluster12->right = dendA[cluster2];
		//	Shift the Dendrogram to the left by 1 index.
		shiftDendrogram(dendA, cluster12, cluster1, cluster2, dendASize);
		dendASize--;
		//	Create a new shifted distArray with dimensions one less than that
		//	of the main distArray. Then fill it will the shifted values from the
		//	original distArray, before copying the contents of the shifted
		//	distArray into the original one.
		double **shiftedDistArray = createArray(dendASize);
		shiftDistArray(distArray, shiftedDistArray, cluster1, cluster2, 
			dendASize, method);
			
		copyDistArray(distArray, shiftedDistArray, dendASize);
		freeDistArray(shiftedDistArray, dendASize);
	}
	Dendrogram finaldendA = dendA[0];
	freeDistArray(distArray, vertices);
	free(dendA);
	return finaldendA;
}

// 	Create a 2D double array which contains the shortest distances between two
//	vertices. If their is not direct edge between the vertices then the distance
//	will equal INFINITY.
double **createDistArray(Graph g) {
	int vertices = GraphNumVertices(g);
	double **distArray = createArray(vertices);
	for (Vertex v = 0; v < vertices; v++) {
		AdjList outIncident = GraphOutIncident(g, v);
		AdjList inIncident = GraphInIncident(g, v);
		updateShortestDistance(distArray, outIncident, v);
		updateShortestDistance(distArray, inIncident, v);
	}
	return distArray;
}

// 	Create a 2D double array and initialise the each value in the array to 
//	INFINITY.
double **createArray(int vertices) {
	double **array = malloc(vertices * sizeof(double *));
	for (int i = 0; i < vertices; i++) {
		array[i] = malloc(vertices * sizeof(double));
	}
	for (int i = 0; i < vertices; i++) {
		for (int j = 0; j < vertices; j++) {
			array[i][j] = INFINITY;
		}
	}
	return array;
}

//	Loop through an AdjList and if distance of the current edge is less than the
//	distance stored in the distArray for that edge then replace distArray 
//	distance for both directed edges between the two vertices.
void updateShortestDistance(double **distArray, AdjList list, Vertex v) {
	struct adjListNode *curr = list;
	while (curr != NULL) {
		double newDist = (double) 1 / curr->weight;
		if ((newDist < distArray[curr->v][v]) || 
			(newDist < distArray[v][curr->v])) 
		{
			distArray[curr->v][v] = newDist;
			distArray[v][curr->v] = newDist;
		}
		curr = curr->next;
	}
}

//	Create a dNode and set the Vertex to v and the left and right pointers to 
//	NULL.
DNode *makeDNode(Vertex v) {
	DNode *dNode = malloc(sizeof(DNode));
	dNode->vertex = v;
	dNode->left = NULL;
	dNode->right = NULL;
	return dNode;
}

//	Shift the Dendrogram to the left and add the new cluster12 into the second
//	right-most index of the Dendrogram array.
void shiftDendrogram(Dendrogram *dendA, DNode *cluster12, Vertex cluster1, 
	Vertex cluster2, int dendASize) 
{
	for (int i = 0; i < dendASize - 1; i++) {
		if (i == dendASize - 2) {
			dendA[i] = cluster12;
		} else {
			if (i >= cluster2 - 1) {
				dendA[i] = dendA[i + 2];
			} else if (i >= cluster1) {
				dendA[i] = dendA[i + 1];
			}
		}
	}
}

//	Update the shiftDistArray using the Lance Williams method, such that the
//	non-merged clusters are shifted to the left and upwards, overwriting the
//	values of the merged clusters. Then, fill the data of the new cluster into 
//	the right most and bottom most indices of the shiftedDistArray.
void shiftDistArray(double **distArray, double **shiftedDistArray, 
	Vertex cluster1, Vertex cluster2, int dendASize, int method) 
{
	for (Vertex i = 0; i < dendASize; i++) {
		for (Vertex j = i; j < dendASize; j++) {
			if (i == j) {
				continue;
			}
			// 	When j reaches the last row calculate the value of the data in
			//	the new cluster indices, based on the chosen method.
			// 	If method is 0, complete the Single Linkage method.
			//	If method is 1, complete the Complete Linkage method.
			if (j == dendASize - 1) {
				if (method == 1) {
					shiftSingleClusteredValues(distArray, shiftedDistArray, 
						i, cluster1, cluster2, dendASize);
				} else {
					shiftCompleteClusteredValues(distArray, shiftedDistArray, 
						i, cluster1, cluster2, dendASize);
				}
			//	When j is not the last row, shift the data in distArray into the
			//	correct position in shiftedDistIndex.
			} else {
				shiftDistArrayValue(distArray, shiftedDistArray, i, j,
				cluster1, cluster2);
			}
		}
	}
}

//	Single Linkage:
//	Calculate the minimum distances for each non-merged vertices, between that
//	vertex and two clusters and fill these distances inside of the 
//	shiftedDistArray.
void shiftSingleClusteredValues(double **distArray, double **shiftedDistArray, 
	Vertex i, Vertex cluster1, Vertex cluster2, int dendASize) 
{
	//	Check shiftDistArrayValue for logic of shiftedI.
	int shiftedI = i;
	if (shiftedI >= cluster2 - 1) {
		shiftedI += 2;
	} else if (i >= cluster1) {
		shiftedI += 1;
	}

	//	Find the distance between the shifted I index and the two clusters.
	double distIToCluster1 = distArray[cluster1][shiftedI];
	double distIToCluster2 = distArray[cluster2][shiftedI];

	//	Find the minium distance between i and the merged cluster.
	double min = INFINITY;
	if (distIToCluster1 < distIToCluster2) {
		min = distIToCluster1;
	} else {
		min = distIToCluster2;
	}
	shiftedDistArray[dendASize - 1][i] = min;
    shiftedDistArray[i][dendASize - 1] = min;
}

//	Complete Linkage:
//	Calculate the maximum distances for each non-merged vertices, between that
//	vertex and two clusters and fill these distances inside of the 
//	shiftedDistArray.
void shiftCompleteClusteredValues(double **distArray, double **shiftedDistArray, 
	Vertex i, Vertex cluster1, Vertex cluster2, int dendASize) 
{
	//	Check shiftDistArrayValue for logic of shiftedI.
	int shiftedI = i;
	if (shiftedI >= cluster2 - 1) {
		shiftedI += 2;
	} else if (i >= cluster1) {
		shiftedI += 1;
	}

	//	Find the distance between the shifted I index and the two clusters.
	double distIToCluster1 = distArray[cluster1][shiftedI];
	double distIToCluster2 = distArray[cluster2][shiftedI];

	//	Find the maximum non-infinite distance between i and the merged cluster.
	double max = 0;
	if (distIToCluster1 == INFINITY) {
		max = distIToCluster2;
	} else if (distIToCluster2 == INFINITY) {
		max = distIToCluster1;
	} else {
		if (distIToCluster1 > distIToCluster2) {
			max = distIToCluster1;
		} else {
			max = distIToCluster2;
		}
	}
	shiftedDistArray[dendASize - 1][i] = max;
    shiftedDistArray[i][dendASize - 1] = max;
}

//	Fill the shiftedDistArray with the corresponding values from distArray.
//	This function will shift the position of the vertices i or j based on 
//	position relative to the two merged clusters. 
//	If the vertex is after cluster2, then it will be shifted by two, if it is 
//	between cluster1 and cluster 2, then it will be shifted by one and if it is 
//	before cluster1 then it will not be shifted.
void shiftDistArrayValue(double **distArray, double **shiftedDistArray, 
	Vertex i, Vertex j, Vertex cluster1, Vertex cluster2) 
{
	if (i >= cluster2 - 1) {
		if (j >= cluster2 - 1) {
			shiftedDistArray[i][j] = distArray[i + 2][j + 2];
		} else if (j >= cluster1) {
			shiftedDistArray[i][j] = distArray[i + 2][j + 1];
		} else {
			shiftedDistArray[i][j] = distArray[i + 2][j];
		}
	} else if (i >= cluster1) {
		if (j >= cluster2 - 1) {
			shiftedDistArray[i][j] = distArray[i + 1][j + 2];
		} else if (j >= cluster1) {
			shiftedDistArray[i][j] = distArray[i + 1][j + 1];
		} else {
			shiftedDistArray[i][j] = distArray[i + 1][j];
		}
	} else {
		if (j > cluster2 - 1) {
			shiftedDistArray[i][j] = distArray[i][j + 2];
		} else if (j > cluster1) {
			shiftedDistArray[i][j] = distArray[i][j + 1];
		} else {
			shiftedDistArray[i][j] = distArray[i][j];
		}
	}
	shiftedDistArray[j][i] = shiftedDistArray[i][j];
}

//	Copy the contents of shiftedDistArray into distArray.
void copyDistArray(double **distArray, double **shiftedDistArray,
	int dendASize) 
{
	for (int i = 0; i < dendASize; i++) {
		for (int j = 0; j < dendASize; j++) {
			distArray[i][j] = shiftedDistArray[i][j];
		}
	}
}

// 	Free the contents of shiftedDistArray.
void freeDistArray(double **shiftedDistArray, int size) {
	for	(int i = 0; i < size; i++){
        free(shiftedDistArray[i]);
    }
    free(shiftedDistArray);
}
