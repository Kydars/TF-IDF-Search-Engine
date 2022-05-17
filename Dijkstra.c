// Dijkstra API implementation
// COMP2521 Assignment 2
// 
// An API implementation of Dijkstra's Algorithm.

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

PredNode *predReplace(PredNode *head, Vertex v);
PredNode *predInsert(PredNode *head, Vertex v);
void freePredNode(PredNode *pred);

// 	Check Dijkstra.h for the description of this function.
NodeData *dijkstra(Graph g, Vertex src) {
	int vertices = GraphNumVertices(g);
	NodeData *new = malloc((vertices) * sizeof(NodeData));
	assert(new != NULL);
	PQ pq = PQNew(); 
	for (int i = 0; i < vertices; i++) {
		new[i].dist = INFINITY;
		new[i].pred = NULL;
	}
	new[src].dist = 0;
	PQInsert(pq, src, new[src].dist);
	// 	While the priority queue is not empty, dequeue the lowest weighted 
	//	element, generate an out incident graph for it and loop through the
	//	incident graph.
	while (!PQIsEmpty(pq)) {
		int deqItem = PQDequeue(pq);
		AdjList curr = GraphOutIncident(g, deqItem);
		while (curr != NULL) {
			int newDist = new[deqItem].dist + curr->weight;
			//	If the distance to the curr vertex through deqItem is less than 
			//	the curr distance in new, replace the curr distance, insert the 
			//	curr vertex into the queue and replace the curr PredNode with a 
			//	new one containing deqItem as the vertex.
			if (newDist < new[curr->v].dist) {
				new[curr->v].dist = newDist;
				PQInsert(pq, curr->v, new[curr->v].dist);
				new[curr->v].pred = predReplace(new[curr->v].pred, deqItem);
				// 	If the distance to the curr vertex through deqItem is equal 
				//	to the curr distance in new, insert a new PredNode 
				//	containing deqItem as the vertex.
			} else if (newDist == new[curr->v].dist) {
				new[curr->v].pred = predInsert(new[curr->v].pred, deqItem);
			}
			curr = curr->next;
		}
	}
	PQFree(pq);
	return new;
}

//	Replace the current PredNode with a new PredNode containing the vertex v.
PredNode *predReplace(PredNode *head, Vertex v) {
	freePredNode(head);
	PredNode *new = malloc(sizeof(PredNode));
	assert(new != NULL);
	new->v = v;
	new->next = NULL;
	return new;
}

//	Insert a new PredNode containing the vertex v into the PredNode linked 
//	list, in ascending order.
PredNode *predInsert(PredNode *head, Vertex v) {
	PredNode *new = malloc(sizeof(PredNode));
	assert(new != NULL);
	new->v = v;
	//	If v is less than the head vertex insert new at the head of the list.
	if (v < head->v) {
		new->next = head;
		return new;
	}
	//	Loop through the linked list and insert new in the correct ascending
	//	position.
	PredNode *curr = head;
	while (curr->next != NULL) {
		if (v < curr->next->v) {
			PredNode *temp = curr->next;
			curr->next = new;
			new->next = temp;
			return head;
		}
		curr = curr->next;
	}
	new->next = NULL;
	curr->next = new;
	return head;
}

// 	Check Dijkstra.h for the description of this function.
void freeNodeData(NodeData *data, int nV) {
	for (int i = 0; i < nV; i++) {
		freePredNode(data[i].pred);
	}
	free(data);
}

//	Loop through the pred linked list and free each node iteratively.
void freePredNode(PredNode *pred) {
	while (pred != NULL) {
		PredNode *tmp = pred;
		pred = pred->next;
		free(tmp);
	}
}
