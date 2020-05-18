#include <stdio.h>
#include <stdlib.h>

#define BASE_SIZE 1500

typedef struct graph
{
    unsigned int nnzSize;
    unsigned int *row;
    unsigned int *col;
    unsigned int *rowPtr;
    unsigned int rowPtrSize;
} Graph;


void done(Graph * g, int size);

int initGraph(Graph * g, char *filename);

int intersect(double *s1, double *s2, int size);

void addEdge(Graph * g, int v1, int v2, double weight, int index);

void createCSRFromCOO(Graph * g, int numRows);

unsigned int *DFSUtil(Graph * g, int v, int visited[], int *size);

unsigned int **connectedComponents(Graph * g, int *returnSize, int *componentSizes);

void printTrussComponents(Graph * g, int k);

void truss_gpu(Graph * g, int k);

void truss_cpu(Graph * g, int k);