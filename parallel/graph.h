#include <stdio.h>
#include <stdlib.h>

#define BASE_SIZE 1500

typedef struct COO
{
    unsigned int cooSize;
    unsigned int *row;
    unsigned int *col;
} COO;

typedef struct CSR
{
    unsigned int *rowPtr;
    unsigned int *col;
    unsigned int rowPtrSize;
    unsigned int colSize;
} CSR;

void done(COO * coo, int size);

int initGraph(COO *coo, char *filename);

int intersect(double *s1, double *s2, int size);

void addEdge(COO *coo, int v1, int v2, double weight, int index);

CSR *createCSRfromCOO(COO *A, int numRows);

int *DFSUtil(COO *coo, CSR *csr, int v, int visited[], int *size);

int **connectedComponents(COO *coo, CSR *csr, int *returnSize, int *componentSizes);

void printTrussComponents(COO *coo, CSR *csr, int k);

void truss_gpu(COO *coo, CSR *csr, int k);
