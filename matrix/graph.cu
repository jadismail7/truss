
#include "graph.h"
#include <string>
__global__ void truss_kernel(Graph * g, int k, int *done)
{
    unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < g->nnzSize) {
        int v1 = g->row[index];
        int v2 = g->col[index];
        if (g->col[index] != UINT_MAX) {
            unsigned int commonNeighbors = 0;
            for (unsigned int i = g->rowPtr[v1]; i < g->rowPtr[v1 + 1]; ++i)
            {
                for (unsigned int j = g->rowPtr[v2]; j < j < g->rowPtr[v2 + 1]; ++j)
                {
                    unsigned int neighbor1 = g->col[i];
                    unsigned int neighbor2 = g->col[j];
                    if (neighbor1 == neighbor2 && g->col[i] != UINT_MAX && g->col[j] != UINT_MAX)
                    {
                        ++commonNeighbors;
                    }
                    if (commonNeighbors >= k-2) {
                        break;
                    }
                }
            }
            if (commonNeighbors <= k)
            {
                g->col[index] = UINT_MAX;
                done[0] = 0;
            }
        }
    }

}

void copyGraph(Graph * source, Graph * destination, cudaMemcpyKind direction) {
    // unsigned int * cooCol, *cooRow;
    // cudaMalloc((void **)&cooCol, h1->cooSize*sizeof(int));
    // cudaMalloc((void **)&cooRow, h1->cooSize*sizeof(int));
    // d1->cooSize = h1->cooSize;
    // unsigned int * cooColH = h1->col, *cooRowH=h1->row;

    // cudaMemcpy(cooCol, cooColH, sizeof(int)*h1->cooSize, direction);
    // cudaMemcpy(cooRow, cooRowH, sizeof(int)*h1->cooSize, direction);

    // d1->col = cooCol;
    // d1->row = cooRow;

    // printf("test2\n");
    // unsigned int *csrRowPtr, *csrCol;
    // cudaMalloc((void **) &csrRowPtr, h2->rowPtrSize*sizeof(int));
    // cudaMalloc((void **) &csrCol, h2->colSize*sizeof(int));
    // unsigned int * csrRowH = h2->rowPtr, *csrColH = h2->col;
    // d2->rowPtrSize = h2->rowPtrSize;
    // d2->colSize = h2->colSize;
    // cudaMemcpy(csrRowPtr, csrRowH, h2->rowPtrSize*sizeof(int), direction);
    // cudaMemcpy(csrCol, csrColH, h2->colSize*sizeof(int), direction);

    // d2->col = csrCol;
    // d2->rowPtr = csrRowPtr;
}
void truss_gpu(Graph * g, int k)
{
    // Graph * g
    // int *done_d;

    // cudaMalloc((void **)&coo_d, sizeof(COO));
    // cudaMalloc((void **)&csr_d, sizeof(CSR));
    // cudaMalloc((void **)&done_d, sizeof(int));

    // copyGraph(coo, coo_d, csr, csr_d, cudaMemcpyHostToDevice);
    // printf("Hello?");

    // cudaDeviceSynchronize();
    // unsigned int numThreads = 1024;
    // unsigned int numBlocks = (coo->cooSize + numThreads - 1) / numThreads;
    // int *done = (int *)malloc(sizeof(int));
    // done[0] = 0;

    // while (!done[0])
    // {

    //     done[0] = 1;
    //     cudaMemcpy(done_d, done, sizeof(int), cudaMemcpyHostToDevice);
    //     cudaDeviceSynchronize();
    //     truss_kernel<<<numThreads, numBlocks>>>(coo_d, csr_d, k, done_d);
    //     cudaDeviceSynchronize();
    //     cudaMemcpy(done, done_d, sizeof(int), cudaMemcpyDeviceToHost);
    //     cudaDeviceSynchronize();


    // }
    // printf("done");

    // copyGraph(coo_d, coo, csr_d, csr, cudaMemcpyDeviceToHost);

    // cudaFree(coo_d);
    // cudaFree(csr_d);
    // cudaFree(done_d);
}

Graph * mult_cpu(Graph * g, int k) {
    Graph * s = (Graph *) malloc(sizeof(Graph));
    s->nnzSize = 0;
    unsigned int index = 0;
    for (unsigned int i = 0 ; i < g->rowPtrSize; ++i) {
        for (unsigned int j = 0 ; j < g->rowPtrSize; ++j) {
            int sum = 0;
            for (unsigned int m1 = g->rowPtr[i]; m1 < g->rowPtr[i+1]; ++m1) {
                for (unsigned int m2 = g->rowPtr[j]; m2 < g->rowPtr[j+1]; ++m2) {
                    sum+= g->values[m1]*g->values[m2];
                }
            }
            addEdge(s, i, j, sum, index++);
        }
    }
    done(s, index);
    return s;
}

void truss_cpu(Graph * g, int k) {
    int changed = false;
    while(!changed) {
        changed = true;
        Graph * S = mult_cpu(g, k);
        for (int i = 0; i < S->nnzSize; ++i) {
            if (S->values[i] < k-2) {
                changed = false;
                g->values[i] = 0;
            }
            S->values[i] *= g->values[i];
        }

    }
    for (int i = 0 ; i < g->nnzSize; ++i) {
        if (g->values[i] >= k-2) {
            printf("%d - %d - %d\n", g->row[i], g->col[i], g->values[i]);
        }
    }

}


void addEdge(Graph * g, int v1, int v2, double weight, int index)
{

    if (g->nnzSize == 0)
    {
        g->row = (unsigned int *)malloc(BASE_SIZE * sizeof(int));
        g->col = (unsigned int *)malloc(BASE_SIZE * sizeof(int));
        g->values = (int *)malloc(BASE_SIZE * sizeof(int));
        g->nnzSize = BASE_SIZE;

    }
    if (index > g->nnzSize)
    {
        g->nnzSize *= 2;
        g->row = (unsigned int *)realloc(g->row, g->nnzSize * sizeof(int));
        g->col = (unsigned int *)realloc(g->col, g->nnzSize * sizeof(int));
        g->values = (int *)realloc(g->col, g->nnzSize * sizeof(int));
    }
    g->row[index] = v1;
    g->col[index] = v2;
    g->values[index] = weight;
}

void done(Graph * g, int size) {
    g->nnzSize = size;
}

void createCSRFromCOO(Graph * g, int numRows)
{
    // Allocate
    unsigned int *rowPtrs = (unsigned int *)calloc(numRows + 1, sizeof(unsigned int));
    unsigned int *colIdxs = (unsigned int *)malloc(g->nnzSize * sizeof(unsigned int));
    unsigned int *rowIdxs = (unsigned int *)malloc(g->nnzSize * sizeof(unsigned int));
    int *values = (int *)malloc(g->nnzSize * sizeof(int));

    // Histogram
    for (unsigned int i = 0; i <g->nnzSize; ++i)
    {
        unsigned int row = g->row[i];
        rowPtrs[row]++;
    }

    // Prefix sum
    unsigned int sum = 0;
    for (unsigned int row = 0; row < numRows; ++row)
    {
        unsigned int val = rowPtrs[row];
        rowPtrs[row] = sum;
        sum += val;
    }
    rowPtrs[numRows] = sum;

    // Binning
    for (unsigned int index = 0; index < g->nnzSize; ++index)
    {
        unsigned int row = g->row[index];
        unsigned int i = rowPtrs[row]++;
        colIdxs[i] = g->col[index];
        rowIdxs[i] = g->row[index];
        values[i] = g->values[index];
    }

    // Restore row pointers
    for (unsigned int row = numRows - 1; row > 0; --row)
    {
        rowPtrs[row] = rowPtrs[row - 1];
    }
    rowPtrs[0] = 0;

    g->rowPtrSize = numRows;
    g->rowPtr = rowPtrs;
    g->col = colIdxs;
    g->row = rowIdxs;
    g->values = values;
}

unsigned int *DFSUtil(Graph * g, int v, int visited[], int *size)
{

    unsigned int *visiting = (unsigned int *)malloc(BASE_SIZE * sizeof(int));

    visited[v] = 1;

    unsigned int i = 0;
    visiting[i++] = v;
    for (int j = g->rowPtr[v]; j < g->rowPtr[v + 1]; ++j)
    {
        if (g->col[j] != UINT_MAX && !visited[g->col[j]] )
        {
            int size2 = BASE_SIZE;
            unsigned int *temp = DFSUtil(g, g->col[j], visited, &size2);
            int index = i;
            for (i; i < index + size2; ++i)
            {
                if (i >= *size)
                {
                    *size *= 2;
                    visiting = (unsigned int *)realloc(visiting, (*size) * sizeof(int *));
                }
                visiting[i] = temp[i - index];
            }
        }
    }
    *size = i;
    return visiting;
}

unsigned int **connectedComponents(Graph * g, int *returnSize, int *componentSizes)
{
    unsigned int **dfs = (unsigned int **)malloc(BASE_SIZE * sizeof(int *));
    int index = 0;
    int size = BASE_SIZE;
    int *visited = (int *)malloc((g->rowPtrSize) * sizeof(int));
    for (int v = 0; v < g->rowPtrSize; v++){
        visited[v] = 0;
    }
    for (int v = 0; v < g->rowPtrSize; v++)
    {
        if (visited[v] == 0)
        {
            if (index >= size)
            {
                size *= 2;
                dfs = (unsigned int **)realloc(dfs, size * sizeof(int *));
                componentSizes = (int *)realloc(dfs, size * sizeof(int));
            }
            int size = BASE_SIZE;
            dfs[index] = (DFSUtil(g, v, visited, &size));
            componentSizes[index++] = size;
        }
    }
    *returnSize = index;
    return dfs;
}
void printTrussComponents(Graph * g, int k)
{
    int size;
    int *componentSizes = (int *)malloc(BASE_SIZE * sizeof(int));
    unsigned int **cc = connectedComponents(g, &size, componentSizes);

    for (int i = 0; i < size; ++i)
    {
        if (componentSizes[i] >= k)
        {
            printf("[");
            for (int j = 0; j < componentSizes[i]; ++j)
            {
                char end[3] = ", ";
                if (j == (componentSizes[i] - 1))
                {
                    strcpy(end, "]\n");
                }
                printf("%d%s", cc[i][j], end);
            }
        }
    }
}
