
#include "graph.h"
#include "../timer.h"
#include <string>
__global__ void truss_kernel(Graph g, int k, int *done)
{
    unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < g.nnzSize) {
        unsigned int v1 = g.row[index];
        unsigned int v2 = g.col[index];
        if (v2 != UINT_MAX) {
            unsigned int commonNeighbors = 0;
            int index1 = g.rowPtr[v1];
            int index2 = g.rowPtr[v2];

            while (index1 < g.rowPtr[v1+1] && index2 < g.rowPtr[v2+1]) {
                if (g.col[index1] == UINT_MAX) {
                    ++index1;
                    continue;
                } else if (g.col[index2] == UINT_MAX) {
                    ++index2;
                    continue;
                }
                if (g.col[index1] < g.col[index2]) {
                    ++index1;
                } else if (g.col[index1] > g.col[index2]) {
                    ++index2; 
                } else {
                    ++commonNeighbors; ++index1; ++index2;
                }
                if (commonNeighbors == k-2) {
                    break;
                }
            }
            if (commonNeighbors < k-2)
            {
                g.col[index] = UINT_MAX;
                *done = 0;
            }
        }
    }

}

void truss_gpu(Graph * g, int k)
{
    Graph graph_d;
    int *done_d;

    Timer timer;
    startTime(&timer);
    cudaMalloc((void**)&graph_d.row, g->nnzSize*sizeof(int));
    cudaMalloc((void**)&graph_d.col, g->nnzSize*sizeof(int));
    cudaMalloc((void**)&graph_d.rowPtr, g->rowPtrSize*sizeof(int));
    cudaMalloc((void **)&done_d, sizeof(int));
    cudaDeviceSynchronize();
    stopTime(&timer);
    printElapsedTime(timer, "Allocation");
    startTime(&timer);
    cudaMemcpy(graph_d.row, g->row, g->nnzSize*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(graph_d.col, g->col, g->nnzSize*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(graph_d.rowPtr, g->rowPtr, g->rowPtrSize*sizeof(int), cudaMemcpyHostToDevice);
    
    graph_d.nnzSize = g->nnzSize;
    graph_d.rowPtrSize = g->rowPtrSize;
    cudaDeviceSynchronize();
    stopTime(&timer);
    printElapsedTime(timer, "Copy");

    unsigned int numThreads = 1024;
    unsigned int numBlocks = (g->nnzSize + numThreads - 1) / numThreads;
    int *done = (int *)malloc(sizeof(int));
    *done = 0;

    startTime(&timer);
    while (!*done)
    {
        *done = 1;
        cudaMemcpy(done_d, done, sizeof(int), cudaMemcpyHostToDevice);
        truss_kernel<<<numThreads, numBlocks>>>(graph_d, k, done_d);
        cudaMemcpy(done, done_d, sizeof(int), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
    }
    stopTime(&timer);
    printElapsedTime(timer, "Kernel");
    
    startTime(&timer);
    cudaMemcpy(g->row, graph_d.row, g->nnzSize*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(g->col, graph_d.col, g->nnzSize*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(g->rowPtr, graph_d.rowPtr, g->rowPtrSize*sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    stopTime(&timer);
    printElapsedTime(timer, "Copy back");

    startTime(&timer);
    cudaFree(graph_d.row);
    cudaFree(graph_d.col);
    cudaFree(graph_d.rowPtr);
    cudaFree(&graph_d);
    cudaFree(done_d);
    cudaDeviceSynchronize();
    stopTime(&timer);
    printElapsedTime(timer, "Free");
}

void truss_cpu(Graph * g, int k) {
    int done = 0;
    while(!done){
        done = 1;
        for (unsigned int i = 0; i < g->nnzSize; ++i) {
            unsigned int v1 = g->row[i];
            unsigned int v2 = g->col[i];
            if (g->col[i] == UINT_MAX){
                continue;
            }
            unsigned int commonNeighbors = 0;
            int index1 = g->rowPtr[v1];
            int index2 = g->rowPtr[v2];

            while (index1 < g->rowPtr[v1+1] && index2 < g->rowPtr[v2+1]) {
                if (g->col[index1] == UINT_MAX) {
                    ++index1;
                    continue;
                } else if (g->col[index2] == UINT_MAX) {
                    ++index2;
                    continue;
                }
                if (g->col[index1] < g->col[index2]) {
                    ++index1;
                } else if (g->col[index1] > g->col[index2]) {
                    ++index2; 
                } else {
                    ++commonNeighbors; ++index1; ++index2;
                }
                if (commonNeighbors == k-2) {
                    break;
                }
            }
            if (commonNeighbors < (k-2)) {
                g->col[i] = UINT_MAX;
                done = 0;
            }
        }
    }
}


void addEdge(Graph * g, int v1, int v2, double weight, int index)
{

    if (g->nnzSize == 0)
    {
        g->row = (unsigned int *)malloc(BASE_SIZE * sizeof(int));
        g->col = (unsigned int *)malloc(BASE_SIZE * sizeof(int));
        g->nnzSize = BASE_SIZE;

    }
    if (index > g->nnzSize)
    {
        g->nnzSize *= 2;
        g->row = (unsigned int *)realloc(g->row, g->nnzSize * sizeof(int));
        g->col = (unsigned int *)realloc(g->col, g->nnzSize * sizeof(int));
    }
    g->row[index] = v1;
    g->col[index] = v2;
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

    //sort by col
    for (unsigned int i = 0 ; i < g->nnzSize; ++i) {
        for (unsigned int j = 0 ; j < g->nnzSize - i - 1; ++j) {
            if (g->col[j] > g->col[j+1]) {
                int rowTemp = g->row[j];
                int colTemp = g->col[j];

                g->row[j] = g->row[j+1];
                g->col[j] = g->col[j+1];
                
                g->row[j+1] = rowTemp;
                g->col[j+1] = colTemp;
                
            }
        }
    }


    // Histogram
    for (unsigned int i = 0; i <g->nnzSize; ++i)
    {
        unsigned int row = g->row[i];
        rowPtrs[row]++;
    }

    // Prefix sum row
    unsigned int sumRow = 0;
    for (unsigned int row = 0; row < numRows; ++row)
    {
        unsigned int val = rowPtrs[row];
        rowPtrs[row] = sumRow;
        sumRow += val;
    }
    rowPtrs[numRows] = sumRow;

    // Binning
    for (unsigned int index = 0; index < g->nnzSize; ++index)
    {
        unsigned int row = g->row[index];
        unsigned int i = rowPtrs[row]++;
        colIdxs[i] = g->col[index];
        rowIdxs[i] = g->row[index];
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
        if (componentSizes[i] > 1)
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
