
#include "graph.h"
#include "../timer.h"
#include <string>
__global__ void truss_kernel(Graph g, int k, int *done)
{
   
}

void copyGraph(Graph * source, Graph * destination, cudaMemcpyKind direction) {
   
}
void truss_gpu(Graph * g, int k)
{
    Graph graph_d;
    cudaMalloc((void**)&graph_d.row, g->nnzSize*sizeof(int));
    cudaMalloc((void**)&graph_d.col, g->nnzSize*sizeof(int));
    cudaMalloc((void**)&graph_d.rowPtr, g->rowPtrSize*sizeof(int));
    
    cudaMemcpy(graph_d.row, g->row, g->nnzSize*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(graph_d.col, g->col, g->nnzSize*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(graph_d.rowPtr, g->rowPtr, g->rowPtrSize*sizeof(int), cudaMemcpyHostToDevice);
    
    int *done_d;
    graph_d.nnzSize = g->nnzSize;
    graph_d.rowPtrSize = g->rowPtrSize;

    cudaMalloc((void **)&done_d, sizeof(int));


    cudaDeviceSynchronize();
    unsigned int numThreads = 1024;
    unsigned int numBlocks = (g->nnzSize + numThreads - 1) / numThreads;
    int *done = (int *)malloc(sizeof(int));
    *done = 0;

    Timer timer;
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

    cudaMemcpy(g->row, graph_d.row, g->nnzSize*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(g->col, graph_d.col, g->nnzSize*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(g->rowPtr, graph_d.rowPtr, g->rowPtrSize*sizeof(int), cudaMemcpyDeviceToHost);
    
    cudaFree(graph_d.row);
    cudaFree(graph_d.col);
    cudaFree(graph_d.rowPtr);
    cudaFree(&graph_d);
    cudaFree(done_d);
}

Graph * mult_cpu(Graph * g, int k) {
    Graph * s = (Graph *) malloc(sizeof(Graph));
    s->nnzSize = 0;
    unsigned int index = 0;
    for (unsigned int i = 0 ; i < g->rowPtrSize; ++i) {
        for (unsigned int j = 0 ; j < g->rowPtrSize; ++j) {
            int sum = 0;
            unsigned int m1 = g->rowPtr[i];
            unsigned int m2 = g->rowPtr[j];
            while (m1 < g->rowPtr[i+1] && m2 < g->rowPtr[j+1]) {
                if (g->col[m1] > g->col[m2]) {
                    ++m2;
                } else if (g->col[m2] > g->col[m1]) {
                    ++m1;
                } else {
                    ++sum; ++m1; ++m2;
                }
            }
            if (sum != 0) {
                addEdge(s, i, j, sum, index++);
            }
        }
    }
    done(s, index);
    return s;
}

Graph * elementWise_mult_cpu(Graph * g1, Graph * g2, int * rowPtrSize) {
    Graph * s = (Graph*) malloc(sizeof(Graph));
    s->nnzSize = 0;
    unsigned int sIndex = 0;
    unsigned int index1 = 0;
    unsigned int index2 = 0;

    while (index1 < g1->nnzSize && index2 < g2->nnzSize) {
        if (g1->row[index1] < g2->row[index2]) {
            ++index1;
        } else if (g1->row[index1] > g2->row[index2]) {
            ++index2;
        } else {
            if (g1->col[index1] < g2->col[index2]) {
                ++index1;
            } else if (g1->col[index1] > g2->col[index2]) {
                ++index2;
            } else {
                if (rowPtrSize != NULL) 
                    *rowPtrSize = g1->row[index1];
                    
                addEdge(s, g1->row[index1], g1->col[index1], g1->values[index1]*g2->values[index2], sIndex++);
                ++index1; ++index2;
            }
        }
    }
    done(s, sIndex);
    return s;
}

Graph * truss_cpu(Graph * g, int k) {
    int changed = true;
    while(changed) {
        changed = false;
        Graph * s = mult_cpu(g, k);
        
        s = elementWise_mult_cpu(s, g, NULL);
        
        
        Graph * m = (Graph *) malloc(sizeof(Graph));
        m->nnzSize = 0;
        int mIndex = 0;
        for (unsigned int i = 0 ; i < s->nnzSize; ++i) {
            if (s->values[i] >= k-2) {
                addEdge(m, s->row[i], s->col[i], 1, mIndex++);
            } else {
                changed = true;
            }
        }
        done(m, mIndex);
        int rowPtrSize;
        g = elementWise_mult_cpu(g, m, &rowPtrSize);
        createCSRFromCOO(g, rowPtrSize);
        
    }
    return g;

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

void sortGraphByCol(Graph * g) {
    //sort by col
    for (unsigned int i = 0 ; i < g->nnzSize; ++i) {
        for (unsigned int j = 0 ; j < g->nnzSize - i - 1; ++j) {
            if (g->col[j] > g->col[j+1]) {
                int rowTemp = g->row[j];
                int colTemp = g->col[j];
                int valTemp = g->values[j];

                g->row[j] = g->row[j+1];
                g->col[j] = g->col[j+1];
                g->values[j] = g->values[j+1];
                
                g->row[j+1] = rowTemp;
                g->col[j+1] = colTemp;
                g->values[j+1] = valTemp;
            }
        }
    }
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





//-------------------------Graph Stuff----------------------------------------------------------------------------------------------------------------

unsigned int *DFSUtil(Graph * g, int v, int visited[], int *size)
{

    unsigned int *visiting = (unsigned int *)malloc(BASE_SIZE * sizeof(int));

    visited[v] = 1;

    unsigned int i = 0;
    visiting[i++] = v;
    for (int j = g->rowPtr[v]; j < g->rowPtr[v + 1]; ++j)
    {
        if (!visited[g->col[j]] )
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
