
#include "graph.h"
#include <string>
__global__ void truss_kernel(COO *coo, CSR *csr, int k, int *done)
{
    unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
    int v1 = coo->row[index];
    int v2 = coo->col[index];

    unsigned int commonNeighbors = 0;
    for (unsigned int i = csr->rowPtr[v1]; i < csr->rowPtr[v1 + 1]; ++i)
    {
        for (unsigned int j = csr->rowPtr[v2]; j < j < csr->rowPtr[v2 + 1]; ++j)
        {
            unsigned int neighbor1 = csr->col[i];
            if (neighbor1 == csr->col[j] && neighbor1 != INT_MAX)
            {
                ++commonNeighbors;
            }
        }
    }
    if (commonNeighbors <= k)
    {
        csr->col[index] = INT_MAX;
        coo->col[index] = INT_MAX;
        done[0] = 0;
    }

}

void truss_gpu(COO *coo, CSR *csr, int k)
{
    COO *coo_d;
    CSR *csr_d;
    int *done_d;

    cudaMalloc((void **)&coo_d, sizeof(COO));
    cudaMalloc((void **)&csr_d, sizeof(CSR));
    cudaMalloc((void **)&done_d, sizeof(int));

    cudaMemcpy(coo_d, coo, sizeof(COO), cudaMemcpyHostToDevice);
    cudaMemcpy(csr_d, csr, sizeof(CSR), cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();

    unsigned int numThreads = 1024;
    unsigned int numBlocks = (coo->cooSize + numThreads - 1) / numThreads;
    int *done = (int *)malloc(sizeof(int));
    done[0] = 0;

    while (!done[0])
    {

        done[0] = 1;
        cudaMemcpy(done_d, done, sizeof(int), cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        truss_kernel<<<numThreads, numBlocks>>>(coo_d, csr_d, k, done_d);
        cudaDeviceSynchronize();
        cudaMemcpy(done, done_d, sizeof(int), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        printf("Done?\t%d\n", done[0]);
        
    }

    cudaMemcpy(coo, coo_d, sizeof(COO), cudaMemcpyDeviceToHost);
    cudaMemcpy(csr, csr_d, sizeof(CSR), cudaMemcpyDeviceToHost);

    cudaFree(coo_d);
    cudaFree(csr_d);
    cudaFree(done_d);
}

int intersect(double *s1, double *s2, int size)
{
    int intersectionSize = 0;
    for (unsigned int i = 0; i < size; ++i)
    {
        if (s1[i] != 0 && s2[i] != 0)
        {
            ++intersectionSize;
        }
    }
    return intersectionSize;
}

void addEdge(COO *coo, int v1, int v2, double weight, int index)
{

    if (coo->cooSize == 0)
    {
        coo->row = (unsigned int *)malloc(BASE_SIZE * sizeof(int));
        coo->col = (unsigned int *)malloc(BASE_SIZE * sizeof(int));
        coo->cooSize = BASE_SIZE;
        
    }
    if (index > coo->cooSize)
    {
        coo->cooSize *= 2;
        coo->row = (unsigned int *)realloc(coo->row, coo->cooSize * sizeof(int));
        coo->col = (unsigned int *)realloc(coo->col, coo->cooSize * sizeof(int));
    }
    coo->row[index] = v1;
    coo->col[index] = v2;
}

void done(COO * coo, int size) {
    coo->cooSize = size;
}

CSR *createCSRfromCOO(COO *A, int numRows)
{
    // Allocate
    unsigned int *rowPtrs = (unsigned int *)calloc(numRows + 1, sizeof(unsigned int));
    unsigned int *colIdxs = (unsigned int *)malloc(A->cooSize * sizeof(unsigned int));

    // Histogram
    for (unsigned int i = 0; i < A->cooSize; ++i)
    {
        unsigned int row = A->row[i];
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
    for (unsigned int index = 0; index < A->cooSize; ++index)
    {
        unsigned int row = A->row[index];
        unsigned int i = rowPtrs[row]++;
        colIdxs[i] = A->col[index];
    }
    A->col = colIdxs;
    // Restore row pointers
    for (unsigned int row = numRows - 1; row > 0; --row)
    {
        rowPtrs[row] = rowPtrs[row - 1];
    }
    rowPtrs[0] = 0;

    CSR *csr = (CSR *)malloc(sizeof(CSR));
    csr->rowPtrSize = numRows;
    csr->colSize = A->cooSize;
    csr->rowPtr = rowPtrs;
    csr->col = colIdxs;

    return csr;
}

int *DFSUtil(COO *coo, CSR *csr, int v, int visited[], int *size)
{
    int *visiting = (int *)malloc(BASE_SIZE * sizeof(int));

    visited[v] = 1;
    
    visiting[0] = v;
    unsigned int i;
    for (int j = csr->rowPtr[v]; j < csr->rowPtr[v + 1]; ++j)
    {
        if (!visited[j])
        {
            int size2 = BASE_SIZE;
            int *temp = DFSUtil(coo, csr, j, visited, &size2);
            for (i = 1; i < 1 + size2; i++)
            {
                if (i >= *size)
                {
                    *size *= 2;
                    visiting = (int *)realloc(visiting, (*size) * sizeof(int *));
                }
                visiting[i] = temp[i - 1];
            }
        }
    }
    *size = i;
    return visiting;
}

int **connectedComponents(COO *coo, CSR *csr, int *returnSize, int *componentSizes)
{
    int **dfs = (int **)malloc(BASE_SIZE * sizeof(int *));
    int index = 0;
    int size = BASE_SIZE;
    int *visited = (int *)malloc(coo->cooSize * sizeof(int));
    for (int v = 0; v < size; v++)
        visited[v] = 0;

    for (int v = 0; v < coo->cooSize; v++)
    {
        if (visited[v] == 0)
        {
            if (index >= size)
            {
                size *= 2;
                dfs = (int **)realloc(dfs, size * sizeof(int *));
                componentSizes = (int *)realloc(dfs, size * sizeof(int));
            }
            int size = BASE_SIZE;
            dfs[index] = (DFSUtil(coo, csr, v, visited, &size));
            componentSizes[index] = size;
            index++;
        }
    }
    *returnSize = index;
    return dfs;
}
void printTrussComponents(COO *coo, CSR *csr, int k)
{
    int size;
    int *componentSizes = (int *)malloc(BASE_SIZE * sizeof(int));
    int **cc = connectedComponents(coo, csr, &size, componentSizes);

    for (int i = 0; i < size; ++i)
    {
        if (componentSizes[i] >= k)
        {
            printf("[");
            for (int j = 0; j < componentSizes[i]; ++j)
            {
                char end[4] = ", \0";
                if (j == (componentSizes[i] - 1))
                {
                    strcpy(end, "]\n\0");
                }
                printf("%d%s", cc[i][j], end);
            }
        }
    }
}
