#include <stdio.h>
#include <stdlib.h>

#define BASE_SIZE 1500

int graphSize;
double ** matrix;


    void Graph(const int size)
    {
        graphSize = size;
        matrix = (double**) malloc(graphSize*sizeof(double*));
        for (unsigned int i = 0; i < graphSize; ++i)
        {
            matrix[i] = (double*) malloc(graphSize*sizeof(double));
            for (unsigned int j = 0; j < graphSize; ++j)
            {
                matrix[i][j] = 0;
            }
        }
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


    void addEdge(int v1, int v2, double weight)
    {
        matrix[v1][v2] = weight;
        matrix[v2][v1] = weight;
    }

    void removeEdge(int v1, int v2)
    {
        matrix[v1][v2] = 0;
        matrix[v2][v1] = 0;
    }

    int* DFSUtil(int v, int visited[], int* size)
    {
        int* visiting = (int*) malloc(BASE_SIZE*sizeof(int));
        
        visited[v] = 1;
        unsigned int i = 0;
        visiting[i] = v;
        ++i;
        for (int j = 0; j < graphSize; ++j)
        {
            if (!visited[j] && matrix[v][j] != 0)
            {
                int size2 = BASE_SIZE;
                int * temp = DFSUtil(j, visited, &size2);
                int index = i;
                for (i; i < index + size2; i++) {
                    if (i >= *size) {
                        *size*=2;
                        visiting = (int*) realloc(visiting, (*size)*sizeof(int*));
                    }
                    visiting[i] = temp[i-index];
                }
            }
        }
        *size = i;
        return visiting;
    }

    int** connectedComponents(int* returnSize, int* componentSizes)
    {
        int** dfs = (int**) malloc(BASE_SIZE*sizeof(int*));
        int index = 0;
        int size = BASE_SIZE;
        int *visited = (int*)malloc(size*sizeof(int));
        for (int v = 0; v < size; v++)
            visited[v] = 0;

        for (int v = 0; v < graphSize; v++)
        {
            if (visited[v] == 0)
            {
                if (index >= size) {
                    size*=2;
                    dfs = (int**)realloc(dfs, size*sizeof(int*));
                    componentSizes = (int*)realloc(dfs, size*sizeof(int));
                }
                int size = BASE_SIZE;
                dfs[index] = (DFSUtil(v, visited, &size));
                componentSizes[index] = size;
                index++;
            }
        }
        *returnSize = index;
        return dfs;
    }


    void truss(int k)
    {
        int done = 0;
        while (!done)
        {
            done = 1;
            for (int v1 = 0; v1 < graphSize; v1++)
            {
                for (int v2 = 0; v2 < graphSize; v2++)
                {
                    if (matrix[v1][v2] != 0)
                    {
                        if (intersect(matrix[v1], matrix[v2], graphSize) < k - 2)
                        {
                            removeEdge(v1, v2);
                            done = 0;
                        }
                    }
                }
            }
        }
    }
    void printTrussComponents(int k)
    {
        int size;
        int *componentSizes = (int*) malloc(BASE_SIZE*sizeof(int));
        int** cc = connectedComponents(&size, componentSizes);

        for (int i = 0; i < size; ++i)
        {
            if (componentSizes[i] >= k)
            {
                printf("[");
                for (int j = 0; j < componentSizes[i]; ++j)
                {
                    char* end = ", ";
                    if (j == (componentSizes[i] - 1))
                    {
                        end = "]\n";
                    }
                    printf("%d%s", cc[i][j], end);
                }
            }
        }
    }

