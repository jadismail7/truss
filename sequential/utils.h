#include "graph.h"



int getSize(char* filename)
{
    int size = 0;
    int v1, v2;
    double w;
    FILE* myfile = fopen(filename, "rt");
    
    char line[100];
    while (fgets(line, 100, myfile) != NULL)
    {
        sscanf(line, "%d %d", &v1, &v2);
        int localMax = v1 > v2 ? v1 : v2;
        size = localMax > size ? localMax : size;
    }
    fclose(myfile);

    return size + 1;
}

void initGraph(char* filename)
{
    int v1, v2;
    double w;

    FILE* myfile = fopen(filename, "r");

    if (myfile)
    {
        while (!feof(myfile))
        {
            fscanf(myfile, "%d %d", &v1, &v2);
            addEdge(v1, v2, 1);
        }
        fclose(myfile);
    }
}