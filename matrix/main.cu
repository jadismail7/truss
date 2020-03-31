
#include "graph.h"
#include <string>
#include <cstring>

int main(int argc, char **argv)
{
	std::string directory = "./graphs/";
	std::string filename = argc>2?argv[1]:"graph.txt";
	std::string name = directory + filename;
	int k = argc>3?atoi(argv[2]):3;
	Graph * g = (Graph *) malloc(sizeof(Graph));
	g->nnzSize = 0;
	char * graphName = new char[name.size() + 1];
	strcpy(graphName, name.c_str());

	int numRows = initGraph(g, graphName);
	createCSRFromCOO(g, numRows);

	truss_cpu(g, k);

	return 0;
}


int initGraph(Graph * g, char * filename)
{
    int v1, v2;
    // double w;
    int numRows = 0;

    FILE* myfile = fopen(filename, "r");
	int index = 0;
    if (myfile)
    {
        while (!feof(myfile))
        {
			fscanf(myfile, "%d %d", &v1, &v2);
			addEdge(g, v1, v2, 1, index++);
            int localMax = v1>v2?v1:v2;
            numRows = numRows>localMax?numRows:localMax;
        }
        fclose(myfile);
	}
	done(g, index);
    return numRows;
}