
#include "graph.h"
#include <string>
#include <cstring>

int main(int argc, char **argv)
{
	printf("Is this working?");
	std::string directory = "./graphs/";
	std::string filename = argc>2?argv[1]:"graph.txt";
	std::string name = directory + filename;
	int k = argc>3?atoi(argv[2]):3;
	COO * coo = (COO*) malloc(sizeof(COO));
	coo->cooSize = 0;
	char * graphName = new char[name.size() + 1];
	strcpy(graphName, name.c_str());

	int numRows = initGraph(coo, graphName);
	printf("Numrows = %d\n", numRows);
	CSR * csr = createCSRfromCOO(coo, numRows);

	truss_gpu(coo, csr, k);
	printTrussComponents(coo, csr, k);
	return 0;
}


int initGraph(COO* coo, char * filename)
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
            addEdge(coo, v1, v2, 1, index++);

            int localMax = v1>v2?v1:v2;
            numRows = numRows>localMax?numRows:localMax;
        }
        fclose(myfile);
	}
	done(coo, index);
    return numRows;
}