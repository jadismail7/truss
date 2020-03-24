
#include "utils.h"
#include <string.h>

int main(int argc, char **argv)
{
	char * directory = "./graphs/";
	char* filename = argv[1];
	char * name = (char *) malloc(strlen(filename) + strlen(directory) + 1);
    strcpy(name, directory);
    strcat(name, filename);
	int k = atoi(argv[2]);
	int size = getSize(name);
	Graph(size);
	initGraph(name);

	truss(k);
	printTrussComponents(k);
	return 0;
}