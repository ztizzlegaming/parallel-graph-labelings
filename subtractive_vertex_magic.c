//CSC 350 - Parallel Computing and MAT 490/491 - Research on Graph Labelings
//
//Read in a graph from a text file or generate a graph of two cycles joined by
//a certain number of vertices and try to find subtractive vertex magic
//labelings of the graph.
//
//To compile:
//    g++ subtractive_vertex_magic.c -O3 -o vertex_magic
//
//To run:
//    ./vertex_magic 3 2
//    ./vertex_magic 3 1
//    ./vertex_magic 4 3
//    ...
//
//Jordan Turley, Jason Pinto, Matthew Ko

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <vector>

#define OUTPUT_FILE_LINE_SIZE 1000
#define DEVELOPMENT 1 //Print out the permutations and write to file
#define PRODUCTION 2 //Only write the permutations to the file

struct Graph readGraph(char* filename);
struct Graph generateGraph(int cycleSize, int connectingVertices);
long factorial(long n);
int* generatePermutation(long idx, int permSize);
void printArr(int arr[], int size, int newline);

//Struct for a graph, holding the number of vertices, edges, and the
//adjacency matrix for the graph
struct Graph {
	int vertices;
	int edges;
	int** graph;
};

int main(int argc, char *argv[]) {
	int cycleSize = 4;
	int connectingVertices = 2;
	if (argc == 3) {
		cycleSize = atoi(argv[1]);
		connectingVertices = atoi(argv[2]);
	}

	printf("Cycle size = %d\nConnecting vertices = %d\n", cycleSize, connectingVertices);

	//Get the graph
	//char* filename = (char *) "graph.txt";
	//struct Graph graph = readGraph(filename);
	struct Graph graph = generateGraph(cycleSize, connectingVertices);
	int vertices = graph.vertices;
	int edges = graph.edges;
	int** matrix = graph.graph;

	//Print out graph
	printf("%d, %d\n", vertices, edges);
	for (int i1 = 0; i1 < vertices; i1++) {
		for (int i2 = 0; i2 < vertices; i2++) {
			printf("%d ", matrix[i1][i2]);
		}
		printf("\n");
	}

	//Generate first permutation
	long permSize = vertices + edges;
	int permutation[permSize];
	for (int i1 = 0; i1 < permSize; i1++) {
		permutation[i1] = i1 + 1;
	}

	long numPermutations = factorial(permSize);

	printf("Permutations to check: (|V| + |E|) = (%d + %d)! = %ld! = %ld\n",vertices, edges, permSize, numPermutations);

	time_t start;
	time_t finish;

	start = time(NULL);

	std::vector<long> worksIdxs;
	std::vector<int> magicNumbers;

	//Loop over all permutations
	//Try to find valid subtractive vertex magic labeling
	int magicNumber, firstVertex, works, vertex, curVertexValue, vertexIn, edge, vertexOut;
	for (long permIdx = 0; permIdx < numPermutations; permIdx++) {
		magicNumber = 0;
		firstVertex = 1;
		works = 1;
		for (vertex = 0; vertex < vertices; vertex++) {
			//Calculate the value for this vertex
			curVertexValue = permutation[vertex];

			//Add the values for the edges in
			for (vertexIn = 0; vertexIn < vertices; vertexIn++) {
				edge = matrix[vertexIn][vertex];
				if (edge) {
					//Add in edge value
					curVertexValue += permutation[vertices + edge - 1];
				}
			}

			//Subtract the values for the edges out
			for (vertexOut = 0; vertexOut < vertices; vertexOut++) {
				edge = matrix[vertex][vertexOut];
				if (edge) {
					//Subtract edge value
					curVertexValue -= permutation[vertices + edge - 1];
				}
			}

			if (firstVertex) {
				//If this is the first vertex, set the magic number
				magicNumber = curVertexValue;
				firstVertex = 0;
			} else {
				//Check if this vertex value is the same as the magic number
				if (magicNumber != curVertexValue) {
					works = 0;
					break;
				}
			}
		}
		
		if (works) {
			worksIdxs.push_back(permIdx);
			magicNumbers.push_back(magicNumber);
		}

		std::next_permutation(permutation, permutation + permSize);
	}

	finish = time(NULL);

	double timeTaken = difftime(finish, start);
	printf("Time: %f seconds\n", timeTaken);

	//Output file to store permutations
	//Write the graph parameters to the first line
	char filename[OUTPUT_FILE_LINE_SIZE];
	sprintf(filename, "output_%d_%d.txt", cycleSize, connectingVertices);
	FILE* outputFile = fopen(filename, "w");
	char firstLine[OUTPUT_FILE_LINE_SIZE];
	sprintf(firstLine, "Graph: Cycle size = %d, connecting vertices = %d\n",
		cycleSize, connectingVertices);
	fputs(firstLine, outputFile);

	//Write the time taken to the file
	char timeLine[OUTPUT_FILE_LINE_SIZE];
	sprintf(timeLine, "Time taken: %f seconds\n", timeTaken);
	fputs(timeLine, outputFile);

	//Write the adjacency matrix to the output file
	char matrixLine[OUTPUT_FILE_LINE_SIZE] = "";
	for (int i1 = 0; i1 < vertices; i1++) {
		for (int i2 = 0; i2 < vertices; i2++) {
			int edge = matrix[i1][i2];
			char part[3];
			sprintf(part, "%d ", edge);
			strcat(matrixLine, part);
		}
		//Add \n to the end
		char end[2] = "\n";
		strcat(matrixLine, end);
	}
	//Write line to file
	fputs(matrixLine, outputFile);

	//Print out all the permutations
	printf("Num worked: %lu\n", worksIdxs.size());
	for (int i1 = 0; i1 < worksIdxs.size(); i1++) {
		long permIdx = worksIdxs[i1];
		int magicNumber = magicNumbers[i1];
		int* perm = generatePermutation(permIdx, permSize);
		printArr(perm, permSize, 0);
		printf(" Magic Number: %d\n", magicNumber);

		//Write the permutation to the output file
		char outputLine[OUTPUT_FILE_LINE_SIZE] = "";
		sprintf(outputLine, "%d: {", (i1 + 1));
		for (int i1 = 0; i1 < permSize; i1++) {
			char part[10]; //Shouldn't ever be bigger than 6, 10 to be save
			if (i1 != permSize - 1) {
				sprintf(part, "%d, ", perm[i1]);
			} else {
				sprintf(part, "%d", perm[i1]);
			}
			strcat(outputLine, part);
		}
		char end[25]; //Shouldn't be bigger than 21, 25 to be save
		sprintf(end, "} Magic Number: %d\n", magicNumber);
		strcat(outputLine, end);
		fputs(outputLine, outputFile);
		
		free(perm);
	}
	printf("\n");

	//Close the file
	if (outputFile != NULL) {
		fclose(outputFile);
	}

	//Free the dynamic array for the graph
	for (int i1 = 0; i1 < vertices; i1++) {
		free(matrix[i1]);
	}
	free(matrix);

	return 0;
}

/**
 * Read in a graph from a given text file.
 * The first line should contain the number of vertices in the graph.
 * Then the file should contain the adjacency matrix of the graph.
 * @param filename The name of the file containing the graph data
 * @return A struct of the number of vertices, edges, and the adjacency matrix
 */
struct Graph readGraph(char* filename) {
	FILE *file;
	file = fopen(filename, "r");

	int vertices;
	int edges = 0;

	//First line of input file is number of vertices
	int result = fscanf(file, "%d\n", &vertices);
	if (result != 1) {
		printf("An error occured getting graph size from input file.\n");
		exit(1);
	}

	printf("Result = %d\n", result);

	//Init adjacency matrix
	int** graph = (int **) malloc(vertices * sizeof(int *));

	char* res;

	//Loop over each vertex
	for (int i1 = 0; i1 < vertices; i1++) {
		graph[i1] = (int *) malloc(vertices * sizeof(int));

		//Read in the line
		char line[256]; //Hope this is big enough
		res = fgets(line, sizeof(line), file);
		if (res == NULL) {
			printf("An error occured reading line from input file.\n");
			exit(1);
		}

		//Split the line into each number
		char* part = strtok(line, " ");
		int c = 0;
		while (part) {
			char numC = part[0];
			int num = numC - '0';
			if (num == 0) {
				graph[i1][c] = 0;
			} else {
				edges += num;
				graph[i1][c] = edges;
			}
			part = strtok(NULL, " ");
			c++;
		}
	}

	fclose(file);

	struct Graph g;
	g.vertices = vertices;
	g.edges = edges;
	g.graph = graph;

	return g;
}

/**
 * Generates a graph of two cycles connected by some vertices
 * @param cycleSize The size of the cycles to be connected
 * @param connectingVertices The number of vertices the cycles are connected by
 * @return A struct of the number of vertices, edges, and the adjacency matrix
 */
struct Graph generateGraph(int cycleSize, int connectingVertices) {
	//Calculate number of vertices and generate adjacency matrix
	int vertices = 2 * cycleSize - connectingVertices;
	int** graph = (int **) malloc(vertices * sizeof(int *));
	for (int i1 = 0; i1 < vertices; i1++) {
		graph[i1] = (int *) malloc(vertices * sizeof(int));
		for (int i2 = 0; i2 < vertices; i2++) {
			graph[i1][i2] = 0;
		}
	}

	//Connect the vertices
	graph[0][cycleSize - 1] = 1;
	graph[cycleSize - connectingVertices][vertices - 1] = cycleSize - connectingVertices + 2;
	int edges = 2;
	for (int i1 = 1; i1 < vertices; i1++) {
		if (i1 < cycleSize - connectingVertices + 1) {
			graph[i1][i1 - 1] = i1 + 1;
		} else {
			graph[i1][i1 - 1] = i1 + 2;
		}
		edges++;
	}

	struct Graph g;
	g.vertices = vertices;
	g.edges = edges;
	g.graph = graph;

	return g;
}

/**
 * Recursively calculate the factorial of the given number
 * @param n The number to calculate the factorial of
 * @return n!
 */
long factorial(long n) {
	if (n == 0) {
		return 1;
	}
	if (n == 1 || n == 2) {
		return n;
	}

	return n * factorial(n - 1);
}

int* generatePermutation(long idx, int permSize) {
	int* permutation = (int *) malloc(permSize * sizeof(int));
	for (int i1 = 0; i1 < permSize; i1++) {
		permutation[i1] = i1 + 1;
	}
	if (idx != 0) {
		//Generate perm for idx using factoradic numbers
		int factoradicLength = 0;
		long factorials[permSize];
		long fac = 1;
		int count = 0;
		while (idx >= fac) {
			fac = factorial(factoradicLength);
			factorials[count] = fac;
			count++;
			factoradicLength++;
		}
		factoradicLength--;
		//printf("Factoradic length: %d\n", factoradicLength);

		long factoradic[permSize];
		for (int i1 = 0; i1 < permSize; i1++) {
			factoradic[i1] = 0;
		}

		//printf("Array inited\n");

		long cur = idx;
		count = 0;
		for (int i1 = factoradicLength - 1; i1 >= 0; i1--) {
			long fNum = cur / factorials[i1];
			//printf("fNum = %ld\n", fNum);
			factoradic[count + permSize - factoradicLength] = fNum;
			count++;
			cur = cur % factorials[i1];
		}

		int permutationTemp[permSize];
		permutationTemp[0] = factoradic[0] + 1;
		permutation[factoradic[0]] = -1;
		for (int i1 = 1; i1 < permSize; i1++) {
			int fNum = factoradic[i1];
			int permIdx = 0;
			int actualIdx = 0;
			while (permIdx != fNum || permutation[actualIdx] == -1) {
				if (permutation[actualIdx] != -1) {
					permIdx++;
				}
				actualIdx++;
			}
			permutationTemp[i1] = permutation[actualIdx];
			permutation[actualIdx] = -1;
		}

		for (int i1 = 0; i1 < permSize; i1++) {
			permutation[i1] = permutationTemp[i1];
		}
	}

	return permutation;
}

void printArr(int arr[], int size, int newline) {
	printf("{");
	for (int i1 = 0; i1 < size; i1++) {
		printf("%d", arr[i1]);
		if (i1 != size - 1) {
			printf(", ");
		}
	}
	printf("}");
	if (newline) {
		printf("\n");
	}
}
