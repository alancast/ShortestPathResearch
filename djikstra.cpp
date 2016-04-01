#include <iostream>
#include <fstream>
#include <string>

using std::cout;
using std::endl;

// Reads in graph from input file
// INPUTS: node_count - initialized to 0 
//         file_name - string name of file 
// MODIFIES: node_count - sets it equal to the number of lines in the file
//                which should be the number of nodes in the graph
// RETURNS: double pointer to 2D array of adjacency graph
int** readInGraph(int &node_count, const std::string &file_name);
// A debugging function that prints out the graph
// INPUTS: graph - 2D array by double pointer of graph
//         node_count - row & column size of graph
// MODIFIES: nothing
// RETURNS: nothing
void printGraph(int **graph, int node_count);
// Prints the constructed distance array from starting point
// INPUTS: dist - int array with distances to every other node
//         node_count = length of dist array
// MODIFIES: nothing
// RETURNS: nothing
void printSolution(int *dist, int node_count);
// Implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
// INPUTS: graph - 2D array by double pointer of graph
//         src - What node we are starting at and thus computing 
//              shortest paths from
//         node_count - row & column size of graph
// MODIFIES: nothing
// RETURNS: nothing
void dijkstra(int **graph, int src, int node_count);
// Finds the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
// INPUTS: dist - int array with distances to every other node
//         sptSet - bool array with same dimensions as dist
//                  represents nodes that have been finished
//                  True if they have, false if they haven't
//         node_count = length of dist array
// MODIFIES: nothing
// RETURNS: the index of the vertex with the minimum distance
int minDistance(int *dist, bool *sptSet, int node_count);

// The main driver of the program
// Finds the shortest path from a selected node to all others
// COMMAND LINE ARGS: currently none, will be input file at some point
int main(int argc, char** argv){
    if (argc < 2){
        cout << "ABORTING: Not enough command line arguments\n";
        cout << "Need Graph filename\n";
        return 1;
    }
    // RoadMap graph used to find shortest path
    int node_count = 0;
    std::string file_name = argv[1];
    int **graph = readInGraph(node_count, file_name);
    printGraph(graph, node_count);
    dijkstra(graph, 0, node_count);
    return 0;
}

int** readInGraph(int &node_count, const std::string &file_name)
{
    std::ifstream infile(file_name);
    if (infile.is_open()){
        infile >> node_count;
        int **graph = new int*[node_count];
        for(int i = 0; i < node_count; ++i){
            graph[i] = new int[node_count];
        }
        for (int i = 0; i < node_count; ++i)
        {
            for (int j = 0; j < node_count; ++j)
            {
                infile >> graph[i][j];
            }
        }
        infile.close();
        return graph;
    }
    else{
        cout << "ABORTING: Unable to open input file\n";
        exit(1);
    }
}


void printGraph(int **graph, int node_count)
{
    for (int i = 0; i < node_count; ++i)
    {
        for (int j = 0; j < node_count; ++j)
        {
            cout << graph[i][j] << " ";
        }
        cout << endl;
    }
}

void printSolution(int *dist, int node_count)
{
    cout<<"Vertex\tDistance from Source\n";
    for (int i = 0; i < node_count; i++)
        cout<<i<<"\t\t"<<dist[i]<<endl;
}

void dijkstra(int **graph, int src, int node_count)
{
    // Distance array. dist[i] will hold the shortest distance from src to i
    int *dist = new int[node_count];
 
    // sptSet[i] will be true if vertex i is included in shortest
    // path tree or shortest distance from src to i is finalized
    bool *sptSet = new bool[node_count];
 
    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0; i < node_count; i++){
        dist[i] = INT_MAX, sptSet[i] = false;
    }
    // Distance of source vertex from itself is always 0
    dist[src] = 0;
    // Find shortest path for all vertices
    for (int i = 0; i < node_count-1; i++)
    {
        // Pick the minimum distance vertex from the set of vertices not
        // yet processed. u is always equal to src in first iteration.
        int u = minDistance(dist, sptSet, node_count);
        // Mark the picked vertex as processed
        sptSet[u] = true;
        // Update dist value of the adjacent vertices of the picked vertex.
        for (int v = 0; v < node_count; v++){
            // Update dist[v] only if is not in sptSet, there is an edge from 
            // u to v, and total weight of path from src to  v through u is 
            // smaller than current value of dist[v]
            if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX 
                                       && dist[u]+graph[u][v] < dist[v]){
                dist[v] = dist[u] + graph[u][v];
            }
        }
    } 
    // print the constructed distance array
    printSolution(dist, node_count);
}

int minDistance(int *dist, bool *sptSet, int node_count)
{
    // Initialize min value
    int min = INT_MAX, min_index;
    for (int v = 0; v < node_count; v++){
        if (sptSet[v] == false && dist[v] <= min){
            min = dist[v], min_index = v;
        }
    }
   return min_index;
}










