#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <map>
#include <boost/algorithm/string.hpp>
#include <chrono>

using std::cout;
using std::cin;
using std::endl;

// Reads in graph from input file
// INPUTS: graph - unordered map with key of source node and value of
//                 list of pairs in format (destination_node, arc distance)
//         node_count - initialized to 0 
//         file_name - string name of file 
// MODIFIES: graph - updates it to have edges from each node 
//           node_count - sets it equal to the number of lines in the file
//                which should be the number of nodes in the graph
// RETURNS: nothing
void readInGraph(std::map<int,std::vector<std::pair<int,int> > > &graph, 
                    int &node_count, const std::string &file_name);

// ARCHIVED FUNCTION
// ----------------------------------------------------------------------
// A debugging function that prints out the graph
// INPUTS: graph - 2D array by double pointer of graph
//         node_count - row & column size of graph
// MODIFIES: nothing
// RETURNS: nothing
// void printGraph(int **graph, int node_count);
// ----------------------------------------------------------------------

// A debugging function that prints out the graph
// INPUTS: graph - unordered map of <int, vector<pairs>>
// MODIFIES: nothing
// RETURNS: nothing
void printMap(std::map<int,std::vector<std::pair<int,int> > > &graph);
// Prints the constructed distance array from starting point
// INPUTS: dist - int array with distances to every other node
//         node_count = length of dist array
// MODIFIES: nothing
// RETURNS: nothing
void printSolution(int *dist, int node_count);
// Implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
// INPUTS: graph - unordered map with key of source node and value of
//                 list of pairs in format (destination_node, arc distance)
//         src - What node we are starting at and thus computing 
//              shortest paths from
//         node_count - row & column size of graph
// MODIFIES: nothing
// RETURNS: nothing
void bidirectional_dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph,
                int node_count, int src, int dest);
// Custom comparator for priority queue
// INPUTS: left - pair of (node, dist) of left part
//         right - pair of (node, dist) of right part
// MODIFIES: nothing
// RETURNS: true if left >= right
bool pair_comparator(std::pair<int, int> left, std::pair<int, int> right);

// The main driver of the program
// Finds the shortest path from a selected node to all others
// COMMAND LINE ARGS: currently none, will be input file at some point
// BEGINNING MAIN
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
int main(int argc, char** argv){
    if (argc < 2){
        cout << "ABORTING: Not enough command line arguments\n";
        cout << "Need Graph filename\n";
        return 1;
    }
    // RoadMap graph used to find shortest path
    int node_count = 0;
    std::string file_name = argv[1];
    std::map<int,std::vector<std::pair<int,int> > > graph;
    readInGraph(graph, node_count, file_name);
    int src, dest;
    char temp_char = 'c';
    while (temp_char == 'c'){
        cout << "Enter starting node: ";
        cin >> src;
        cout << "Enter ending node: ";
        cin >> dest;
        cout << "Starting Dijkstras algorithm" << endl;
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        bidirectional_dijkstra(graph, node_count, src, dest);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        typedef std::chrono::duration<int,std::milli> millisecs_t ;
        millisecs_t duration( std::chrono::duration_cast<millisecs_t>(end-start));
        std::cout << "That took: " << duration.count() << " milliseconds.\n";
        cout << "Do you want to do another pair?\n";
        cout << "c to continue, q to quit:";
        cin >> temp_char;
    }
    return 0;
}
// ENDING MAIN
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

void readInGraph(std::map<int,std::vector<std::pair<int,int> > > &graph, 
                    int &node_count, const std::string &file_name)
{
    cout << "Reading in graph from file: " << file_name << endl;
    char temp_char;
    std::string line;
    std::ifstream infile(file_name);
    if (infile.is_open()){
        // Read whole file
        while (!infile.eof()){
            // Get first character of line
            infile >> temp_char;
            // Text line
            if (temp_char == 'c'){
                // waste the line
                getline(infile,line);
                // cout << line << endl;
            }
            // Node count and arc count line
            else if (temp_char == 'p'){
                getline(infile,line);
                // cout << line << endl;
                // currently line is in format (" sp node_count arc_cout")
                std::vector<std::string> strs;
                boost::split(strs, line, boost::is_any_of(" "));
                node_count = atoi(strs[2].c_str());
            }
            // Arc info line (node1 node2 distance)
            else if (temp_char == 'a'){
                int node1, node2, distance;
                infile >> node1 >> node2 >> distance;
                // node1 already exists in map so just insert next arc
                if (graph.find(node1) != graph.end()){
                    graph[node1].push_back(std::make_pair(node2, distance));
                }
                // node1 doesn't exist in map so create it and add arc
                else{
                    std::pair <int,int> temp_pair;
                    temp_pair = std::make_pair(node2, distance);
                    std::vector<std::pair<int, int> > temp_vector;
                    temp_vector.push_back(temp_pair);
                    graph[node1] = temp_vector;
                }
            }
            // Some unkown starting character
            else{
                cout << "ABORTING: Saw unknown character in file\n";
                cout << "Character was:\t" << temp_char << endl;
                exit(2);
            }
        }
        infile.close();
        return;
    }
    else{
        cout << "ABORTING: Unable to open input file\n";
        exit(1);
    }
}

/*
ARCHIVED FUNCTION
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
*/

void bidirectional_dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph,
                int node_count, int src, int dest)
{
    // Distance priority queue. src_dist_node[0] will hold the node with the minimum distance
    // This one is from src to destination
    std::priority_queue<std::pair<int,int>, 
            std::vector<std::pair<int, int> >, 
            std::function<bool(std::pair<int,int>, std::pair<int,int>)> > src_dist_node(pair_comparator);
    src_dist_node.push(std::make_pair(src, 0));
    // Distance priority queue. dest_dist_node[0] will hold the node with the minimum distance
    // This one is from destination to src
    std::priority_queue<std::pair<int,int>, 
            std::vector<std::pair<int, int> >, 
            std::function<bool(std::pair<int,int>, std::pair<int,int>)> > dest_dist_node(pair_comparator);
    dest_dist_node.push(std::make_pair(dest, 0));
    // Distance array. dist[x][i] will hold the shortest distance from start to i
    int **dist = new int*[2];
    dist[0] = new int[node_count+1];
    dist[1] = new int[node_count+1];
    // How you got to that node array. path_info[x][i] will hold what node got us to i
    int **path_info = new int*[2];
    path_info[0] = new int[node_count+1];
    path_info[1] = new int[node_count+1];
    // visited[x][i] will be true if we have already computed the shortest path to it
    bool **visited = new bool*[2];
    visited[0] = new bool[node_count+1];
    visited[1] = new bool[node_count+1];
    // Initialize all distances as INFINITE and visited[x][i] as false
    for (int i = 0; i < node_count+1; i++){
        dist[0][i] = INT_MAX, visited[0][i] = false, path_info[0][i] = -1;
        dist[1][i] = INT_MAX, visited[1][i] = false, path_info[1][i] = -1;
    }
    // Distance of source vertex from itself is always 0
    dist[0][src] = 0;
    dist[1][dest] = 0;
    // Find shortest path to the destination vertex
    // Will be set to true when overlapped distance < all new paths
    bool src_finished = false;
    bool dest_finished = false;
    // Will be 0 when going from src to dest
    // Will be 1 when going from dest to src
    int dir = 0;
    // Shortest known distance of path 
    // Combination of dest to point and source to point
    int combined_dist = INT_MAX;
    // Used for vertex we are analyzing in the loop
    int u = 0;
    // Will be set to true if we "visited" a node coming from each direction
    bool found_in_each = false;
    while (!src_finished || !dest_finished){
        // As soon as you have found one in each you know that you can't
        // Add anything else to the queue so just check through all them
        // Until you find the shortest path
        if (found_in_each){
            while (!src_dist_node.empty()){
                int node_num = src_dist_node.top().first;
                if (dist[1][node_num] != INT_MAX){
                    if (dist[0][node_num] + dist[1][node_num] < combined_dist){
                        combined_dist = dist[0][node_num] + dist[1][node_num];
                    }
                }
                src_dist_node.pop();
            }
            while (!dest_dist_node.empty()){
                int node_num = dest_dist_node.top().first;
                if (dist[0][node_num] != INT_MAX){
                    if (dist[0][node_num] + dist[1][node_num] < combined_dist){
                        combined_dist = dist[0][node_num] + dist[1][node_num];
                    }
                }
                dest_dist_node.pop();
            }
            src_finished = true;
            dest_finished = true;
            continue;
        }
        // Pick the minimum distance vertex not visited (from either queue)
        // From SRC to Dest
        while (!src_dist_node.empty() && 
            (visited[0][src_dist_node.top().first] 
                || src_dist_node.top().second > combined_dist)){
            src_dist_node.pop();
        }
        // From dest to src
        while (!dest_dist_node.empty() && 
            (visited[1][dest_dist_node.top().first] 
                || dest_dist_node.top().second > combined_dist)){
            dest_dist_node.pop();
        }
        if (src_dist_node.empty()){
            // All possible paths from src have been exhausted
            if (src_finished){
                dir = 1;
            }
            else if (combined_dist != INT_MAX){
                cout << "SRC searching is done!" << endl;
                src_finished = true;
                dir = 1;
            }
            else{
                cout << "ABORTING: No Possible path from " << src << " to " << dest << endl;
                exit(1);
            }
        }
        if (dest_dist_node.empty()){
            // All possible paths from dest have been exhausted
            if (dest_finished){
                dir = 0;
            }
            else if (combined_dist != INT_MAX){
                cout << "Dest searching is done!" << endl;
                dest_finished = true;
                dir = 0;
            }
            else{
                cout << "ABORTING: No Possible path from " << dest << " to " << src << endl;
                exit(1);
            }
        }
        // Figure out which scenario we are in (what queue to pop from)
        // Pick from whichever is lower (this will happen 99% of time)
        if (!src_finished && !dest_finished){
            // Source has lower distance than destination so choose from it
            if (src_dist_node.top().second <= dest_dist_node.top().second){
                u = src_dist_node.top().first;
                src_dist_node.pop();
                dir = 0;
            }
            // Destination has lower distance than source so choose from it
            else{
                u = dest_dist_node.top().first;
                dest_dist_node.pop();
                dir = 1;
            }
        }
        // Dest is finished so only go from source
        else if (!src_finished && dest_finished){
            u = src_dist_node.top().first;
            src_dist_node.pop();
            dir = 0;
        }
        // Source is finished so only go from dest
        else if (src_finished && !dest_finished){
            u = dest_dist_node.top().first;
            dest_dist_node.pop();
            dir = 1;
        }
        // both are done so continue
        else{
            continue;
        }
        // Mark the picked vertex as visited
        visited[dir][u] = true;
        // As soon as you have found one in each you know that you can't
        // Add anything else to the queue so just check through all them
        // Until you find the shortest path
        if (visited[!dir][u]){
            if (dist[!dir][u] + dist[dir][u] < combined_dist){
                combined_dist = dist[!dir][u] + dist[dir][u];
            }
            found_in_each = true;
            continue;
        }
        // Update dist value of the adjacent vertices of the picked vertex.
        for (int i = 0; i < graph[u].size(); ++i)
        {
            std::pair<int, int> v = graph[u][i];
            // Update dist[dir][v.first] only if there is an edge from 
            // u to v, and total weight of path from src to v through u is 
            // smaller than current value of dist[dir][v]
            if (dist[dir][u] + v.second < dist[dir][v.first]
                && dist[dir][u] + v.second < combined_dist){
                dist[dir][v.first] = dist[dir][u]+v.second;
                path_info[dir][v.first] = u;
                // Potentially reducing combined distance
                if (dist[!dir][v.first] != INT_MAX &&
                    (dist[dir][v.first] + dist[!dir][v.first])<combined_dist){
                    combined_dist = dist[dir][v.first] + dist[!dir][v.first];
                }
                // From SRC to Dest
                if (dir == 0){
                    src_dist_node.push(std::make_pair(v.first, dist[dir][v.first]));
                }
                // From dest to src
                else{
                    dest_dist_node.push(std::make_pair(v.first, dist[dir][v.first]));
                }
            }
        }
    }
    // Print the constructed distance array
    // printSolution(dist, node_count);
    // Print the shortest path from src to dest
    cout << "The shortest path from node " << src << " to node " << dest;
    cout << " is " << combined_dist << endl;
    // cout << "Working backwards the path was: " << dest << " ";
    // int path_finder = dest;
    // while (path_info[path_finder] != src){
    //     cout << path_info[path_finder] << " ";
    //     path_finder = path_info[path_finder];
    // }
    // cout << src << endl;
}

bool pair_comparator(std::pair<int, int> left, std::pair<int, int> right)
{
    return left.second >= right.second;
}

void printSolution(int *dist, int node_count)
{
    cout<<"Vertex\tDistance from Source\n";
    for (int i = 0; i < node_count; i++)
        cout<<i<<"\t\t"<<dist[i]<<endl;
}

void printMap(std::map<int,std::vector<std::pair<int,int> > > &graph)
{
    for (std::map<int,std::vector<std::pair<int,int> > >::iterator it=graph.begin(); 
                    it!=graph.end(); ++it){
        cout << "Outgoing nodes for node " << it->first << " include: ";
        for (std::vector<std::pair<int,int> >::iterator it2=it->second.begin(); 
                it2!=it->second.end(); ++it2)
        {
            cout << "(" << it2->first << ", " << it2->second << "), ";
        }
        cout << endl;
    }
}









