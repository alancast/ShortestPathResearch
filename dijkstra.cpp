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

// Reads in coordinate graph from input file
// INPUTS: graph_coords - vector of double pairs of coordinates
//         node_count - number of nodes
//         file_name - string name of file 
// MODIFIES: graph_coords - puts coordinates of node i in entry i
//           node_count - sets it equal to the number of lines in the file
//                which should be the number of nodes in the graph
// RETURNS: nothing
void readInGraphCoordinates(std::vector<std::pair<double, double> > &graph_coords, 
                            int &node_count, const std::string &file_name);

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
void printDistanceArray(int *dist, int node_count);

// Prints the solution (total distance as well as path to get there)
// INPUTS: src - source node #
//         dest - destination node #
//         distance - total distance of path
//         path_info = array of how each node was gotten to
// MODIFIES: nothing
// RETURNS: nothing
void printSolution(int src, int dest, int distance, int *path_info,
                    const std::string &outfile);

// Implements Dijkstra's single source shortest path algorithm
// for a graph represented using edge map representation
// INPUTS: graph - unordered map with key of source node and value of
//                 list of pairs in format (destination_node, arc distance)
//         src - What node we are starting at and thus computing 
//              shortest paths from
//         node_count - row & column size of graph
//         outfile - name of output file we are writing to
// MODIFIES: nothing
// RETURNS: nothing
void dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph,
                int node_count, int src, int dest, const std::string &outfile);

// Implements Bidirectional Dijkstra's single source shortest path algorithm
// for a graph represented using edge map representation
// This assumes that all roads are two directional (hopefully this won't be a problem)
// INPUTS: graph - unordered map with key of source node and value of
//                 list of pairs in format (destination_node, arc distance)
//         src - What node we are starting at and thus computing 
//              shortest paths from
//         node_count - row & column size of graph
//         outfile - name of output file we are writing to
// MODIFIES: nothing
// RETURNS: nothing
void bidirectional_dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph,
                int node_count, int src, int dest, const std::string &outfile);

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
std::vector<std::pair<double,double> > graph_coords;
int main(int argc, char** argv){
    if (argc < 4){
        cout << "ABORTING: Not enough command line arguments, need 3\n";
        cout << "<distance_graph> <coordinates_graph> <outfile>\n";
        return 1;
    }
    // RoadMap graph used to find shortest path
    int node_count = 0;
    std::string distance_graph = argv[1];
    std::string coordinates_graph = argv[2];
    std::string outfile = argv[3];
    std::map<int,std::vector<std::pair<int,int> > > graph;
    readInGraph(graph, node_count, distance_graph);
    readInGraphCoordinates(graph_coords, node_count, coordinates_graph);
    int src, dest;
    char temp_char = 'c';
    while (temp_char == 'c'){
        char which_algos_char;
        cout << "Do you want to do both Dijkstras or only one?\n";
        cout << "1 for 1 directional only, 2 for bidirectional only, b for both:";
        cin >> which_algos_char;
        cout << "Enter starting node: ";
        cin >> src;
        cout << "Enter ending node: ";
        cin >> dest;
        if (which_algos_char == '1' || which_algos_char == 'b'){
            cout << "Starting Dijkstras algorithm" << endl;
            dijkstra(graph, node_count, src, dest, outfile);
        }
        if (which_algos_char == '2' || which_algos_char == 'b'){
            cout << "Starting Bidirectional Dijkstras algorithm" << endl;
            bidirectional_dijkstra(graph, node_count, src, dest, outfile);
        }
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

void readInGraphCoordinates(std::vector<std::pair<double, double> > &graph_coords, int &node_count,
                        const std::string &file_name)
{
    cout << "Reading in graph from file: " << file_name << endl;
    char temp_char;
    std::string line;
    std::ifstream infile(file_name);
    node_count = 0;
    graph_coords.resize(1);
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
                // node_count = atoi(strs[4].c_str());
            }
            // Arc info line (node1 node2 distance)
            else if (temp_char == 'v'){
                node_count++;
                int node_num;
                double lng, lat;
                infile >> node_num >> lng >> lat;
                // node1 already exists in map so just insert next arc

                // node1 doesn't exist in map so create it and add arc
                std::pair <double,double> temp_pair;
                temp_pair = std::make_pair(lng/1000000.0, lat/1000000.0);
                graph_coords.push_back(temp_pair);
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

void dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph,
                int node_count, int src, int dest, const std::string &outfile)
{
    // For writing to output file
    std::ofstream output;
    output.open(outfile);
    output << "c Starting one directional Dijkstras.\n";
    // Output starting coordinate and destination coordinate
    output << "s " << graph_coords[src].first << " " << graph_coords[src].second << endl;
    output << "d " << graph_coords[dest].first << " " << graph_coords[dest].second << endl;
    // For timing purposes only
    typedef std::chrono::duration<int,std::milli> millisecs_t;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    // Distance priority queue. dist_node[0] will hold the node with the minimum distance
    std::priority_queue<std::pair<int,int>, 
            std::vector<std::pair<int, int> >, 
            std::function<bool(std::pair<int,int>, std::pair<int,int>)> > dist_node(pair_comparator);
    dist_node.push(std::make_pair(src, 0));
    // Distance array. dist[i] will hold the shortest distance from src to i
    int *dist = new int[node_count+1];
    // How you got to that node array. path_info[i] will hold what node got us to i
    int *path_info = new int[node_count+1];
    // visited[i] will be true if we have already computed the shortest path to it
    bool *visited = new bool[node_count+1];
    // Initialize all distances as INFINITE and visited[] as false
    for (int i = 0; i < node_count+1; i++){
        dist[i] = INT_MAX, visited[i] = false, path_info[i] = -1;
    }
    // Distance of source vertex from itself is always 0
    dist[src] = 0;
    int u = 0;
    // Find shortest path to the destination vertex
    while (!visited[dest]){
        // Pick the minimum distance vertex not visited
        while (!dist_node.empty() && visited[dist_node.top().first]){
            dist_node.pop();
        }
        if (dist_node.empty()){
            cout << "ABORTING: No Possible path from " << src << " to " << dest << endl;
            exit(1);
        }
        u = dist_node.top().first;
        output << "u " << graph_coords[u].first << " " << graph_coords[u].second << endl;
        dist_node.pop();
        // Mark the picked vertex as visited
        visited[u] = true;
        // Update dist value of the adjacent vertices of the picked vertex.
        for (int i = 0; i < graph[u].size(); ++i)
        {
            std::pair<int, int> v = graph[u][i];
            // Update dist[v.first] only if there is an edge from 
            // u to v, and total weight of path from src to   
            // v through u is smaller than current value of dist[v]
            if (dist[u]+v.second < dist[v.first]){
                dist[v.first] = dist[u]+v.second;
                path_info[v.first] = u;
                output << "a " << graph_coords[v.first].first << " " << graph_coords[v.first].second << endl;
                dist_node.push(std::make_pair(v.first, dist[v.first]));
            }
        }
    }
    output.close();
    printSolution(src, dest, dist[dest], path_info, outfile);
    // Again for timing purposes only
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end-start));
    std::cout << "That took: " << duration.count() << " milliseconds.\n";
}

void bidirectional_dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph,
                int node_count, int src, int dest, const std::string &outfile)
{
    // For writing to output file
    std::ofstream output;
    output.open(outfile);
    output << "c Starting bi-directional Dijkstras (bijkstras).\n";
    // Output starting coordinate and destination coordinate
    output << "s " << graph_coords[src].first << " " << graph_coords[src].second << endl;
    output << "d " << graph_coords[dest].first << " " << graph_coords[dest].second << endl;
    // For timing purposes only
    typedef std::chrono::duration<int,std::milli> millisecs_t;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
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
    int *final_path_info = new int[node_count+1];
    // visited[x][i] will be true if we have already computed the shortest path to it
    bool **visited = new bool*[2];
    visited[0] = new bool[node_count+1];
    visited[1] = new bool[node_count+1];
    // Initialize all distances as INFINITE and visited[x][i] as false
    for (int i = 0; i < node_count+1; i++){
        dist[0][i] = INT_MAX, visited[0][i] = false, path_info[0][i] = -1;
        dist[1][i] = INT_MAX, visited[1][i] = false, path_info[1][i] = -1;
        final_path_info[i] = -1;
    }
    // Distance of source vertex from itself is always 0
    dist[0][src] = 0;
    dist[1][dest] = 0;
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
    // What node was the crossover node
    int crossover_node = -1;
    while (!found_in_each){
        // Pick the minimum distance vertex not visited (from either queue)
        // From SRC to Dest
        while (!src_dist_node.empty() && visited[0][src_dist_node.top().first]){
            src_dist_node.pop();
        }
        // From dest to src
        while (!dest_dist_node.empty() && visited[1][dest_dist_node.top().first]){
            dest_dist_node.pop();
        }
        // Again this assumes that roads go both directions
        if (src_dist_node.empty() || dest_dist_node.empty()){
            // All possible paths have been exhausted and none was found
            cout << "ABORTING: No Possible path from " << src << " to " << dest << endl;
            exit(1);
        }
        // Pick from whichever is lower
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
        // Output updated node we are searching from
        output << "u " << graph_coords[u].first << " " << graph_coords[u].second << endl;
        // Mark the picked vertex as visited
        visited[dir][u] = true;
        // As soon as you have found one in each you know that you can't
        // Add anything else to the queue so just check through all them
        // Until you find the shortest path
        if (visited[!dir][u]){
            combined_dist = dist[!dir][u] + dist[dir][u];
            found_in_each = true;
            crossover_node = u;
            continue;
        }
        // Update dist value of the adjacent vertices of the picked vertex.
        for (int i = 0; i < graph[u].size(); ++i)
        {
            std::pair<int, int> v = graph[u][i];
            // Update dist[dir][v.first] only if there is an edge from 
            // u to v, and total weight of path from src to v through u is 
            // smaller than current value of dist[dir][v]
            if (dist[dir][u] + v.second < dist[dir][v.first]){
                dist[dir][v.first] = dist[dir][u]+v.second;
                path_info[dir][v.first] = u;
                // From SRC to Dest
                if (dir == 0){
                    src_dist_node.push(std::make_pair(v.first, dist[dir][v.first]));
                    // Output forward search
                    output << "f " << graph_coords[v.first].first << " " << graph_coords[v.first].second << endl;
                }
                // From dest to src
                else{
                    dest_dist_node.push(std::make_pair(v.first, dist[dir][v.first]));
                    // Output backward search
                    output << "b " << graph_coords[v.first].first << " " << graph_coords[v.first].second << endl;
                }
            }
        }
    }
    // Find the shortest path from whatever is remaining in the queues
    while (!src_dist_node.empty()){
        int node_num = src_dist_node.top().first;
        src_dist_node.pop();
        // This path can't be shorter because both distances are longer
        if (!visited[1][node_num]){
            continue;
        }
        // No need to evaluate rest of queue cuz it won't produce a shorter path
        // Assuming all positive edge weights 
        if (dist[0][node_num] > combined_dist){
            break;
        }
        if (dist[1][node_num] != INT_MAX){
            if (dist[0][node_num] + dist[1][node_num] < combined_dist){
                combined_dist = dist[0][node_num] + dist[1][node_num];
                crossover_node = node_num;
            }
        }
    }
    while (!dest_dist_node.empty()){
        int node_num = dest_dist_node.top().first;
        dest_dist_node.pop();
        // This path can't be shorter because both distances are longer
        if (!visited[0][node_num]){
            continue;
        }
        // No need to evaluate rest of queue cuz it won't produce a shorter path
        // Assuming all positive edge weights 
        if (dist[1][node_num] > combined_dist){
            break;
        }
        if (dist[0][node_num] != INT_MAX){
            if (dist[0][node_num] + dist[1][node_num] < combined_dist){
                combined_dist = dist[0][node_num] + dist[1][node_num];
                crossover_node = node_num;
            }
        }
    }
    // Update final_path_info with paths from both ends
    int path_finder = crossover_node;
    // Get from u to src
    while (path_info[0][path_finder] != -1){
        final_path_info[path_finder] = path_info[0][path_finder];
        path_finder = path_info[0][path_finder];
    }
    // Get from dest to u
    path_finder = crossover_node;
    while (path_info[1][path_finder] != -1){
        final_path_info[path_info[1][path_finder]] = path_finder;
        path_finder = path_info[1][path_finder];
    }
    output.close();
    printSolution(src, dest, combined_dist, final_path_info, outfile);
    // Again for timing purposes only
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end-start));
    std::cout << "That took: " << duration.count() << " milliseconds.\n";
}

bool pair_comparator(std::pair<int, int> left, std::pair<int, int> right)
{
    return left.second >= right.second;
}

void printDistanceArray(int *dist, int node_count)
{
    cout<<"Vertex\tDistance from Source\n";
    for (int i = 0; i < node_count; i++)
        cout<<i<<"\t\t"<<dist[i]<<endl;
}

void printSolution(int src, int dest, int distance, int *path_info,
                    const std::string &outfile)
{
    // Print the shortest path from src to dest
    cout << "The shortest path from node " << src << " to node " << dest;
    cout << " is " << distance << endl;
    // cout << "Working backwards the path was: " << dest << " ";
    // int path_finder = dest;
    // while (path_info[path_finder] != src){
    //     cout << path_info[path_finder] << " ";
    //     path_finder = path_info[path_finder];
    // }
    // cout << src << endl;
    std::ofstream output;
    output.open(outfile, std::ios::app);
    output << "c Starting output of actual shortest path (from dest to source).\n";
    // Output destination node
    output << "n " << graph_coords[dest].first << " " << graph_coords[dest].second << endl;
    int path_finder = dest;
    // Output all path nodes
    while (path_info[path_finder] != src){
        int node_num = path_info[path_finder];
        output << "n " << graph_coords[node_num].first << " " << graph_coords[node_num].second << endl;
        path_finder = path_info[path_finder];
    }
    // Output source node
    output << "n " << graph_coords[src].first << " " << graph_coords[src].second << endl;
    output.close();
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









