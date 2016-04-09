#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#include <chrono>
#include <math.h>


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
// Prints the solution (total distance as well as path to get there)
// INPUTS: src - source node #
//         dest - destination node #
//         distance - total distance of path
//         path_info = array of how each node was gotten to
// MODIFIES: nothing
// RETURNS: nothing
void printSolution(int src, int dest, int distance, int *path_info);
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

// Implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
// INPUTS: graph - unordered map with key of source node and value of
//                 list of pairs in format (destination_node, arc distance)
//         src - What node we are starting at and thus computing 
//              shortest paths from
//         node_count - row & column size of graph
// MODIFIES: nothing
// RETURNS: nothing

void dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph, std::map<int, std::vector<std::pair<int, double> > >  &land_dist,
                int node_count, int src, int dest);

//reads in graph_coords, graph_coords will be 1-indexed first int is lng, 2nd is lat
void readInGraphCoordinates(std::vector<std::pair<int, int> > &graph_coords, int &node_count,
                        const std::string &file_name);

void print_all_coords(std::vector<std::pair<int,int> > &graph_coords);


//calculates the distance between two points based on lng lat coords
double distance_btw_coords(std::pair<int,int> coord_1, std::pair<int,int> coord_2);


void get_landmarks(int k, std::vector<std::pair<int,int> > &graph_coords,
    std::set<int> &landmarks);

void get_dist_btw_landmarks(std::vector<std::pair<int,int> > &graph_coords, 
                            std::map<int,std::vector<std::pair<int,int> > > &graph,
                            std::set<int> &landmarks, std::map<int,std::vector<std::pair<int, double> > >  &land_dist);

void print_land_mark_distances(std::map<int, std::vector<std::pair<int, double> > >  &land_dist);


void print_landmarks(std::set<int> &landmarks);

std::map<int,std::vector<std::pair<int,int> > > graph;
std::vector<std::pair<int,int> > graph_coords; //after reading in coordinates will be 1-indexed



class ALT_Class
{
    private:
        std::set<int> landmarks;
        static std::map<int, std::vector<std::pair<int, double> > >  land_dist;
        static int dest;

    public:
        // void set_dest(int input_dest)
        // {
        //     dest = input_dest;
        // }
        ALT_Class()
        {
            dest = 0;
        }
        void get_landmarks(int k);
        void get_dist_btw_landmarks();
        void print_landmarks();

        void print_land_mark_distances();
        // bool alt_comparator(std::pair<int, int> left, std::pair<int, int> right);

        void alt_alg(int node_count, int src, int dest);

        class alt_comparator
        {
            public:
                bool operator() (std::pair<int, int> left, std::pair<int, int> right)
                {
                    int left_node = left.first;
                    int right_node = right.first;
                    double left_est = left.second;
                    double right_est = right.second;

                    double max_left_heur = INT_MIN;
                    double max_right_heur = INT_MIN;

                    for (int i = 0; i < land_dist[left_node].size(); ++i)
                    {
                        int land_mark_node = land_dist[left_node][i].first;
                        int left_dist_land = land_dist[left_node][i].second;

                        int right_dist_land = land_dist[right_node][i].second;

                        int dest_dist_land = land_dist[dest][i].second;

                        int total_left = left_dist_land + dest_dist_land;
                        int total_right = right_dist_land + dest_dist_land;

                        if(total_left > max_left_heur)
                        {
                            max_left_heur = total_left;
                        }

                        if(total_right > max_right_heur)
                        {
                            max_right_heur = total_right;
                        }
                    }
                    return max_left_heur >= max_right_heur;
                }
        };

        private:
            std::priority_queue<std::pair<int,int>, 
            std::vector<std::pair<int, int> >, 
            alt_comparator > dist_node;




};


std::map<int, std::vector<std::pair<int, double> > > ALT_Class::land_dist;
int ALT_Class::dest;

int main(int argc, char** argv){
    if (argc < 3){
        cout << "ABORTING: Not enough command line arguments\n";
        cout << "Need Graph filename\n";
        return 1;
    }
    // RoadMap graph used to find shortest path
    int node_count = 0;
    std::string file_name_arcs = argv[1];
    std::string file_name_coords = argv[2];
  
    
    readInGraph(graph, node_count, file_name_arcs);
    readInGraphCoordinates(graph_coords, node_count, file_name_coords);
    // printMap(graph);
    // print_all_coords(graph_coords);
    ALT_Class alt_inst;

    // alt_inst.set_dest(45);

    alt_inst.get_landmarks(8);

    alt_inst.print_landmarks();

    alt_inst.get_dist_btw_landmarks();

    alt_inst.print_land_mark_distances();

    alt_inst.alt_alg(graph_coords.size(), 5, 8000);

    int src, dest;

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

void readInGraphCoordinates(std::vector<std::pair<int, int> > &graph_coords, int &node_count,
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
                int node_num, lng, lat;
                infile >> node_num >> lng >> lat;
                // node1 already exists in map so just insert next arc

                // node1 doesn't exist in map so create it and add arc
                std::pair <int,int> temp_pair;
                temp_pair = std::make_pair(lng, lat);
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

void ALT_Class::print_landmarks()
{
    cout << "Landmarks are:" << endl;
    for (auto it = landmarks.begin(); it != landmarks.end(); it++)
    {
        cout << *it << endl;
    }
}


void ALT_Class::get_landmarks(int k)
{
    int num_nodes = graph_coords.size();
    int start = num_nodes / 2;

    landmarks.insert(start);

    while(landmarks.size() < k)
    {
        int cur_furthest_node = -1;
        double cur_furthest_avg = INT_MIN;

        for(int j = 1; j < graph_coords.size(); ++j)
        {
            //if landmark has already been added
            if(landmarks.find(j) != landmarks.end())
            {
                continue;
            }

            double total_dist = 0;
            // cout << j << endl;
            for(auto it = landmarks.begin(); it != landmarks.end(); it++)
            {
                int cur_node = *it;
                total_dist += distance_btw_coords(graph_coords[cur_node], graph_coords[j]);
            }
            double average_dist = total_dist/landmarks.size();
            if(average_dist > cur_furthest_avg)
            {
                cout << "avg distance_to_node is: " << average_dist << endl;
                cout << "cur_furthest_dist is: " << cur_furthest_avg << endl;
                cur_furthest_avg = average_dist;
                cur_furthest_node = j;
            }
        }


        cout << "Node " << cur_furthest_node << " is being inserted" << endl;
        landmarks.insert(cur_furthest_node);
    }
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



void print_all_coords(std::vector<std::pair<int,int> > &graph_coords)
{
    for(int i = 1; i < graph_coords.size(); ++i)
    {
        cout << "Coordinates of node " << i << " is: " << graph_coords[i].first << ", " << graph_coords[i].second;
        cout << endl;
    }
}

void ALT_Class::print_land_mark_distances()
{
    for(int i = 1; i < land_dist.size(); ++i)
    {
        for(int j = 0; j < land_dist[i].size(); ++j)
        {
            int land_mark_node = land_dist[i][j].first;
            double dist_to_land = land_dist[i][j].second;
            cout << "Distance Between " << i << " and node " << land_mark_node
                << " is: " << dist_to_land << endl;
        }
    }
}

void printSolution(int src, int dest, int distance, int *path_info)
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
}


double distance_btw_coords(std::pair<int,int> coord_1, std::pair<int,int> coord_2)
{
    double x1 = coord_1.first;
    double x2 = coord_1.second;

    double y1 = coord_2.first;
    double y2 = coord_2.second;

    double x = x1 - x2;
    double y = y1 - y2;
    double dist;

    dist = pow(x,2)+pow(y,2);           //calculating distance by euclidean formula
    dist = sqrt(dist);                  //sqrt is function in math.h

    return dist;
}

void ALT_Class::get_dist_btw_landmarks()
{
    // land_dist.resize(1);
    // for(int i = 1; i < graph_coords.size(); ++i)
    // {
    //     std::vector<std::pair<int, double> > temp_vector;
    //     for (auto it = landmarks.begin(); it != landmarks.end(); it++)
    //     {
    //         int land_mark_node = *it;
    //         double dist_between = distance_btw_coords(graph_coords[land_mark_node], graph_coords[i]);
    //         temp_vector.push_back(std::make_pair(land_mark_node, dist_between));
    //     }
    //     land_dist.push_back(temp_vector);
    // }

    for (auto it = landmarks.begin(); it != landmarks.end(); it++)
    {
        int land_mark_node = *it;
        dijkstra(graph, land_dist, graph_coords.size(), land_mark_node, 8);
    }
}

bool pair_comparator(std::pair<int, int> left, std::pair<int, int> right)
{
    return left.second >= right.second;
}


void dijkstra(std::map<int,std::vector<std::pair<int,int> > > &graph, 
                std::map<int, std::vector<std::pair<int, double> > >  &land_dist,
                int node_count, int src, int dest)
{
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
    while (!dist_node.empty()){
        // Pick the minimum distance vertex not visited
        while (!dist_node.empty() && visited[dist_node.top().first]){
            dist_node.pop();
        }
        if (dist_node.empty()){
            cout << "ABORTING: No Possible path from " << src << " to " << dest << endl;
            exit(1);
        }
        u = dist_node.top().first;
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
                dist_node.push(std::make_pair(v.first, dist[v.first]));
            }
        }
    }
    printSolution(src, dest, dist[dest], path_info);

    for(int i = 1; i < node_count + 1; ++i)
    {
        land_dist[i].push_back(std::make_pair(src, dist[i]));
    }
}

// bool ALT_Class::alt_comparator(std::pair<int, int> left, std::pair<int, int> right)
// {
//     int left_node = left.first;
//     int right_node = right.first;
//     double left_est = left.second;
//     double right_est = right.second;

//     double max_left_heur = INT_MIN;
//     double max_right_heur = INT_MIN;

//     for (int i = 0; i < land_dist[left_node].size(); ++i)
//     {
//         int land_mark_node = land_dist[left_node][i].first;
//         int left_dist_land = land_dist[left_node][i].second;

//         int right_dist_land = land_dist[right_node][i].second;

//         int dest_dist_land = land_dist[dest][i].second;

//         int total_left = left_dist_land + dest_dist_land;
//         int total_right = right_dist_land + dest_dist_land;

//         if(total_left > max_left_heur)
//         {
//             max_left_heur = total_left;
//         }

//         if(total_right > max_right_heur)
//         {
//             max_right_heur = total_right;
//         }
//     }


//     return max_left_heur >= max_right_heur;
// }


void ALT_Class::alt_alg(int node_count, int src, int dest)
{
    // For timing purposes only
    typedef std::chrono::duration<int,std::milli> millisecs_t;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    // Distance priority queue. dist_node[0] will hold the node with the minimum distance

    this->dest = dest;

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
                dist_node.push(std::make_pair(v.first, dist[v.first]));
            }
        }
    }
    printSolution(src, dest, dist[dest], path_info);
    // Again for timing purposes only
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end-start));
    std::cout << "That took: " << duration.count() << " milliseconds.\n";

}


