#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <map>
#include <set>
#include <unordered_set>
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
        std::unordered_set<int> landmarks;
        static std::map<int, std::vector<std::pair<int, double> > >  land_dist;
        static int dest;
        // std::unordered_set<int> open_set;
        std::unordered_set<int> closed_set;
        std::unordered_set<int> closed_set_f;
        std::unordered_set<int> closed_set_b;
        int landmark_index = 0;
        int *g_score;
        int *f_score;
        int *path_info;

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

        void bi_alt_alg(int node_count, int src, int dest);

        double heuristic_cost_estimate(int start, int dest);

        void choose_landmark_index(int src, int dest);


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

                    for (int i = 0; i < 1; ++i)
                    {
                        int land_mark_node = land_dist[left_node][i].first;
                        int left_dist_land = land_dist[left_node][i].second;

                        int right_dist_land = land_dist[right_node][i].second;

                        int dest_dist_land = land_dist[dest][i].second;

                        int total_left = abs(left_dist_land - dest_dist_land);
                        int total_right = abs(right_dist_land + dest_dist_land);

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
            // std::priority_queue<std::pair<int,int>, 
            // std::vector<std::pair<int, int> >, 
            // alt_comparator > dist_node;



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

    alt_inst.get_landmarks(8);

    alt_inst.get_dist_btw_landmarks();

    int src, dest;

    char temp_char = 'c';
    while (temp_char == 'c'){
        char which_algos_char;
        cout << "Do you want to do both Dijkstras or only one?\n";
        cout << "1 for 1 directional only, 2 for bidirectional only, b for both:, 3 for ALT \n"
        <<"4 for bidirectional ALT: ";
        cin >> which_algos_char;
        cout << "Enter starting node: ";
        cin >> src;
        cout << "Enter ending node: ";
        cin >> dest;
        if (which_algos_char == '1' || which_algos_char == 'b'){
            cout << "Starting Dijkstras algorithm" << endl;
            // dijkstra(graph, node_count, src, dest);
        }
        if (which_algos_char == '2' || which_algos_char == 'b'){
            cout << "Starting Bidirectional Dijkstras algorithm" << endl;
            // bidirectional_dijkstra(graph, node_count, src, dest);
        }
        if (which_algos_char == '3')
        {
            cout << "Starting ALT algorithm" << endl;
            alt_inst.alt_alg(graph_coords.size(), src, dest);
            // alt_inst.bi_alt_alg(graph_coords.size(), src, dest);
        }
        if (which_algos_char == '4')
        {
            cout << "Starting Bi-Directional ALT algorithm" << endl;
            alt_inst.bi_alt_alg(graph_coords.size(), src, dest);
            // alt_inst.bi_alt_alg(graph_coords.size(), src, dest);
        }
        cout << "Do you want to do another pair?\n";
        cout << "c to continue, q to quit:";
        cin >> temp_char;
    }

    // alt_inst.set_dest(45);



    // alt_inst.print_landmarks();


    // alt_inst.print_land_mark_distances();

    
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

        for(int j = 1; j < graph_coords.size() + 1; ++j)
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
                // cout << "avg distance_to_node is: " << average_dist << endl;
                // cout << "cur_furthest_dist is: " << cur_furthest_avg << endl;
                cur_furthest_avg = average_dist;
                cur_furthest_node = j;
            }
        }


        // cout << "Node " << cur_furthest_node << " is being inserted" << endl;
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
    for(int i = 1; i < land_dist.size() + 1; ++i)
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

    for (auto it = landmarks.begin(); it != landmarks.end(); it++)
    {
        int land_mark_node = *it;
        // for(int i = 1; i < graph_coords.size() + 1; ++i)
        // {
        //     double dist_btwn = distance_btw_coords(graph_coords[i], graph_coords[land_mark_node]);
        //     land_dist[i].push_back(std::make_pair(land_mark_node, dist_btwn));
        // }
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


void ALT_Class::alt_alg(int node_count, int src, int dest)
{

    //For timing purposes only
    typedef std::chrono::duration<int,std::milli> millisecs_t;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    g_score = new int[node_count+1];
    f_score = new int[node_count+1];
    path_info = new int[node_count+1];
    for (int i = 0; i < node_count+1; i++){
        g_score[i] = INT_MAX;
        f_score[i] = INT_MAX;
        path_info[i] = -1;
    }

    choose_landmark_index(src, dest);

    std::unordered_set<int> nodes_in_open_set;
    g_score[src] = 0;
    f_score[src] = heuristic_cost_estimate(src, dest);

    std::priority_queue<std::pair<int,int>, 
            std::vector<std::pair<int, int> >, 
            std::function<bool(std::pair<int,int>, std::pair<int,int>)> > open_set(pair_comparator);

    open_set.push(std::make_pair(src, f_score[src]));

    nodes_in_open_set.insert(src);
    
    closed_set.clear();

    while(!open_set.empty())
    {
        int current = open_set.top().first;
        int cur_f_score = open_set.top().second;

        if(closed_set.find(current) != closed_set.end())
        {
            open_set.pop();
            continue;
        }
        
        if(current == dest)
        {
            // Again for timing purposes only
            printSolution(src, dest, g_score[dest], path_info);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end-start));
            std::cout << "That took: " << duration.count() << " milliseconds.\n";
            return;
        }

        open_set.pop();
        nodes_in_open_set.erase(current);
        closed_set.insert(current);

        for (int i = 0; i < graph[current].size(); ++i)
        {
            std::pair<int, int> v = graph[current][i];

            int neighbor = v.first;
            int dist_to_neighb = v.second;
            // Update dist[v.first] only if there is an edge from 
            // u to v, and total weight of path from src to   
            // v through u is smaller than current value of dist[v]

            if(closed_set.find(neighbor) != closed_set.end())
            {
                continue;
            }

            int tentative_g_score = g_score[current] + dist_to_neighb;

            if(tentative_g_score >= g_score[neighbor])
            {
                continue;
            }
            
            path_info[neighbor] = current;
            g_score[neighbor] = tentative_g_score;
            f_score[neighbor] = g_score[neighbor] + heuristic_cost_estimate(neighbor, dest);


            open_set.push(std::make_pair(neighbor, f_score[neighbor]));
            nodes_in_open_set.insert(neighbor);            
        }
    }
}

void ALT_Class::choose_landmark_index(int src, int dest)
{
    double max_heur = INT_MIN;
    double max_index = 0;
    for (int i = 0; i < land_dist[src].size(); ++i)
    {
        double start_dist_land = land_dist[src][i].second;
        double dest_dist_land = land_dist[dest][i].second; 
        double total_dist = std::abs(start_dist_land - dest_dist_land);
        if(total_dist > max_heur)
        {
            max_heur = total_dist;
            max_index = i;
        }

    }

    landmark_index = max_index;

}


double ALT_Class::heuristic_cost_estimate(int start, int dest)
{
    double max_heur = INT_MIN;
    for (int i = 0; i < land_dist[start].size(); ++i)
    {
        double start_dist_land = land_dist[start][i].second;
        double dest_dist_land = land_dist[dest][i].second; 
        double total_dist = std::abs(start_dist_land - dest_dist_land);
        if(total_dist > max_heur)
        {
            max_heur = total_dist;
        }

    }
    
    return max_heur;


    // double start_dist_land = land_dist[start][landmark_index].second;
    // double dest_dist_land = land_dist[dest][landmark_index].second; 
    // double total_dist = std::abs(start_dist_land - dest_dist_land);
    // return total_dist;
}

void ALT_Class::bi_alt_alg(int node_count, int src, int dest)
{

    //For timing purposes only
    typedef std::chrono::duration<int,std::milli> millisecs_t;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    int *g_score_f = new int[node_count+1];
    int *f_score_f = new int[node_count+1];
    int *path_info_f = new int[node_count+1];
    int *g_score_b = new int[node_count+1];
    int *f_score_b = new int[node_count+1];
    int *path_info_b = new int[node_count+1];

    int *final_path_info = new int[node_count+1];
    for (int i = 0; i < node_count+1; i++){
        g_score_f[i] = INT_MAX;
        f_score_f[i] = INT_MAX;
        path_info_f[i] = -1;
        g_score_b[i] = INT_MAX;
        f_score_b[i] = INT_MAX;
        path_info_b[i] = -1;
        final_path_info[i] = -1;
    }

    choose_landmark_index(src, dest);

    
    g_score_f[src] = 0;
    f_score_f[src] = heuristic_cost_estimate(src, dest);

    g_score_b[dest] = 0;
    f_score_b[dest] = heuristic_cost_estimate(dest, src);

    std::priority_queue<std::pair<int,int>, 
            std::vector<std::pair<int, int> >, 
            std::function<bool(std::pair<int,int>, std::pair<int,int>)> > open_set_f(pair_comparator);

    std::priority_queue<std::pair<int,int>, 
            std::vector<std::pair<int, int> >, 
            std::function<bool(std::pair<int,int>, std::pair<int,int>)> > open_set_b(pair_comparator);

    std::unordered_set<int> nodes_in_open_set_f;
    std::unordered_set<int> nodes_in_open_set_b;

    open_set_f.push(std::make_pair(src, f_score_f[src]));
    open_set_b.push(std::make_pair(dest, f_score_b[dest]));

    nodes_in_open_set_f.insert(src);
    nodes_in_open_set_b.insert(dest);

    closed_set_f.clear();
    closed_set_b.clear();
    
    bool found_in_each = false;

    bool forward = false;

    int combined_dist = INT_MAX;
    int crossover_node = -1;


    while(!found_in_each)
    {
        int current = 0;
        int cur_f_score = 0;
        if(open_set_f.top().second <= open_set_b.top().second)
        {
            current = open_set_f.top().first;
            cur_f_score = open_set_f.top().second;
            if(closed_set_b.find(current) != closed_set_b.end())
            {
                combined_dist = g_score_f[current] + g_score_b[current];
                crossover_node = current;
                break;
            }
            open_set_f.pop();
            nodes_in_open_set_f.erase(current);
            closed_set_f.insert(current);
            forward = true;
        }
        else
        {
            current = open_set_b.top().first;
            cur_f_score = open_set_b.top().second;
            if(closed_set_f.find(current) != closed_set_f.end())
            {
                combined_dist = g_score_f[current] + g_score_b[current];
                crossover_node = current;
                break;
            }
            open_set_b.pop();
            nodes_in_open_set_b.erase(current);
            closed_set_b.insert(current);
            forward = false;
        }

        // if(current == dest)
        // {
        //     // Again for timing purposes only
        //     printSolution(src, dest, g_score_f[dest], path_info_f);
        //     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        //     millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end-start));
        //     std::cout << "That took: " << duration.count() << " milliseconds.\n";
        //     return;
        // }

        
        if(forward)
        {
            for (int i = 0; i < graph[current].size(); ++i)
            {
                std::pair<int, int> v = graph[current][i];

                int neighbor = v.first;
                int dist_to_neighb = v.second;
                // Update dist[v.first] only if there is an edge from 
                // u to v, and total weight of path from src to   
                // v through u is smaller than current value of dist[v]

                if(closed_set_f.find(neighbor) != closed_set_f.end())
                {
                    continue;
                }

                int tentative_g_score = g_score_f[current] + dist_to_neighb;

                if(tentative_g_score >= g_score_f[neighbor])
                {
                    continue;
                }
                
                path_info_f[neighbor] = current;
                g_score_f[neighbor] = tentative_g_score;
                f_score_f[neighbor] = g_score_f[neighbor] + heuristic_cost_estimate(neighbor, dest);
                // if(nodes_in_open_set_f.find(neighbor) == nodes_in_open_set_f.end())
                // {
                    open_set_f.push(std::make_pair(neighbor, f_score_f[neighbor]));
                    nodes_in_open_set_f.insert(neighbor);
                // }
            }
        }
        else
        {
            for (int i = 0; i < graph[current].size(); ++i)
            {
                std::pair<int, int> v = graph[current][i];

                int neighbor = v.first;
                int dist_to_neighb = v.second;
                // Update dist[v.first] only if there is an edge from 
                // u to v, and total weight of path from src to   
                // v through u is smaller than current value of dist[v]

                if(closed_set_b.find(neighbor) != closed_set_b.end())
                {
                    continue;
                }

                int tentative_g_score = g_score_b[current] + dist_to_neighb;

                if(tentative_g_score >= g_score_b[neighbor])
                {
                    continue;
                }
                
                path_info_b[neighbor] = current;
                g_score_b[neighbor] = tentative_g_score;
                f_score_b[neighbor] = g_score_b[neighbor] + heuristic_cost_estimate(neighbor, src);
                // if(nodes_in_open_set_b.find(neighbor) == nodes_in_open_set_b.end())
                // {
                    open_set_b.push(std::make_pair(neighbor, f_score_b[neighbor]));
                    nodes_in_open_set_b.insert(neighbor);
                // }
            }
        }
    }

    while (!open_set_f.empty()){
        int node_num = open_set_f.top().first;
        if (g_score_b[node_num] != INT_MAX){
            if (g_score_f[node_num] + g_score_b[node_num] < combined_dist){
                combined_dist = g_score_f[node_num] + g_score_b[node_num];
                crossover_node = node_num;
            }
        }
        open_set_f.pop();
    }
    while (!open_set_b.empty()){
        int node_num = open_set_b.top().first;
        if (g_score_f[node_num] != INT_MAX){
            if (g_score_f[node_num] + g_score_b[node_num] < combined_dist){
                combined_dist = g_score_f[node_num] + g_score_b[node_num];
                crossover_node = node_num;
            }
        }
        open_set_b.pop();
    }
    // Update final_path_info with paths from both ends
    int path_finder = crossover_node;
    // Get from u to src
    while (path_info_f[path_finder] != -1){
        final_path_info[path_finder] = path_info_f[path_finder];
        path_finder = path_info_f[path_finder];
    }
    // Get from dest to u
    path_finder = crossover_node;
    while (path_info_b[path_finder] != -1){
        final_path_info[path_info_b[path_finder]] = path_finder;
        path_finder = path_info_b[path_finder];
    }

    printSolution(src, dest, combined_dist, final_path_info);
    // Again for timing purposes only
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end-start));
    std::cout << "That took: " << duration.count() << " milliseconds.\n";
}
