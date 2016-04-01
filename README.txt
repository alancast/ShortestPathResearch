Repo Rules:
1. ALWAYS PULL BEFORE YOU PUSH. Honestly just use the desktop
    client GUI. We all have Macs and it has a really nice 
    Mac GUI, so use it.
    -Link: https://desktop.github.com/
2. Put datasets in folder at top of repo titled "Datasets".
    -I have already added to the .gitignore "Datasets/*"
    so that none of these potentially large files get pushed
    to git (as you don't want the repo to grow too large).
    -For reference I have downloaded Bay Area, Colorado, 
    Great lakes, Eastern US, and all of USA in all 3 formats
    (distance graph, travel time graph, and coordinates).
    In total this folder is already 4.5GB so DON'T PUSH IT TO GIT
    DON'T PUSH IT TO GIT. DON'T PUSH IT TO GIT
    -link: http://www.dis.uniroma1.it/challenge9/download.shtml
3. Before each function comment what parameters it takes 
    as input and what it returns as output. 
    Example:
    // Reads in graph from input file
    // INPUTS: node_count - initialized to 0 
    //         file_name - string name of file 
    // MODIFIES: node_count - sets it equal to the number of lines in the file
    //                which should be the number of nodes in the graph
    // RETURNS: double pointer to 2D array of graph
    int** readInGraph(int &node_count, const std::string &file_name);
4. Try to keep lines to <60 characters in case we open file
    in VIM or emacs or something
5. Try not to use namespace std (good career habit to get into)

Repo structure:
Directories:
    -Datasets/
        This is a directory that is pretty self explanatory.
        It is comprised of the graphs for the road networks

How to run files:
-djikstras.cpp: ./<exec> <graph_input_file>