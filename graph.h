#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <set>
#include <queue>
#include <math.h>


using namespace std;


struct Vertices //the vertices in the graph
{
    int index;
    int xPos;
    int yPos;
    bool visited;
    vector <int> neighboars;
};

struct Edges //the edges as given by the MST
{
    int fromIndex;
    int toIndex;
    double length;
    bool revised;
};


class Graph
{
public:
    Graph(std::vector< pair<int, int> > initialsolution, std::vector<std::pair<int,int> > coordinates);
    //constructor: construct graph from the file
    
    ~Graph();
private:
    
    std::vector<Vertices> allVert;
    std::vector<Edges> allEdge;
    
};

