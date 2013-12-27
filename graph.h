#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>


using namespace std;


typedef struct Vertices //the vertices in the graph
{
    int index;
    int xPos;
    int yPos;
    bool visited;
    vector <int> neighbors;
}Vertices;

typedef struct Edges //the edges as given by the MST
{
    int fromIndex;
    int toIndex;
    //double length;
    bool revised;
}Edges;


class Graph
{
public:
    ~Graph();
    Graph(){
        std::vector<Vertices*> allVert;
        std::vector<Edges*> allEdge;
    }
    
    void createGraph(std::vector< pair<int, int> > initialsolution, std::vector<std::pair<int,int> > coordinates);
    bool isConnected(int,int);
    void printGraphContent();
    std::vector<Vertices*> allVert;
    std::vector<Edges*> allEdge;
    Vertices* findVertex(int ind);
    void DFS(int, vector<int>*);
};


