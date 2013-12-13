#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <set>
#include <queue>
#include <math.h>


using namespace std;


class Vertices //the vertices in the graph
{
public:
    Vertices(int xPos, int yPos, int ind);   //constructor
    ~Vertices();            //destructor
    int GetIndex();         //function that returns the index
    std::set<int> GetNeighbors(); //function that retuns a set of neighbors
    void VisitVert(); //function that turns the vertex into visited
    bool GetVisit(); //function that tells whether this vertex has been visited
    void addNeighbor(int n); //add an additional neibor to the set
    
private:
    int vertInd; //the index of the vertex
    int xPos, yPos; //the vertex location
    std::set<int> neighbors; //a set of neighbors that the vertex points to (in case we need them)
    double weight; //the weight of the vertex (in case weight is needed)
    bool visited;
    int degree;
};



class Edge //the vertices in the graph
{
public:
    Edge(int index, int fromV , int toV, double weight);   //constructor
    ~Edge();            //destructor
    int GetIndex();         //function that returns the index
    double GetWeight(); //function that retuns the weight of the edge
    bool getTraversed();
    void VisteEdge();
    std::pair<int, int> getEnds();
    
private:
    int edgeInd; //the index of the edge
    int from, to; //the vertices from which the edge goes to the vertex that the edge ends
    double weight; //the weight of the edge (in our case the distance)
    bool traversed; //whether the edge has been traversed
};


class Graph
{
public:
    //initialsolution is a vector of airs of edges
    //coordinates1 records the (x,y) positions of the ith point at index i
    //filename will be the name of the file from which the graph is imported.
    Graph(std::vector< pair<int, int> > initialsolution, std::vector<std::pair<int,int> > coordinates1, string filename);  //constructor: construct graph from the file
    ~Graph();
private:
    std::set<Vertices> allVert;
    std::set<Edge> allEdge;
    std::vector< pair<int, int> > mst_start;
    std::vector<std::pair<int,int> > coordinates;
    string name;
};

