#include "graph.h"
#include <iostream>
#include <ctime>
using namespace std;


//After constructor for the graph, we now have create graph, which puts in necessary parts
void Graph::createGraph(vector< pair<int, int> > initialsolution, vector<pair<int,int> > coordinates){
    
    //construct edges vector
    for (vector<pair<int,int> >::iterator it = initialsolution.begin(); it!= initialsolution.end(); it++) {
        Edges* newEdge = new Edges;
        newEdge->fromIndex = it->first;
        newEdge->toIndex = it->second;
        newEdge->revised = false;
        allEdge.push_back(newEdge);
    }
    //construct vertex vector
    int count = 0;
    for (vector<pair<int,int> >::iterator it = coordinates.begin(); it!= coordinates.end(); it++) {
        Vertices* newVert = new Vertices;
        newVert->index = count;
        newVert->xPos = it->first;
        newVert->yPos = it->second;
        newVert->visited = false;
        allVert.push_back(newVert);
        count ++;
    }
    //Now find the neighbors of the current edges
    for (vector<pair<int,int> >::iterator it = initialsolution.begin();it!= initialsolution.end(); it++) {
        Vertices *p1,*p2;
        p1 = findVertex(it->first);
        p2 = findVertex(it->second);
        p1 -> neighbors.push_back(p2->index);
        p2 -> neighbors.push_back(p1->index);
    }
}

//Destructor for the graph
Graph::~Graph(){
    //delete all vertices
    for (vector<Vertices*>::iterator it = allVert.begin(); it!=allVert.end(); it++) {
        delete (*it);
    }
    //delete all edges
    for (vector<Edges*>::iterator it=allEdge.begin(); it!=allEdge.end(); it++) {
        delete (*it);
    }
}

//given an index to the vertex, find the pointer to the vertex
Vertices* Graph:: findVertex(int ind) {
    for (vector<Vertices*>::iterator it = allVert.begin(); it!= allVert.end(); it++) {
        if ((*it)->index == ind) {
            return (*it);
        }
    }
    printf("Vertex not found.\n");
    return NULL;
}


//Recursive depth first search to get the graph
//algorithm is ismplified as we know that the graph is MST
void Graph::DFS(int x, vector<int>* order) {
    Vertices* v;
    v = findVertex(x);
    v -> visited = true;
    order->push_back(x);
    
    for (vector<int>::iterator it = v->neighbors.begin();
         it != v->neighbors.end(); it++) {
        Vertices* w = findVertex(*it);
        if (!w->visited) {
            DFS(*it,order);
        }
    }
}





int main(){
    vector< pair<int, int> > initialsolution;
    vector<pair<int,int> > coordinates;
    
    initialsolution.push_back(pair<int,int>(0,2));
    initialsolution.push_back(pair<int,int>(1,2));
    initialsolution.push_back(pair<int,int>(3,2));
    initialsolution.push_back(pair<int,int>(3,4));
    initialsolution.push_back(pair<int,int>(4,5));
    
    coordinates.push_back(pair<int,int>(0,0)); //the #0 item
    coordinates.push_back(pair<int,int>(1,0)); //the #1 item
    coordinates.push_back(pair<int,int>(2,0)); //the #2 item
    coordinates.push_back(pair<int,int>(4,3)); //the #3 item
    coordinates.push_back(pair<int,int>(5,3)); //the #4 item
    coordinates.push_back(pair<int,int>(10,10)); //the #5 item

    Graph* thisGraph = new Graph();
    thisGraph->createGraph(initialsolution,coordinates);
    thisGraph->printGraphContent();

    
    vector<int>* order;
    order = new vector<int>();
    thisGraph->DFS(0,order);
    
    for (vector<int>::iterator it=order->begin(); it!=order->end(); it++)
    {
        printf("%d\n",*it);
    }
    
    delete thisGraph;
    delete order;
    
    
}




//This function prints out the graph content
void Graph::printGraphContent(){
    for (vector<Vertices*>::iterator it = allVert.begin(); it!= allVert.end(); it++) {
        printf("Vertex %d, at (%d,%d)\n",(*it)->index, (*it)->xPos, (*it)->yPos);
    }
    for (vector<Edges*>::iterator it = allEdge.begin(); it!= allEdge.end(); it++) {
        printf("Edge (%d,%d)\n",(*it)->fromIndex, (*it)->toIndex);
    }
    
    printf("END\n");
}


