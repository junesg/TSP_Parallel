#include "graph.h"

//constructor for graph
Graph::Graph(std::vector< pair<int, int> > initialsolution, std::vector<std::pair<int,int> > coordinates1, string filename){
    init_solution = initialsolution;
    coordinates = coordinates1;
    name = filename;
    
    //create vertices of the graph
    std::vector< pair<int, int> > ::iterator it;
    for (it=init_solution.begin(); it!=init_solution.end(); it++) {
        std::set<Vertices>::iterator viter;
        int found =0;
        for (viter = allVert.begin(); viter!= allVert.end(); viter++) {
            if (it->first == viter->vertInd) {
                //this vertex already exists
                (*it).addNeighbor(it->second);
                found = found + 1 ; //found ==1
            }
            if (it->second == viter->vertInd) {
                //this vertex already exists
                (*it).addNeighbor(it->first);
                found = found  + 2; //found ==2 or 3
            }
            if (found == 3) {
                break;
            }
        }
        
        switch (found) {
            case 0://neither is found
                //create vertex 1 and 2, add to allVert
                
                Vertices vec1 = new Vertices();
                break;
                
            case 1://only first one found
                
                break;
            case 2://only second one is found
                
                break;
                
            case 3://both found
                
            default:
                break;
        }

        
    }
}






/*
 * The Vertices Class of the Graph
 *
 */



//constructor of vertices
Vertices::Vertices(int x, int y, int i ){
    xPos = x;
    yPos = y;
    weight= 0;   //default right now.
    vertInd = i;
    visited = false;
    degree = 0;
    //neighboars are null set
}

//destructor of vertices
Vertices::~Vertices(){
    delete[];
}

//function that returns the index
int Vertices::GetIndex(){
    return vertInd;
};

//function that retuns a set of neighbors
std::set<int>  Vertices::GetNeighbors(){
    return neighbors;
}

void Vertices::VisitVert(){
    visited = true;
}

bool Vertices::GetVisit(){
    return visited;
}

void addNeighbor(int n){
    this.neighbors.push_back(n);
}

int getDegree(){
    return degree;
}

void setDegree(int i){
    degree = i;
}


/*
 * The Edge Class of the Graph
 *
 */


//constructor of edge
Edge::Edge(int index, int fromV , int toV, double weight2){
    from = fromV;  //index of vertex that this edge is from
    to = toV;      //index of vertex that this edge points to
    weight= weight2; //the weight of the edge
    edgeInd = index;
    traversed = false;
}

//destructor of edge
Edge::~Edge(){
    delete[];
}

//function that returns the index
int Edge::GetIndex(){
    return edgeInd;
};

//function that retuns the weight
double Edge::GetWeight(){
    return weight;
}

//visit this edge
void Edge::VisitEdge(){
    traversed = true;
}
//function that returns whether this edge has been visited
bool Edge::getTraversed(){
    return traversed;
}

std::pair<int,int> getEnds() {
    return pair<int,int>(fromV, toV);
    
}





