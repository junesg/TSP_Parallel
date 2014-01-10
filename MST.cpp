#include "MST.hpp"

doublylinkedlist* DLLFromMST(std::vector<double> edgeWeight,std::vector<std::pair<int,int> > coordinates, std::vector<std::pair<int,int> > vertexPair )
{

	//now sort the edgeWeight array
	double* arr = &edgeWeight[0];
	quickSort(arr, &edgeWeight, &vertexPair, 0, edgeWeight.size()-1);

    DisjointSets set((coordinates).size());
    //create a vector of edges for the Kruskal's algorithm
	std::vector< std::pair<int,int> > mstsolution = kurskalsAlgo(&edgeWeight, &vertexPair, &coordinates, &set);

    
    Graph* solutionGraph = new Graph();
    vector<int>* order = new vector<int>();
    solutionGraph->createGraph(mstsolution,coordinates);
    solutionGraph->DFS(0,order); //starting from the 0th item
    
    int n = order->size();
    int xPos[n], yPos[n],ind[n];
    int count = 0;
    for (vector<int>::iterator it = order->begin(); it != order->end(); it++) {
        ind[count] = *it;
        xPos[count] = coordinates.at(*it).first;
        yPos[count] = coordinates.at(*it).second;
        count ++;
    }
    
    doublylinkedlist* newDLL = new doublylinkedlist();
    newDLL->createList(ind, xPos, yPos, n);
    //Constructed
    /*
    cout<<"SHOW DOUBLYLINKEDLIST: &*********"<<endl;
    newDLL-> displayforward();
    cout<<"END SHOW DOUBLYLINKEDLIST: &*********"<<endl;
    */
    
    delete solutionGraph;
    delete order;
    return newDLL;
}






/*
 * This function takes in
 */
std::vector< std::pair<int,int> > kurskalsAlgo(std::vector<double>* edgeWeight, std::vector<std::pair<int,int> >* vertexPair, std::vector<std::pair<int,int> >* coordinates, DisjointSets* set) {
    
    //create set of empty edges
    std::vector< std::pair<int,int> > MSTedges;
    std::vector< std::pair<int,int> >::iterator MSTIt;

    
	for (int i=0; i< (*vertexPair).size(); i++) {
        int pair1 =(*vertexPair).at(i).first;
        int pair2 =(*vertexPair).at(i).second;

        if ((*set).FindSet(pair1) != (*set).FindSet(pair2)) {

            (*set).Union((*set).FindSet(pair1),(*set).FindSet(pair2));
            MSTedges.push_back(pair<int,int> (pair1, pair2));
        }
    }
    

    return MSTedges;
}




//This function gets the edgeWeight matrix, coupled with a vertexPair matrix (same index in vertexPair gives the edgeWeight.
//It also stores into coordinate matrix the coordinates of each points. (all matrices, and vertex numbers starts from 0)
void getEdgeWeight(std::vector<double> *edgeWeight, std::vector<std::pair<int,int> > *coordinates, std::vector<std::pair<int,int> > *vertexPair, string filename) {
	
    //first get the information of points in the file
	std::string line;
	ifstream file;
	file.open(filename.c_str()); //read file, file name needs to change

    //read each line and save the data into the coordinate vector
	if (!file.eof()){
		while (getline(file,line)){
			int xPos, yPos, index;
			sscanf(line.c_str(),"%d %d %d", &index, &xPos, &yPos);
			(*coordinates).push_back(std::pair<int,int>(xPos,yPos));
		}
		file.close();
	}
	else cout<<"Unable to open file"<<endl;
    
    //size of the vertices set
    int pSize= (*coordinates).size(); 
	
	
    //transfer eldge data into Edge set and corresponding vertice pairs set
    int edgeInd = 0;
    for (int row=0; row<pSize; row++) {
        for (int col=0; col<row; col++) {
            //add the edge into the system
			(*edgeWeight).push_back(sqrt(pow((*coordinates).at(row).first - (*coordinates).at(col).first,2)+ 
					pow((*coordinates).at(row).second - (*coordinates).at(col).second,2)));
            //add the pair of vertex to the relative position in vertexPair
			(*vertexPair).push_back(pair<int,int>(row, col));
        }
    }
}



//quickSort to sort an array 
//the index array originally is (1,2,3,4...) 
//the arr stores the weight of the edges, 
//after sorting, it tells which one is where
void quickSort(double *arr, std::vector<double>* edgeWeight, std::vector<std::pair<int,int> >* vertexPair, int left, int right)
 {
  int i = left, j = right;
  double tmp,tmpEdge;
  std::pair<int,int> tmpPair;
  double pivot = arr[(int)((left + right) / 2)];

  /* partition */
  while (i <= j) {
        while (arr[i] < pivot)
              i++;
        while (arr[j] > pivot)
              j--;
        if (i <= j) {
			//store temporary values
              tmp = arr[i];
			  tmpEdge = (*edgeWeight).at(i);
			  tmpPair = (*vertexPair).at(i);
			//switch to j 
			  arr[i]=arr[j];
			  (*edgeWeight).at(i)=(*edgeWeight).at(j);
			  (*vertexPair).at(i)=(*vertexPair).at(j);
			//using temp to restor j
              arr[j] = tmp;
		      (*edgeWeight).at(j)= tmpEdge;
			  (*vertexPair).at(j)= tmpPair;
              i++;
              j--;
    }
}
/* recursion */
if (left < j)
    quickSort(arr, edgeWeight, vertexPair, left, j);
if (i < right)
   quickSort(arr, edgeWeight,  vertexPair, i, right);
}

