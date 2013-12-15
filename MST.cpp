#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <set>
#include <queue>
#include <math.h>
#include <string.h>


#include "doublylinkedlist.h"
#include "DisjointSets.h"


#define DEBUG
using namespace std;

void getEdgeWeight(std::vector<double>*, std::vector<std::pair<int,int> >*,std::vector<std::pair<int,int> >*, string);
void quickSort(double arr[], std::vector<double>*, std::vector<std::pair<int,int> >*, int left, int right);
void printElementSets(const DisjointSets &s);
std::vector< std::pair<int,int> > kurskalsAlgo(std::vector<double>*, std::vector<std::pair<int,int> >*, std::vector<std::pair<int,int> >*,DisjointSets*);

int main(){

	std::vector<double> edgeWeight; //edgeWeight is coupled wth the vertexPair function
	std::vector<std::pair<int,int> > coordinates; //later expanded in getEdgeWeight function
	std::vector<std::pair<int,int> > vertexPair; //later expanded in getEdgeWeight function
	string filename = "testDist.txt";
	
	getEdgeWeight(&edgeWeight, &coordinates, &vertexPair, filename);

#ifdef DEBUG
	cout<<"Generated coordinates and edges: "<<endl;
	for (int i=0; i<edgeWeight.size();i++){
		cout<<"<"<<vertexPair.at(i).first<<","<<vertexPair.at(i).second<<">  ";
		cout<<"edgeWegith: "<<edgeWeight.at(i)<<endl;
	}
#endif

	//now sort the edgeWeight array
	double* arr = &edgeWeight[0];
	quickSort(arr, &edgeWeight, &vertexPair, 0, edgeWeight.size()-1);


#ifdef DEBUG
	cout<<"After sorting: "<<endl;
	for (int i=0; i<edgeWeight.size();i++){
		cout<<"<"<<vertexPair.at(i).first<<","<<vertexPair.at(i).second<<">   ";
		cout<<"edgeWegith: "<<edgeWeight.at(i)<<endl;
	}
#endif


    DisjointSets set((coordinates).size());
    //create a vector of edges for the Kruskal's algorithm
	std::vector< std::pair<int,int> > mstsolution = kurskalsAlgo(&edgeWeight, &vertexPair, &coordinates, &set);

    
    
    
    
    

   
	return 0;
}

/*
 * This function takes in
 */
std::vector< std::pair<int,int> > kurskalsAlgo(std::vector<double>* edgeWeight, std::vector<std::pair<int,int> >* vertexPair, std::vector<std::pair<int,int> >* coordinates, DisjointSets* set) {
    
    //create set of empty edges
    std::vector< std::pair<int,int> > MSTedges;
    std::vector< std::pair<int,int> >::iterator MSTIt;

    //create disjoint set, and use disjoint set to
    
    
    printElementSets(*set);
    
    
	for (int i=0; i< (*vertexPair).size(); i++) {
        int pair1 =(*vertexPair).at(i).first;
        int pair2 =(*vertexPair).at(i).second;
#ifdef DEBUG
        cout<<"pair 1: "<<pair1<<" , pair 2: "<<pair2<<endl;
#endif
        if ((*set).FindSet(pair1) != (*set).FindSet(pair2)) {
#ifdef DEBUG
            cout<<"Set Ids: "<<(*set).FindSet(pair1)<<" and "<<(*set).FindSet(pair2)<<endl;
            cout<<"Num of sets: "<<(*set).NumSets()<<endl;
#endif
            (*set).Union((*set).FindSet(pair1),(*set).FindSet(pair2));
            MSTedges.push_back(pair<int,int> (pair1, pair2));
        }
    }
   
#ifdef DEBUG
    printElementSets(*set);
    for (std::vector<pair<int,int> >::iterator It = MSTedges.begin();
         It != MSTedges.end(); It++)
    {
        cout<<"("<<It->first<<","<<It->second<<")\t";
        
    }
#endif
    
    
    
    
    return MSTedges;
    
    
    /*
    doublylinkedlist<int> aList;
    int index = 0;
    int x[1];
    x[0]=0;
    //created the first node of the linked list
    aList.createList(x,1);

    //store the previous finished linked list
    int tempInd = 0;
    
    for (index = 1; index < (*coordinates).size(); index++) {
        for (MSTIt = MSTedges.begin(); MSTIt!=MSTedges.end(); MSTIt++) {
            if (MSTIt->first == tempInd) {
                aList.insertAfter(MSTIt->second,tempInd);
                tempInd =MSTIt->second;
            }
            if (MSTIt->second == tempInd) {
                aList.insertAfter(MSTIt->first,tempInd);
                tempInd =MSTIt->second;
            }
        
        }
    }
   
#ifdef DEBUG
    aList.displayforward();
#endif
     */

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
#ifdef DEBUG
            cout<<"vertex pair insert: ("<<row+1<<","<<col+1<<")"<<endl;
#endif
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






#ifdef DEBUG
void printElementSets(const DisjointSets & s)
{
	for (int i = 0; i < s.NumElements(); ++i)
		cout << s.FindSet(i) << "  ";
	cout << endl;
}
#endif
