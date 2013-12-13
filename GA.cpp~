#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <stdio.h>

//#ifndef DLL_H
//#define DLL_H
#include "doublylinked.h"
//#endif

#define CROSSK 0.30

using namespace std;
doublylinkedlist<int> GenerateOneSpecies(std::vector<pair<int,int> > coordinates,int seed);
vector<doublylinkedlist<int> >  GenerateInitPopulation(std::vector <pair<int,int> > coordinates);
void sortPopDistance(vector< doublylinkedlist<int> >  *list, vector<float> *distances, int left, int right);
int myrandom (int i);
#define POPULATION 100;
#define breedPop 30;
#define mutation_perc 10;
#define Iteration 100;

//Overal GA function
int main(){


	int n=4;
	std::vector<pair<int,int> >  coordinates;
	for(int i = 0; i<n; i++) {
		coordinates.push_back(pair<int,int>(i*i,i-1));
	}
	
	vector<doublylinkedlist<int> > group =  GenerateInitPopulation(coordinates);
	//now generate a distance matrix
	vector<float> distances;
	for (vector<doublylinkedlist<int> >::iterator it=group.begin(); it!= group.end(); it++) {
		(*it).displayforward();	
		cout<<endl;	
		(*it).showDistance();
		cout<<endl;
		cout<<(*it).getDistance()<<endl;
		distances.push_back((*it).getDistance());
	}	
	
	sortPopDistance(&group, &distances,0,n);

	//combine for breading ##############3
    




/*
	cout<<"***********"<<endl;
	for (vector<float>::iterator it=distances.begin(); it!= distances.end(); it++) {
		cout<<*it<<" ";
	}	
	for (vector<doublylinkedlist<int> >::iterator it=group.begin(); it!= group.end(); it++) {
		cout<<endl;		
		(*it).displayforward();	
	}
		cout<<endl<<"***********"<<endl;
*/
	
}



doublylinkedlist<int> crossOver1(doublylinkedlist<int> p1,doublylinkedlist<int> p2){
	
    int length = p1.countNodes;
	int start = 0;
	int breakPt = (int)length*CROSSK;
	int end = length-1;
	
    
    doublylinkedlist<int> t1 = p1.copyList(0,breakPt);
    doublylinkedlist<int> t2 = p1.copyList(breakPt+1,end);
    
}


doublylinkedlist<int> crossOver2(doublylinkedlist<int> p1,doublylinkedlist<int> p2){
    
    
    
    
    
}


//Genetic Algorithm

	//establish terminal condition
	
	//randomly genearte a Population of initial solutions
	
	/////recursion (if not terminal)
		//sort the population based on their distances
	//pick breedPop of population and combine to breed 
	//the children are mutated
	//new children are inserted to substitute the bottom breedPop/Population of the whole population
	
	////////recursion




//Generate a group of linked lists as the initial population of the system
//This takes in the 
vector<doublylinkedlist<int> > GenerateInitPopulation(std::vector< pair<int,int> > coordinates){
	vector<doublylinkedlist<int> > population;
	for(int j=0 ; j< 5; j++){
		population.push_back(GenerateOneSpecies(coordinates,j));
	}
	return population;
}

//function to create one species using the seed as a random generator seed
doublylinkedlist<int> GenerateOneSpecies(std::vector< pair<int,int> > coordinates,int seed){
	std::srand ( unsigned ( seed ) );
	int n = coordinates.size(); //number of items
	vector<int> xPos;
	vector<int> yPos;
	//create the array of indexs
	std::vector<int> arr;
	for(int i = 0; i<n; i++) {
		arr.push_back(i);
	}
	//shuffle the array using built-in random generator:
  	std::random_shuffle (arr.begin(), arr.end(), myrandom);

/*
	// print out content:
  	std::cout << "arr contains:";
  	for (std::vector<int>::iterator it=arr.begin(); it!=arr.end(); ++it)
    	std::cout << ' ' << *it;
	std::cout << '\n';
*/

	//create the doublylinkedlist
	
	int *xArr, *yArr, *index;
	for(int j=0; j<n; j++) {
		xPos.push_back( coordinates.at(arr[j]).first);
		yPos.push_back( coordinates.at(arr[j]).second);
	}
	xArr  = &xPos[0];
	yArr  = &yPos[0];
	index = &arr[0];
	doublylinkedlist<int> species;
	species.createList( index , xArr , yArr, n); 
	//species.displayforward();
	return species;
}


// random generator function:
int myrandom (int i) { return std::rand()%i;}



//QuickSort the population based on their distance function
void sortPopDistance(vector< doublylinkedlist<int> >  *list, vector<float> *distances, int left, int right){

  	int i = left, j = right;
  	float tmp;
	doublylinkedlist<int> tmpList;
  	float pivot = (*distances).at((int)(left + right) / 2);

  	/* partition */
  while (i <= j) {
        while ((*distances).at(i) < pivot)
              i++;
        while ((*distances).at(j) > pivot)
              j--;
        if (i <= j) {
              	tmp = (*distances).at(i);
				tmpList = (*list).at(i);
              	(*distances).at(i) = (*distances).at(j);
				(*list).at(i) = (*list).at(j);
              	(*distances).at(j) = tmp;
				(*list).at(j) = tmpList;
             	i++;
              	j--;
    }
}
/* recursion */
if (left < j)
    sortPopDistance(list, distances, left, j);
if (i < right)
    sortPopDistance(list, distances, i, right);
}








