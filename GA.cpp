#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <stdio.h>
#include "TSP_LK.h"
#include "doublylinked.h"

//#include <random>       // std::default_random_engine


#define CROSSK 0.40  //percentage at where we cross over
#define POPULATION 10	//the initial population generated
	//criterial for population  <= all combinations (n-1)!
#define breedPop 8 //the size of the breeding population
#define MUTATION 3 //how many links we mutate
#define LISTSIZE 10 //size of the list --> will be replaced in the future by automatic size detection



using namespace std;

doublylinkedlist* crossOver1(doublylinkedlist *p1,doublylinkedlist* p2);
doublylinkedlist* GenerateOneSpecies(std::vector<pair<int,int> > coordinates,int seed, int* ind);
vector<doublylinkedlist*>* GenerateInitPopulation(std::vector <pair<int,int> > coordinates, int* ind);
double sortPopDistance(vector< doublylinkedlist* >  *list, vector<float> *distances, int left, int right);
int myrandom (int i);
void PopulationBreeding(std::vector<doublylinkedlist*>* group, double fitDistance );
doublylinkedlist* mutate(doublylinkedlist *) ;

//Overal GA function
int main(){

	int ind[LISTSIZE] = {0,1,2,3,4,5,6,7,8,9};
    int x[LISTSIZE]   = {0,0,1,1,2,2,4,2,1,5};
    int y[LISTSIZE]   = {10,0,10,0,10,0,5,6,7,3};

	int n=LISTSIZE; //test population size
	std::vector<pair<int,int> > coordinates;
	for(int i = 0; i<n; i++) {
		coordinates.push_back(pair<int,int>(x[i],y[i]));
	}
	
	
	cout<<"Coordinates size = "<<coordinates.size()<<endl;	

	vector<doublylinkedlist*>* group;
	group =  GenerateInitPopulation(coordinates,ind);
	
    cout<<"FINISH GENERATION"<<endl;
	
    //DEFINE CONVERGENCE
//    
//    int numberOfIteration = 10;
//    
//    while (numberOfIteration > 0) {
//        numberOfIteration --;
//        //now generate a distance matrix
//        vector<float> distances;
//        for (vector<doublylinkedlist* >::iterator it=group.begin(); it!= group.end(); it++) {
//    //		(*it).displayforward();
//    //		cout<<endl;
//    //		cout<<(*it).getDistance()<<endl;
//            distances.push_back((*it)->getDistance());
//        }
//        
//        
//        
//        
//        cout<<"Distancs has "<<distances.size()<<" elements"<<endl;
//        double fitDistance = sortPopDistance(&group, &distances,0,distances.size()-1);
//        cout<<"fit Distance: "<<fitDistance<<endl;
//
//        for (vector<doublylinkedlist* >::iterator it=group.begin(); it!= group.end(); it++)
//        {
//             (*it)->rearrangeList(0);
//             (*it)->displayforward();
//             cout<<"   with distance: "<<(*it)->getDistance()<<endl;
//             cout<<endl;
//         }
//            
//        
//        //combine for breeding ##############3
//        //NEED TO DO: to select a pair randomly from the top breedProp population.
//        
//        PopulationBreeding(&group, fitDistance );
//
//    }
//    
//    
//    
//    
	//cout<<endl<<" ****** After one round of breeding"<<endl;
	  for (vector<doublylinkedlist *>::iterator it=group->begin(); it!= group->end(); it++) {
			(*it)->rearrangeList(0);	       		       	
			(*it)->displayforward();
			cout<<"   with distance: "<<(*it)->getDistance()<<endl;
	       	cout<<endl;
	 }
//
//   
//
//

/*
	cout<<"***********"<<endl;
	for (vector<float>::iterator it=distances.begin(); it!= distances.end(); it++) {
		cout<<*it<<" ";
	}
	*/
	//Now destroying the group that we have stored.
	//delete[] group;
	cout<<endl<<"***********"<<endl;

}



void PopulationBreeding(std::vector<doublylinkedlist*>* group, double fitDistance ){
	//we substitute the bottom breedPop with the new population
	
	for(int i = 0; i <  POPULATION-breedPop; ) {
		//std::srand ( time(0) );
		//call two pairs randomly from the breeding population
		int random1 = rand()%breedPop;
		//srand(random1);
		int random2 = rand()%breedPop; 
	
		//first the two pairs cross-over
		doublylinkedlist* newlist = crossOver1(group->at(random1) , group->at(random2));

		cout<<"##############New List############## "<<endl; newlist->displayforward(); cout<<endl;
		doublylinkedlist *aList = mutate(newlist);	
		cout<<"After new list###############"<<endl; aList->displayforward(); cout<<endl;
		cout<<"rand1 = "<<random1<<", rand2 = "<<random2<<" and fit distance = "<< fitDistance<<" mutated distance: "<<newlist->getDistance()<<endl;

		//delete newlist;

		if (aList->getDistance() < fitDistance) { //IF THE OFFSPRING IS SUITED (FIT)
			//destroy the current worst
			int size = group->size();
			delete group->at(breedPop + i);
			//add new element to the list
			group->at(breedPop + i)  = copyList(aList,0,aList->countNodes()-1);
			//replace the lower ranked population
			cout<<" is accepted"<<endl;
			i++;
		}
        delete aList;
	}
}



doublylinkedlist* mutate(doublylinkedlist *aList) {
	//first define number of iterations for twoopt
	return TwoOpt(aList,MUTATION);
}



//ordercross over, from Davis 85, Oliver,Smith and Holland (1987) sited in Heinrich Braun(1990) as algorithm1
doublylinkedlist* crossOver1(doublylinkedlist* p1,doublylinkedlist* p2){
	
    int length = p1->countNodes();
	int start = 0;
	int breakPt = (int)length*CROSSK;
	int end = length-1;
	
    //t1 is the part used, t2 is the part thrown away
    doublylinkedlist* t1 = copyList(p1, 0,breakPt);
    doublylinkedlist* t2 = copyList(p2,0, length-1);

    int indices[breakPt+1];
    //int content[3];
    t1->extractIndices(indices);
	cout<<"t1 extracted is "<<endl;
	t1->displayforward();cout<<endl;

    //before dlete
   cout<< "BEFOE DELETE t2 = ";
   t2->displayforward();
    for (int i = 0 ; i <= breakPt; i++) {
        t2->deleteNode(indices[i]);
    }
    cout<<"t2 is now ";
    t2->displayforward();
    cout<<endl;

    cout<<endl;
    t1->appendList(*t2);
    cout<<"ater appending to t1 : "<<endl;
    t1->displayforward();
    cout<<endl;

	delete t2;
    return t1;
}




//Generate a group of linked lists as the initial population of the system
//This takes in the ###WHERE THE PROBLEM IS
vector<doublylinkedlist*>* GenerateInitPopulation(std::vector< pair<int,int> > coordinates, int* ind){
	vector<doublylinkedlist*> *population;
    population = new vector<doublylinkedlist*>(POPULATION);
    
    for (int i=0; i< POPULATION;) {
        population->push_back(GenerateOneSpecies(coordinates,i,ind));
        population->at(i)->displayforward();
        cout<<" pushed back"<<endl;
        bool repeat = false;
        for (int j=0; j<i; j++) {
            //if this list already exists
            if (population->at(j)->compareList(*(population->at(i)))) {
                cout<<" already exisits as the "<<j<<"th element"<<endl;
                delete population->at(i);
                population->erase(population->begin()+j-1);
                repeat = true;
            }
        }
        if (!repeat) {
            i++;
        }
    }
	
    
    /*
    for(int j=0 ; j< POPULATION; ){
		doublylinkedlist* newList;
		newList = GenerateOneSpecies(coordinates,j, ind); //newly created doublylinkedlist
		cout<<"#####ind["<<j<<"] is "<<ind[j]<<";"<<endl;
		newList->rearrangeList(ind[0]); //the index data ==1 node is the start of the list
		//cout<<"generated: "; newList.displayforward();cout<<endl;
		bool exists = false;
		cout<<"&&&&&GENERATED: "; newList->displayforward();cout<<"&&&&&&&&&&&&&"<<endl;

		if(j>0) {
			for(int i=0; i< j; i++){
				cout<<"showing population("<<i<<"):"<<endl;
				population->at(i)->displayforward();
                (*newList).displayforward();
                cout<<"TRYING BEFORE FINISHED:"<<endl;
                bool test = p->compareList(q);
                q.~doublylinkedlist();
                

				exists = (exists || test);

				cout<<"Loop2! exist = "<<exists<<endl;
			}
		}

		if(!exists) {
			population->push_back(newList);
			cout<<"$$$$$$$1$$$$$$$$";
			newList ->displayforward(); cout<<endl;
			population->at(j)->displayforward();cout<<"pushed back1"<<endl;
			cout<<"$$$$$$$2$$$$$$$$  "<<endl;
			//newList->rearrangeList(0);
			//newList->displayforward();
			//cout<<"pushed Backllkk2"<<endl;
			j++;
		}
        cout<<"extra display:";
		population->at(j-1)->displayforward();
	}
    */
    
	return population;
}

//function to create one species using the seed as a random generator seed
doublylinkedlist* GenerateOneSpecies(std::vector< pair<int,int> > coordinates,int seed, int* ind){
	//std::srand ( unsigned ( seed ) );

	int n = coordinates.size(); //number of items
	cout<<" In one speces, the coordinate size = "<<n<<endl;
	vector<int> xPos;
	vector<int> yPos;
	//create the array of indexs
	std::vector<int> arr;

	for(int i = 0; i< n; i++) {
		arr.push_back(ind[i]);
	}
	 random_shuffle (arr.begin(), arr.end());
	//create the doublylinkedlist
	
	int *xArr, *yArr, *index;
	for(int j=0; j<n; j++) {
		xPos.push_back( coordinates.at(arr[j]).first);
		yPos.push_back( coordinates.at(arr[j]).second);
	}
	xArr  = &xPos[0];
	yArr  = &yPos[0];
	index = &arr[0];
	doublylinkedlist* species;
	species = new doublylinkedlist();
	species->createList( index , xArr , yArr, n); 
	species->displayforward(); cout<< "  created"<<endl;
	return species;
}


// random generator function:
int myrandom (int i) { srand(time(NULL));return std::rand()%i;}



//QuickSort the population based on their distance function
//the double it returns is the threshold fitness based on our population  - breeding population
//we keep the breedPop (number of population to be kept) and discard the rest based on their fitness (distance) rank
//This double is the highest distance among the discarded
double sortPopDistance(vector< doublylinkedlist* >  *list, vector<float> *distances, int left, int right){

  	int i = left, j = right;
  	float tmp;
	doublylinkedlist* tmpList;
  	float pivot = (*distances).at((int)(left + right) / 2);

  	/* partition */
  while (i <= j) {
        while ((*distances).at(i) < pivot)
              i++;
        while ((*distances).at(j) > pivot)
              j--;
        if (i <= j) {
              	tmp = distances->at(i);
				tmpList = list->at(i);
              	distances->at(i) = distances->at(j);
				list->at(i) = list->at(j);
				distances->at(j) = tmp;
				list->at(j) = tmpList;
             	i++;
              	j--;
    }
}
/* recursion */
if (left < j)
    sortPopDistance(list, distances, left, j);
if (i < right)
    sortPopDistance(list, distances, i, right);

	
if(breedPop > POPULATION) {
	cout<<"ERROR: please readjust your population and breeding population variables!"<<endl;	
}
	return list->at(breedPop)->getDistance();
}








