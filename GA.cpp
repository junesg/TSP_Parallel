#include "GA.hpp"


void GA_produceGroup(vector<pair<int,int> > coordinates, vector<doublylinkedlist*>* groupGA){
        
    int n_size = coordinates.size();
    int ind[n_size];
    for (int i=0; i<n_size; i++) {
        ind[i]=i;
    }
    
	GenerateInitPopulation(coordinates,ind, groupGA);

}


vector<doublylinkedlist*>* GA_function(vector<doublylinkedlist*>* group, int numberOfIteration){
		
    //DEFINE CONVERGENCE
    while (numberOfIteration > 0) {
        numberOfIteration --;
        //now generate a distance matrix
        vector<float> distances;
        for (vector<doublylinkedlist* >::iterator it=group->begin(); it!= group->end(); it++) {
            distances.push_back((*it)->getDistance());
        }
        
        //cout<<"Distancs has "<<distances.size()<<" elements"<<endl;
        double fitDistance = sortPopDistance(group, &distances,0,distances.size()-1);
       // cout<<"fit Distance: "<<fitDistance<<endl;

/*
        for (vector<doublylinkedlist* >::iterator it=group->begin(); it!= group->end(); it++)
        {
             (*it)->rearrangeList(0);
             (*it)->displayforward();
             //cout<<"   with distance: "<<(*it)->getDistance()<<endl;
             //cout<<endl;
         }
         */
        PopulationBreeding(group, fitDistance);
    }
    return group;
/*
//	//cout<<endl<<" ****** After one round of breeding"<<endl;
//	  for (vector<doublylinkedlist *>::iterator it=group->begin(); it!= group->end(); it++) {
//			(*it)->rearrangeList(0);	       		       	
//			(*it)->displayforward();
//			cout<<"   with distance: "<<(*it)->getDistance()<<endl;
//	       	cout<<endl;
//           delete *it;
//          // (*it)->~doublylinkedlist();
//	 }
 */
}



void PopulationBreeding(std::vector<doublylinkedlist*>* group, double fitDistance ){
	//we substitute the bottom breedPop with the new population
	
    int breedCount = 0;
	for(int i = 0; i <  POPULATION-breedPop && breedCount < MAXBREEDITERATION; ) {
		//std::srand ( time(0) );
		//call two pairs randomly from the breeding population
		int random1 = rand()%breedPop;
		//srand(random1);
		int random2 = rand()%breedPop; 
	
		//first the two pairs cross-over
		doublylinkedlist* newlist = crossOver1(group->at(random1) , group->at(random2));

		//cout<<"##############New List############## "<<endl; newlist->displayforward(); cout<<endl;
		doublylinkedlist *aList = mutate(newlist);	
		//cout<<"After new list###############"<<endl; aList->displayforward(); cout<<endl;
		//cout<<"rand1 = "<<random1<<", rand2 = "<<random2<<" and fit distance = "<< fitDistance<<" mutated distance: "<<newlist->getDistance()<<endl;
        delete newlist;

		//delete newlist;

		if (aList->getDistance() <= fitDistance) { //IF THE OFFSPRING IS SUITED (FIT)
			//destroy the current worst
			int size = group->size();
//			delete
            group->at(breedPop + i)->~doublylinkedlist();
			//add new element to the list
			group->at(breedPop + i)  = copyList(aList,0,aList->countNodes()-1);
			//replace the lower ranked population
			//cout<<" is accepted"<<endl;
			i++;
		}
        
        
        breedCount ++;
        delete aList;
	}
	
	//after breeding, the group is again arranged according to the population distances
	vector<float> distances;
    for (vector<doublylinkedlist* >::iterator it=group->begin(); it!= group->end(); it++) {
        distances.push_back((*it)->getDistance());
    }
	sortPopDistance(group, &distances,0,distances.size()-1);
}


//memory allocation notes: Mutate will always create a new list
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
	//cout<<"t1 extracted is "<<endl;
	//t1->displayforward();cout<<endl;

    //before dlete
   //cout<< "BEFOE DELETE t2 = ";
   //t2->displayforward();
    for (int i = 0 ; i <= breakPt; i++) {
        t2->deleteNode(indices[i]);
    }
    //cout<<"t2 is now ";
   // t2->displayforward();
   // cout<<endl;
   // t2 -> ~doublylinkedlist();
    
    
    //cout<<endl;
    t1->appendList(*t2);
    //cout<<"ater appending to t1 : "<<endl;
    //t1->displayforward();
//    cout<<endl;

    
    return t1;
}




//Generate a group of linked lists as the initial population of the system
//This takes in the ###WHERE THE PROBLEM IS
void GenerateInitPopulation(std::vector< pair<int,int> > coordinates, int* ind, vector<doublylinkedlist*>* groupGA ){
    
    for (int i=0; i< POPULATION;) {
        groupGA->push_back(GenerateOneSpecies(coordinates,i,ind));
        //population->at(i)->displayforward();
//        cout<<" pushed back"<<endl;
        bool repeat = false;
        for (int j=0; j<i; j++) {
            //if this list already exists
            //cout<<"the "<<j<<"th item is "<<endl;
            //population->at(j)->displayforward();
            
//            cout<<"***** BEGIN COMPARE*******"<<endl;
            if (groupGA->at(j)->compareList(groupGA->at(i))) {
//                cout<<" already exisits as the "<<j<<"th element"<<endl;
                groupGA->at(i)->~doublylinkedlist();
                groupGA->erase(groupGA->begin()+(i-1));
                repeat = true;
            }
//            cout<<"***** AFTER COMPARE*******"<<endl;
        }
        if (!repeat) {
            i++;
        }
    }
	
}

//function to create one species using the seed as a random generator seed
doublylinkedlist* GenerateOneSpecies(std::vector< pair<int,int> > coordinates,int seed, int* ind){
	//std::srand ( unsigned ( seed ) );

	int n = coordinates.size(); //number of items
//	cout<<" In one speces, the coordinate size = "<<n<<endl;
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
//	species->displayforward(); cout<< "  created"<<endl;
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




//end of file



