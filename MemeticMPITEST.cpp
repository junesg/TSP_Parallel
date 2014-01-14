/*
 *  MemeticMPITEST.c
 *
 
 TEST FILE WRITTEN TO TEST THE NON-PARALLEL part
 *  Created by Yuan Shangguan
 *  Copyright 2013 
 *  University of Michigan Ann Arbor. All rights reserved.
 *  method: 0: MST, 1:GA, 2:TwoOpt, 3:RayOpt, 4:KL, 5:StarOpt 
 *	message coding: (double)convergence, (double)time, (double)lengthOfMethod, method array, method iteration array
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iterator>
#include "mpi.h"
#include "doublylinked.h"
#include "MST.hpp"
#include "GA.hpp"
#include "TSP_LK.h"
#include "HashTable.hpp"

#define proc_root 0
#define WORKTAG 1
#define DIETAG 2
#define ITERATION 100 //each round of individual island development, we have this number of iterations
#define MEMETICFREQUENCY  0.4 //best method propagates to these



using namespace std;

#define DEBUG

void singleRoundImprovement(doublylinkedlist* solutionDLL,
                            int methodCode, string filename,vector<doublylinkedlist*>* groupGA,
                            vector<double> *edgeWeight, vector<std::pair<int,int> > *coordinates,
                            vector<std::pair<int,int> > *vertexPair);

static void master();
void quickSortProperties( double *convergence,  double *timeTaken,  double*index, int left, int right);
double conver_time_measure (double* converg, double* time, int pivot) ;
vector<double> extractStrategy( vector<double> *incomingMessage);
void mixedStrategy(vector<double>* s1, vector<double>* s2) ;
void retrieveStrategy(double* incomingMessage, vector<double> *MethodSequence, vector<double> *MethodIteration);
doublylinkedlist* startingDLL(string filename);


//for each round, we use the method as presented in the method Sequence vector
//Each method loops through method iterations
void singleRoundImprovement(doublylinkedlist* solutionDLL, 
			int methodCode, string filename,vector<doublylinkedlist*>* groupGA, 
			vector<double> *edgeWeight, vector<std::pair<int,int> > *coordinates, 
			vector<std::pair<int,int> > *vertexPair) {
    double convergence;
 	doublylinkedlist* newSolution;
 	int NUMITERATIONS = 100; /**This number is changeable**/

	/* switch method to run the method for only once */
    switch (methodCode) {
        case 0: //the MST method
            newSolution= DLLFromMST(*edgeWeight, *coordinates, *vertexPair);
            //MST solution is unique, so we will preserve the new solution if it is better
            break;
        
        case 1: //the GA method
            if(groupGA->size() ==0) { //if the group is empty 
         		groupGA = GA_produceGroup(*coordinates); //initial population is created
            } 
            groupGA =  GA_function(groupGA, 1);//one iterations of breeding only
            newSolution = groupGA->at(0);
            break;
            
        case 2://the 2opt method
            newSolution  = TwoOpt(solutionDLL, NUMITERATIONS);
            break;
            
        case 3: //the ray-opt method
            newSolution  = rayOpt(solutionDLL, NUMITERATIONS);
            break;
        case 4://LK
            newSolution = TSP_LK(solutionDLL, NUMITERATIONS);
            break;
        case 5://Star Opt 
        	//here, first input is the current doubly linked list solution, 
        	//the second input is the number of edge exchanges
        	//the third is the number of iterations we want to do 
            newSolution = starOpt(solutionDLL, 3 ,NUMITERATIONS);
            break;
        default:
        	newSolution = solutionDLL;
            break;
    }
    
    convergence = (double)(solutionDLL->getDistance() - 
    						newSolution->getDistance())/(double)(solutionDLL->getDistance());
    if (convergence > 0){ //if the appointed method produces solution that has better path
        solutionDLL->~doublylinkedlist();
        solutionDLL = newSolution;
    }
    return;
}


/* 
 * The Master program is in charge of allocating the methods to each slave, 
 * It then receives the resulting convergence rate from the slaves. 
 */
static void master() {
  	int sizeWorld; sizeWorld = 6;
    int messageLength = 6;

   	HashMap *historyOfCommands;
    
#ifdef DEBUG
    printf("Enter master\n");
#endif
    
  	/* find out the number of processors in the common world communicator */
  	historyOfCommands = new HashMap(sizeWorld); //the table size in the hashmap is fixed to the size of the world
  	LinkedHashEntry* nextRoundMethods = (LinkedHashEntry*)malloc(sizeof(LinkedHashEntry)*sizeWorld); //this linkedHashEntry stores the methods for next rounds
  	
  	/* Initialize the methods */

    double *incomingMessage;
	
  
    
  /* Loop over getting new work requests until there is no more work
     to be done */
    

  		
    	/* Receive results from a slave */
    		incomingMessage = (double *)malloc(sizeof(double)*messageLength);
    incomingMessage[0]=2.0;//time
    incomingMessage[1]= 0.4;//convergence
    incomingMessage[2]= 3;
    incomingMessage[3]= 2;
    incomingMessage[4]= 50;
     incomingMessage[5]= 50;


    
    
    
            vector<double> message;
			for(int i=0; i< messageLength; i++)
				message.push_back (incomingMessage[i]);
    cout<<"message is:";
    for(int i=0; i< messageLength; i++)
        cout<<message.at(i)<<",";
    cout<<endl;
    
    
    incomingMessage[0]=2.0;//time
    incomingMessage[1]= 0.4;//convergence
    incomingMessage[2]= 1;
    incomingMessage[3]= 4;
    incomingMessage[4]= 50;
    incomingMessage[5]= 50;
    
    vector<double> message2;
    for(int i=0; i< messageLength; i++)
        message2.push_back (incomingMessage[i]);
    cout<<"message2 is:";
    for(int i=0; i< messageLength; i++)
        cout<<message2.at(i)<<",";
    cout<<endl;
    
			historyOfCommands->put(1, &message);
            historyOfCommands->put(2, &message);
			vector<double> strategyContent = extractStrategy(&message);
    printf("strategy content:\n");
    for (int i=0; i<strategyContent.size(); i++) {
        printf("%f, ",strategyContent.at(i));
    }
    printf("\n");
			nextRoundMethods[0].setValue(&strategyContent);
            vector<double> strategyContent2 = extractStrategy(&message2);
    printf("strategy content2:\n");
    for (int i=0; i<strategyContent2.size(); i++) {
        printf("%f, ",strategyContent2.at(i));
    }
    printf("\n");
            nextRoundMethods[1].setValue(&strategyContent2);
    
    printf("strategy content: take two\n");
    for (int i=0; i<strategyContent.size(); i++) {
        printf("%f, ",nextRoundMethods[0].getValue()->at(i));
    }
    printf("\n");
    
    printf("strategy content2: take two\n");
    for (int i=0; i<strategyContent.size(); i++) {
        printf("%f, ",nextRoundMethods[1].getValue()->at(i));
    }
    printf("\n");
    
    
		free(incomingMessage);
		/*sort the converge and timing performance */
    double convergenceAR[5] = {5.0,4.0,2.0,3.0,1.0};
    double timeTaken[5]={1.2,2.0,1.4,1.2,3.0};
    double index[5] = {1,2,3,4,5};
    
    
    printf("convergence   | time| index\n");
    for (int i=0; i<5; i++) {
        printf("%f       |%f   | %f  \n",convergenceAR[i], timeTaken[i], index[i]);
    }
    printf("\n");
    
    
	quickSortProperties( convergenceAR, timeTaken, index, 0, 4 );
    
    
    
		/*store the smallest convergence value == fastest rate of convergence */
		printf("smallest convergence %f",convergenceAR[0]);
		cout<<"after sorting: "<<endl;
    
/*for debugging printing the sorted sequence */
		for	(int i=0; i < sizeWorld-1; i++) {
			cout<<"index: "<< index[i]<<", conv: "<< convergenceAR[i]<<", time: "<< timeTaken[i]<<endl;
		}
/*for debugging */
		
    
        mixedStrategy(nextRoundMethods[0].getValue(), nextRoundMethods[1].getValue());

    printf("strategy content: take 3\n");
    for (int i=0; i<strategyContent.size(); i++) {
        printf("%f, ",nextRoundMethods[0].getValue()->at(i));
    }
    printf("\n");
    
    printf("strategy content2: take 3\n");
    for (int i=0; i<strategyContent.size(); i++) {
        printf("%f, ",nextRoundMethods[1].getValue()->at(i));
    }
    printf("\n");
    
    
		
		for(int sourceI = 1; sourceI < 3; sourceI++) {
			/* Send a new round : (double)NumberInMethod, method array, and send the method iteration array */
    		/* Send the slave a new work unit */
    		vector<double>* test = nextRoundMethods[sourceI-1].getValue();
            cout<<"DEBUG2"<<endl;
            for (int i=0; i<test -> size(); i++) {
                cout<< test->at(i)<<",";
            }
            cout<<endl;
        }
    
    vector<double> MethodSequence, MethodIteration;
    double *messgeABC;
    
    int vectorSize = (int)(nextRoundMethods[0].getValue()->size())+1;
    printf("vecotr size is %d \n",vectorSize);
    

    messgeABC = (double*) malloc(sizeof(double)*vectorSize);
    messgeABC[0] = (double)(vectorSize -1)/(double)2;
    for (int i = 0; i<vectorSize-1; i++) {
        messgeABC[i+1] =nextRoundMethods[0].getValue()->at(i);
    }
    
    
    retrieveStrategy(messgeABC, &MethodSequence, &MethodIteration);
    printf("MethodSequence:\n");
    for (int i=0; i<MethodSequence.size(); i++) {
        printf("%f, ",MethodSequence.at(i));
    }
    printf("\n");
    printf("MethodIterations:\n");
    for (int i=0; i<MethodIteration.size(); i++) {
        printf("%f, ",MethodIteration.at(i));
    }
    printf("\n");

    
}


/* quickSort allows the */
void quickSortProperties( double *convergence,  double *timeTaken,  double *index, int left, int right) {
  int i = left, j = right;
  double tmpConv, tmpTime,tmpInd;

  
  double pivot = conver_time_measure(convergence, timeTaken, ((left + right) / 2));

  /* partition */
  while (i <= j) {
        while (conver_time_measure(convergence, timeTaken, i) <
        	pivot)
              i++;
        while (conver_time_measure (convergence, timeTaken, j)> 
       		pivot)
              j--;
        if (i <= j) {
			tmpInd = index[i];
			tmpConv = convergence[i];
			tmpTime = timeTaken[i];
			index[i]=index[j];
			convergence[i]=convergence[j];
			timeTaken[i] = timeTaken[j];
			index[j] = tmpInd;
			convergence[j] = tmpConv;
			timeTaken[j] = tmpTime;
            i++;
            j--;
        }
  }
/* recursion */
if (left < j)
    quickSortProperties(convergence, timeTaken, index, left, j);
if (i < right)
        quickSortProperties(convergence, timeTaken, index, i, right);
}


/* helper function to define the criteria for sorting */
double conver_time_measure (double* converg, double* time, int pivot) {
	return (converg[pivot])*(time[pivot]);
}

/* Extracts the strategy array and the method iteration array from the diven message */
/* returns a pointer to a vector */
vector<double> extractStrategy( vector<double> *incomingMessage) {
	vector<double> thisVector;
	for (int i = 0; i < incomingMessage->size()-2; i++) {
		thisVector.push_back(incomingMessage->at(i+2));	 
		//only remains in the strategy: (double)numberOfMethods, Method Array, method frequency array
	}
	return thisVector;
}

/* mixes the strategy from s1 and s2 together */
/** preset at 1:1 mixing */
void mixedStrategy(vector<double>* s1, vector<double>* s2) {

	double count1, count2;
	count1 = (double)s1->size()/2; count2 = (double)s2-> size()/2;
	   
	/* First check that if s1 and s2 are no longer divisible*/
	for(int i=count1-1; i < s1->size(); i++) {
		if (s1->at(i) <= 1)
			return;	
	}
	for(int i=count2-1; i < s2->size(); i++) {
		if (s2->at(i) <= 1)
			return;	
	}
	
	vector<double> method1, method2;
	vector<double> iter1, iter2;
	
	for(int i = 0; i< count1; i++) {
		method1.push_back (s1->at(i));
		iter1.push_back ((s1->at(count1+i))/2);
	}
	for(int i = 0; i< count2; i++) {
		method2.push_back (s2->at(i));
		iter2.push_back ((s2->at(count2+i))/2);
	}
	
	s1-> clear();
	s2-> clear();
	
	
	s1-> insert( s1->end(), method1.begin(), method1.end());
	s1-> insert( s1->end(), method2.begin(), method2.end());
	s1-> insert( s1->end(), iter1.begin(), iter1.end());
	s1-> insert( s1->end(), iter2.begin(), iter2.end());

	s2-> insert( s2->end(), method2.begin(), method2.end());
	s2-> insert( s2->end(), method1.begin(), method1.end());
	s2-> insert( s2->end(), iter2.begin(), iter2.end());
	s2-> insert( s2->end(), iter1.begin(), iter1.end());
}




/*
static void slave(string filename,
                  vector<doublylinkedlist*>* groupGA,
                  std::vector<double> *edgeWeight,
                  std::vector<std::pair<int,int> > *coordinates,
                  std::vector<std::pair<int,int> > *vertexPair)
{
  	MPI_Status status;
  	double *incomingMessage;
  	int messageLength;
  	double t1, t2;
  	double convergence;
  	vector<double> MethodSequence, MethodIteration;
  	vector<double> outgoingMessage;
	double oldDist, newDist;
    //first we initialize the filename and get initial doublylinkedlist
    doublylinkedlist* solutionDLL = startingDLL(filename);
    oldDist = solutionDLL->getDistance();
#ifdef DEBUG
    //printf("Starting Solution \n");
    //solutionDLL -> displayforward();
#endif

    

	
#ifdef DEBUG
    printf("right before the loop\n");
#endif

  while (1) {
    // Probe and Receive a message from the master
	MPI_Probe(proc_root, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	// Check the tag of the received message.
    if (status.MPI_TAG == DIETAG) {
      return;
    }

    MPI_Get_count(&status, MPI_DOUBLE, &messageLength);
	incomingMessage = (double *)malloc(sizeof(double)*messageLength);
#ifdef DEBUG
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      printf("%d probe message from root with tag %d and count %d \n", rank, status.MPI_TAG, messageLength);
#endif
      
      
    //receive the message from the root
    MPI_Recv(incomingMessage, 
    		messageLength, 
    		MPI_DOUBLE, 
    		0,  //receive message from root
    		WORKTAG,
            MPI_COMM_WORLD, 
            &status);
#ifdef DEBUG
      for(int i=0; i< messageLength; i++)
          printf("%d,",(int)incomingMessage[i]);
      printf("\n");
#endif
      
    t1 = MPI_Wtime();
    
    retrieveStrategy(incomingMessage, &MethodSequence, &MethodIteration);
#ifdef DEBUG
      for (int i =0; i<MethodSequence.size(); i++) {
          printf("%d,",(int)MethodSequence.at(i));
      }
      printf(" end method sequence\n");
      for (int i =0; i<MethodIteration.size(); i++) {
          printf("%d,",(int)MethodIteration.at(i));
      }
      printf(" end MethodIteration\n");
#endif
      
      
    // start the loop of work
	for(int method=0; method < MethodSequence.size(); method++){
        for (int iter = 0; iter < MethodIteration.at(method); iter++) {
            //do the work
            double oldDist = solutionDLL-> getDistance();
            int methodCode = MethodSequence.at(method);
            singleRoundImprovement(solutionDLL, 
            					methodCode, 
            					filename,
            					groupGA,
            					edgeWeight, 
            					coordinates, vertexPair);    
    	}		
	}
      
	t2 = MPI_Wtime();
	newDist = solutionDLL-> getDistance();
	double convergence = (oldDist - newDist )/oldDist;
	outgoingMessage.push_back(t2-t1);
	outgoingMessage.push_back(convergence);
	outgoingMessage.insert(outgoingMessage.end(), MethodSequence.begin(), MethodSequence.end());
	outgoingMessage.insert(outgoingMessage.end(), MethodIteration.begin(), MethodIteration.end());
	
    // Send the result back
      MPI_Send(&outgoingMessage,
               outgoingMessage.size(),
               MPI_DOUBLE,
               0,
               WORKTAG,
               MPI_COMM_WORLD);
      
      free(incomingMessage);
      outgoingMessage.clear();
      oldDist = newDist;
    
      MethodSequence.clear();
      MethodIteration.clear();
  }
  
}
*/


void retrieveStrategy(double* incomingMessage, vector<double> *MethodSequence, vector<double> *MethodIteration){
	int count = (int) incomingMessage[0];
	for(int i=0; i< count; i++) {
		MethodSequence->push_back(incomingMessage[i+1]);
		MethodIteration->push_back(incomingMessage[i+count+1]);
	}
}

/*
 * Functions below are helper functions for initializing the master
 */
//find the initial doublylinkedlist from the problem file
doublylinkedlist* startingDLL(string filename)
{
    vector<double> edgeWeight; //edgeWeight is coupled wth the vertexPair function
    vector<std::pair<int,int> > coordinates; //later expanded in getEdgeWeight function
    vector<std::pair<int,int> > vertexPair; //later expanded in getEdgeWeight function
    getEdgeWeight(&edgeWeight, &coordinates, &vertexPair, filename);
    
    int n = coordinates.size();
    int xPos[n], yPos[n],ind[n];
    int count = 0;
    for (vector<std::pair<int,int> >::iterator it = coordinates.begin(); it != coordinates.end(); it++) {
        ind[count] = count;
        xPos[count] = (*it).first;
        yPos[count] = (*it).second;
        count ++;
    }
    //finished initializing the element of doublylinkedlist
    doublylinkedlist* newDLL = new doublylinkedlist();
    newDLL->createList(ind, xPos, yPos, n);
    return newDLL;
}


int main(int argc, char** argv)
{
    
    
    
	/* Initializing the solution in all of the processors*/
    string filename = "testDist.txt";
    
    
    //initialize all the vectors for storing information of the problem.
    vector<doublylinkedlist*>* groupGA;
    //initialize parameters to store the problem
    std::vector<double> edgeWeight;
	std::vector<std::pair<int,int> > coordinates;
	std::vector<std::pair<int,int> > vertexPair;
	//function in MST.cpp, puts value into these variables
    getEdgeWeight(&edgeWeight, &coordinates, &vertexPair, filename);
    
    
    
    
    master();
	
	//slave(filename, groupGA, edgeWeight,