/*
 *  MemeticMPI.c
 *
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
#include "MemeticMPI.hpp"

using namespace std;



int main(int argc, char** argv)
{
	/* Initialize MPI environment */
    int rankWorld;
    int sizeWorld;
	MPI_Init(&argc, &argv);
    
    //get the rank of this processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
    //get the size of the world, find out how many processors are executing at the same time
    MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

	/* Initializing the solution in all of the processors*/
    string filename = "testDist.txt";
    //first we initialize the filename and get initial doublylinkedlist
    doublylinkedlist* solutionDLL = startingDLL(filename);
	//initialize all the vectors for storing information of the problem.
    vector<doublylinkedlist*>* groupGA;
    //initialize parameters to store the problem
    std::vector<double> *edgeWeight;
	std::vector<std::pair<int,int> > *coordinates;
	std::vector<std::pair<int,int> > *vertexPair;
	//function in MST.cpp, puts value into these variables
  	getEdgeWeight(edgeWeight, coordinates, vertexPair, filename); 

	/* Separate the work division */
	if (rankWorld == proc_root) {
		//receive messages from various processors 
		master();
	}
	else {
    	//message code: double Timetaken, double convergence, double numberOfmethod , double iteration of methods
		slave();
	}
	    
    /* Shut down MPI */
    MPI_Finalize();
    
    return 0;
}


//for each round, we use the method as presented in the method Sequence vector
//Each method loops through method iterations
void singleRoundImprovement(doublylinkedlist* solutionDLL, 
			int methodCode, string filename,vector<doublylinkedlist*>* groupGA, 
			vector<double> *edgeWeight, vector<std::pair<int,int> > *coordinates, 
			vector<std::pair<int,int> > vertexPair) {
    double convergence;
 	doublylinkedlist* newSolution;
 	int NUMITERATIONS = 100; /**This number is changeable**/
 	
	/* switch method to run the method for only once */
    switch (methodCode) {
        case 0: //the MST method
            newSolution= DLLFromMST(edgeWeight, coordinates, vertexPair );
            //MST solution is unique, so we will preserve the new solution if it is better
            break;
        
        case 1: //the GA method
            if(groupGA->isEmpty()) { //if the group is empty 
         		groupGA = GA_produceGroup(coordinates); //initial population is created
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
    
    double convergence = (double)(solutionDLL->getDistance() - 
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
  	int sizeWorld, messageLength, source;
  	MPI_Status status;
  	double t_begin, t_end;
   	double overallConvergence = 1;
   	HashMap *historyOfCommands;
   	
   	
    t_begin = MPI_Wtime();
  	/* find out the number of processors in the common world communicator */
  	MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
  	historyOfCommands = new HashMap(sizeWorld); //the table size in the hashmap is fixed to the size of the world
  	LinkedHashEntry* nextRoundMethods = new LinkedHashEntry[sizeWorld] ; //this linkedHashEntry stores the methods for next rounds

  	
  	/* Initialize the methods */
	int *MethodSequence;
    int *MethodIteration;
    double *incomingMessage;
	
    t_begin = MPI_Wtime();
    
  /* Loop over getting new work requests until there is no more work
     to be done */
  while (overallConvergence == 1) {
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &messageLength);
  		int receiveCount = 0;
  		double convergence[sizeWorld];
  		double timeTaken[sizeWorld];
  		double index[sizeWorld];
  		
  		/* Loop through all results from slaves have been received */
  		while(receiveCount < sizeWorld ){
    	/* Receive results from a slave */
    		incomingMessage = (double *)malloc(sizeof(double)*messageLength);
    		MPI_Recv(incomingMessage,           /* message buffer */
             	messageLength,                 /* one data item */
             	MPI_DOUBLE,        /* of type double real */
             	MPI_ANY_SOURCE,    /* receive from any sender */
             	MPI_ANY_TAG,       /* any type of message */
             	MPI_COMM_WORLD,    /* default communicator */
             	&status);          /* info about the received message */
			source = status.Get_source();
			convergence[source] = (*incomingMessage)[0];
			timeTaken[source] = (*incomingMessage)[1];
			index[source] = source;
			historyOfCommands->put(source, incomingMessage);
			nextRoundMethods->at(source)->setValue(extractStrategy(incomingMessage));
			receiveCount ++;
		}
		incomingMessage->clear();
		
		/*sort the converge and timing performance */
		quickSortProperties( &convergence, &timeTaken, &index, 0, sizeWorld-1 );
		/*store the smallest convergence value == fastest rate of convergence */
		overallConvergence = convergence[0];
		if (onverallConvergence < 1/20) break; //exit while loop if convergence has reached the standard.
		cout<<"after sorting: "<<endl;
/*for debugging printing the sorted sequence */
		for	(int i=0; i < sizeWorld; i++) {
			cout<<"index: "<<(*index)[i]<<", conv: "<<(*convergence)[i]<<", time: "<<(*timeTaken)[i]<<endl;
		}
/*for debugging */
		
		/* Now process the results of this round, prepare for crossing of methods */
		for	(int i = 0; i < (int)sizeWorld/2;  i++) {
			//mixed strategy of the best and the worst, then goes to the center
			mixedStrategy(nextRoundMethods->at(i), nextRoundMethods->at(sizeWorld-1-i));
		}
		
		for(int sourceI = 0; sourceI < sizeWorld; sourceI++) {
			/* Send a new round : (double)NumberInMethod, method array, and send the method iteration array */
    		/* Send the slave a new work unit */
    		MPI_Send(nextRoundMethods->at(sourceI),             /* message buffer */
             		nextRoundMethods->at(sourceI)->at(0),                 /* one data item */
             		MPI_DOUBLE,           /* data item is an integer */
             		sourceI, /* to who we just received from */
             		WORKTAG,           /* user chosen message tag */
             		MPI_COMM_WORLD);   /* default communicator */
  		}
}

  t_end = MPI_Wtime();
  
  	/* Tell all the slaves to exit by sending an empty message with the DIETAG. */
  	for (rank = 1; rank < ntasks; ++rank) {
    	MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  	}
  
	cout<<"Time Spent: "<<t_end-t_begin<<endl;
   
}

/* quickSort allows the */
void quickSortProperties( double *convergence,  double *timeTaken,  double *index, int left, int right) {
  int i = left, j = right;
  double tmpConv, tmpTime,tmpInd;
  
  int pivot = arr[(left + right) / 2];

  /* partition */
  while (i <= j) {
        while (conver_time_measure (convergence, timeTaken, i) < 
        	conver_time_measure (convergence, timeTaken, pivot))
              i++;
        while (conver_time_measure (convergence, timeTaken, j)> 
       		conver_time_measure (convergence, timeTaken, pivot))
              j--;
        if (i <= j) {
			tmpInd = (*index)[i];
			tmpConv = (*convergence)[i];
			tmpTime = (*timeTaken)[i];
			(*index)[i]=(*index)[j];
			(*convergence)[i]=(*convergence)[j];
			(*timeTaken)[i] = (*timeTaken)[j];
			(*index)[j] = tmpInd;
			(*convergence)[j] = tmpConv;
			(*timeTaken)[j] = tmpTime;
              i++;
              j--;
    }
}
/* recursion */
if (left < j)
    quickSort(convergence, timeTaken, left, j);
if (i < right)
        quickSort(convergence, timeTaken, i, right);
}

/* helper function to define the criteria for sorting */
double conver_time_measure (double* converg, double* time, int pivot) {
	return ((*converg)[pivot])*((*time)[pivot]);
}

/* Extracts the strategy array and the method iteration array from the diven message */
/* returns a pointer to a vector */
vector<double>* extractStrategy(incomingMessage) {
	vector<double> thisVector;
	for (int i = 0; i < incomingMessage.size()-2; i++) {
		thisVector.push_back(incomingMessage(i+2));	 
		//only remains in the strategy: (double)numberOfMethods, Method Array, method frequency array
	}
	return &thisVector;
}

/* mixes the strategy from s1 and s2 together */
/** preset at 1:1 mixing */
void mixedStrategy(vector<double>* s1, vector<double>* s2) {

	double count1, count2;
	count1 = s1->at(0); count2 = s2-> at(0);
	
	/* First check that if s1 and s2 are no longer divisible*/
	for(int i=count1+1; i < s1->size(); i++) {
		if (s1->at(i) <= 1)
			return;	
	}
	for(int i=count2+1; i < s2->size(); i++) {
		if (s2->at(i) <= 1)
			return;	
	}
	
	
	vector<double> method1, method2;
	vector<double> iter1, iter2;
	
	for(int i = 0; i< count1; i++) {
		method1.at(i) = s1->at(i+1);
		iter1.at(i) = (s1->at(count1+1+i))/2;
	}
	for(int i = 0; i< count2; i++) {
		method2.at(i) = s2->at(i+1);
		iter2.at(i) = (s2->at(count2+1+i))/2;
	}
	
	s1-> clear();
	s2-> clear();
	
	s1->push_back ( count1+count2);
	s2->push_back (count1+count2);
	
	s1-> insert( s1->end(), method1.begin(), method1.end());
	s1-> insert( s1->end(), method2.begin(), method2.end());
	s1-> insert( s1->end(), iter1.begin(), iter1.end());
	s1-> insert( s1->end(), iter2.begin(), iter2.end());

	s2-> insert( s2->end(), method2.begin(), method2.end());
	s2-> insert( s2->end(), method1.begin(), method1.end());
	s2-> insert( s2->end(), iter2.begin(), iter2.end());
	s2-> insert( s2->end(), iter1.begin(), iter1.end());
}


static void slave() {
  MPI_Status status;
  vector<double> *incomingMessage;
  int messageLength;
  double t1, t2;
  double convergence;
  vector<double> *MethodSequence, *MethodIteration;
  vector<double> outgoingMessage;

  while (1) {
    /* Probe and Receive a message from the master */
	MPI_Probe(proc_root, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	/* Check the tag of the received message. */
    if (status.MPI_TAG == DIETAG) {
      return;
    }
	
    MPI_Get_count(&status, MPI_DOUBLE, &messageLength);
	incomingMessage = (double *)malloc(sizeof(double)*messageLength);
    
    /*receive the message from the root */
    MPI_Recv(incomingMessage, 
    		messageLength, 
    		MPI_DOUBLE, 
    		0,  //receive message from root
    		WORKTAG,
            MPI_COMM_WORLD, 
            &status);

    t1 = MPI_Wtime();
    
    retrieveStrategy(incomingMessage, MethodSequence, MethodIteration);
 
    /* start the loop of work */
	for(int method=0; method < MethodSequence->size(); method++){
        for (int iter = 0; iter < MethodIteration->at(method); iter++) {
            //do the work
            double oldDist = solutionDLL-> getDistance();
            int methodCode = MethodSequence->at(method);
            singleRoundImrovement(solutionDLL, 
            					methodCode, 
            					filename,
            					groupGA,
            					edgeWeight, 
            					coordinates, vertexPair);    
    	}		
	}
	t2 = MPI_Wtime();
	double newDist = solutionDLL-> getDistance();
	double convergence = (oldDist- newDist )/oldDist;
	outgoingMessage.push_back(t2-t1);
	outgoingMessage.push_back(convergence);
	outgoingMessage.push_back(outgoingMessage.end(), MethodSequence.begin(), MethodSequence.end());
	outgoingMessage.push_back(outgoingMessage.end(),MethodIteration.begin(), MethodIteration.end());
	
    /* Send the result back */
    MPI_Send(&outgoingMessage, outgoingMessage.length(), MPI_DOUBLE, 0, WORKTAG, MPI_COMM_WORLD);
    incomingMessage -> clear();
    outgoingMessage -> clear();
  }
  
}



void retrieveStrategy(vector<double>* incomingMessage, vector<double>* MethodSequence, vector<double>* MethodIteration){
	double count = incomingMessage -> at(0);
	vector<double> method; 
	vector<double> iter;
	for(int i=0; i< count; i++) {
		method.push_back(incomingMessage->at(i+1);
		iter.push_back(incomingMessage-> at(i+count+1);
	}
	MethodSequence = &method;
	MethodIteration = &iter;
}

/*
 * Functions below are helper functions for initializing the master
 */
//find the initial doublylinkedlist from the problem file
doublylinkedlist* startingDll(string filename)
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
    return newDll;
}

//end of the file
