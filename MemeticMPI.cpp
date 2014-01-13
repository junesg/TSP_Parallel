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

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iterator>
#include "doublylinked.h"
#include "MST.hpp"
#include "GA.hpp"
#include "TSP_LK.h"

#define proc_root 0
#define WORKTAG 1
#define DIETAG 2
#define ITERATION 100 //each round of individual island development, we have this number of iterations
#define MEMETICFREQUENCY 2 //after every ITERATION time, the processor send out to this number of other processors

using namespace std;

double singleRoundImprovement(doublylinkedlist* solutionDLL, 
			int methodCode, string filename,vector<doublylinkedlist*>* groupGA, 
			vector<double> *edgeWeight, vector<std::pair<int,int> > *coordinates, 
			vector<std::pair<int,int> > vertexPair);



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
    doublylinkedlist* solutionDLL = startingDll(string filename);
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
    	//we then send communication to the zeroth the processors
    	//message code: double Timetaken, double convergence, int numberOfmethod , int 
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
   	LinkedHashEntry* nextRoundMethods; //this linkedHashEntry stores the methods for next rounds
   	
   	
    t_begin = MPI_Wtime();
  	/* find out the number of processors in the common world communicator */
  	MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
  	historyOfCommands = new HashMap(sizeWorld); //the table size in the hashmap is fixed to the size of the world
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
    		MPI_Recv(&incomingMessage,           /* message buffer */
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
			receiveCount ++;
		}
		
		quickSortProperties( &convergence, &timeTaken, &index);
		cout<<"after sorting: "<<endl;
/*for debugging */
		for	(int i=0; i < sizeWorld; i++) {
			cout<<"index: "<<(*index)[i]<<", conv: "<<(*convergence)[i]<<", time: "<<(*timeTaken)[i]<<endl;
		}
/*for debugging */
		
		/* Now process the results of this round, prepare for crossing of methods */
		
		
		/* Send a new round : (double)NumberInMethod, method array, and send the method iteration array */
		
		
		



    /* Send the slave a new work unit */
    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
    /* Get the next unit of work to be done */
  }

	t_end = MPI_Wtime();
  /* 
  //There's no more work to be done, so receive all the outstanding
     results from the slaves. 
  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }
*/
  	/* Tell all the slaves to exit by sending an empty message with the DIETAG. */
  	for (rank = 1; rank < ntasks; ++rank) {
    	MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  	}
  
	cout<<"Time Spent: "<<t_end-t_begin<<endl;

  	
 //update method
    //update solution ? --> maybe only MST
    
    //spit out the best solution from the 0 processor, and the distance
    
    //repeat the loop till converge --> convergence criteria

    //spit out the best method and best solution from the 0 processor and the distance
   
}


/* quickSort allows the */
void quickSortProperties( double *convergence,  double *timeTaken,  double *index, int left, int left) {
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



static void 
slave(void)
{
  unit_of_work_t work;
  unit_result_t results;
  MPI_Status status;

  while (1) {

    /* Receive a message from the master */

    MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    /* Check the tag of the received message. */

    if (status.MPI_TAG == DIETAG) {
      return;
    }

    /* Do the work */
        t1 = MPI_Wtime();
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
	double convergence = (newDist- oldDist)/oldDist;

    /* Send the result back */

    MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}



//end of the file
