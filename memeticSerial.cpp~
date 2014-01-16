/*
 *  MemeticMPI.c
 *
 *  Created by Yuan Shangguan
 *  Copyright 2013 
 *  University of Michigan Ann Arbor. All rights reserved.
 *  method: 0: MST, 1:GA, 2:TwoOpt, 3:RayOpt, 4:KL, 5:StarOpt 
 *	message coding: (double)convergence, (double)time, (double)lengthOfMethod, method array, method iteration array
 *  although works with less than 5 processors, recommend working with more than 6 processors to reap the benefits of all functionalities of this library
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iterator>
#include "MemeticMPI.hpp"

using namespace std;


#define proc_root 0
#define WORKTAG 1
#define DIETAG 2
#define EXITTAG 3
#define ITERATION 1024 //each round of individual island development, we have this number of iterations
#define MEMETICFREQUENCY  0.8 //These proportion of communication has been going on
#define DEBUG



//#define DEBUG
int main(int argc, char** argv)
{
	/* Initialize MPI environment */
    int rankWorld;
    int sizeWorld;
	MPI_Init(&argc, &argv);
    double convergence;
    
    //get the rank of this processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
    //get the size of the world, find out how many processors are executing at the same time
    MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

   // chdir("/TSPLIB");
	/* Initializing the solution in all of the processors*/
    string filename = "kroC100.tsp";
    
    //initialize parameters to store the problem
    //initialize all the vectors for storing information of the problem.
    vector<doublylinkedlist*> groupGA;
    std::vector<double> edgeWeight;
	std::vector<std::pair<int,int> > coordinates;
	std::vector<std::pair<int,int> > vertexPair;
    doublylinkedlist* solutionDLL;
    
    

    for(int i=0; i<sizeWorld; i++) {
        if(rankWorld == i) { //to do this is to prevent multiple processors accessing the same file
            getEdgeWeight(&edgeWeight, &coordinates, &vertexPair, filename);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    
	/* Separate the work division */
	if (rankWorld == proc_root) {
		//receive messages from various processors
		convergence = master();
	}
	else {
    	//message code: double Timetaken, double convergence, double numberOfmethod , double iteration of methods
		solutionDLL = slave(filename, &groupGA, &edgeWeight, &coordinates, &vertexPair);
	}
	   
    
    for(int i = 0; i <= sizeWorld; i++ ){
        if (rankWorld == i && i != 0) {
            printf("\n########################\n");
            printf("solution from processor %d is ", rankWorld);
            solutionDLL->rearrangeList(0);
            solutionDLL->displayforward();
            printf("########################\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    for(int i = 0; i <= sizeWorld; i++ ){
        if (rankWorld == i && i != 0) {
            printf("\ndistance from processor %d is ", rankWorld);
            printf("distance = %f\n", solutionDLL->getDistance());
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    printf("Final Convergence is %f\n",convergence);
    if(rankWorld == 0) printf(" filename is "); cout<<filename<<endl;
    
    /* Shut down MPI */
    MPI_Finalize();
    
    return 0;
}


//for each round, we use the method as presented in the method Sequence vector
//Each method loops through method iterations
void singleRoundImprovement(doublylinkedlist* solutionDLL, 
			int methodCode, string filename,
            vector<doublylinkedlist*>* groupGA,
			vector<double> *edgeWeight,
            vector<std::pair<int,int> > *coordinates,
			vector<std::pair<int,int> > *vertexPair,
            int NUMITERATIONS) {
    
    double convergence;
 	doublylinkedlist* newSolution;
 	

	/* switch method to run the method for only once */
    switch (methodCode) {
        case 0: //the MST method
            newSolution= DLLFromMST(*edgeWeight, *coordinates, *vertexPair);
            if (!solutionDLL) {
                if (newSolution->getDistance() > solutionDLL->getDistance()) {
                    newSolution = solutionDLL;
                }
            }
            //MST solution is unique, so we will preserve the new solution if it is better
            break;
        
        case 1: //the GA method
            if(groupGA->empty()) { //if the group is empty
         		GA_produceGroup(*coordinates, groupGA); //initial population is created
            }
            
            GA_function(groupGA, NUMITERATIONS);//one iterations of breeding only
            newSolution =copyList( groupGA->at(0), 0 , groupGA->at(0)->countNodes()-1) ; //chose the top (shortest distance as the new solution
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
        delete solutionDLL;
        solutionDLL = copyList(newSolution, 0, newSolution->countNodes()-1);
        delete newSolution;
    }
    //else solutionDLL does not change

    return;
}


/* 
 * The Master program is in charge of allocating the methods to each slave, 
 * It then receives the resulting convergence rate from the slaves. 
 */
static double master() {
  	int sizeWorld, messageLength;
  	int source;
  	MPI_Status status;
    MPI_Request request;
  	double t_begin, t_end;
   	double overallConvergence = 1;
   	HashMap *historyOfCommands;
    double *incomingMessage;
    
    
  	/* find out the number of processors in the common world communicator */
  	MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
  	historyOfCommands = new HashMap(sizeWorld-1); //the table size in the hashmap is fixed to the size of the world
  	LinkedHashEntry* nextRoundMethods[sizeWorld-1]; //this linkedHashEntry stores the pointer to the methods for next rounds
    for (int i=0; i<sizeWorld-1; i++) {
        nextRoundMethods[i] = new LinkedHashEntry(i);
    }
    

	
    /*first we send out first round of method allocation*/
    for (int rank = 1; rank < sizeWorld; rank ++) {
        double Initialmessage[3] = {1.0, rank%6, (double)ITERATION}; //we only have 6 methods
        MPI_Isend(Initialmessage,             /* message buffer */
                 3,                 /* one data item */
                 MPI_DOUBLE,           /* data item is an integer */
                 rank,
                 WORKTAG,           /* user chosen message tag */
                 MPI_COMM_WORLD,
                 &request);   /* default communicator */
    }

    // Starts to time
    t_begin = MPI_Wtime();
  /* Loop over getting new work requests until there is no more work
     to be done */
  while (overallConvergence > 1/10000) {
		//printf("DEBG1!\n");
  		//int receiveCount = 0;
  		double convergenceAR[(int)sizeWorld-1];
  		double timeTaken[(int)sizeWorld-1];
  		double index[(int)sizeWorld-1];
       // printf("DEBG2!\n");

  		/* Loop through all results from slaves have been received */
        for(int receiveCount =1; receiveCount < sizeWorld; receiveCount ++){
    	/* Receive results from a slave */
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &messageLength);
    		incomingMessage = (double *)malloc(sizeof(double)*messageLength);
    		MPI_Recv(incomingMessage,           /* message buffer */
             	messageLength,                 /* one data item */
             	MPI_DOUBLE,        /* of type double real */
             	MPI_ANY_SOURCE,    /* receive from any sender */
             	MPI_ANY_TAG,       /* any type of message */
             	MPI_COMM_WORLD,    /* default communicator */
             	&status);          /* info about the received message */

			source = (int)status.MPI_SOURCE;

            //MEssage: time, convergence, series
			convergenceAR[(int)source-1] = (double)incomingMessage[1];
			timeTaken[(int)source-1] = (double)incomingMessage[0];
			index[(int)source-1] = source;
			
			vector<double>* message = new vector<double>((int)messageLength);
            
            
            printf("incoming message from %d\n",source);
			for(int i=0; i< messageLength; i++) {
				message->at(i)= incomingMessage[i];
                printf("%f,", message->at(i));
            }
            printf("\n");
            
			historyOfCommands->put(source-1, message);
            vector<double>* strategyContent = new vector<double>((int)messageLength-2);
            extractStrategy(message, strategyContent);
            nextRoundMethods[((int)source-1)]->setValue(strategyContent);

#ifdef DEBUG
            printf("strategy content for %d processor: ", source);
            for (int i = 0 ; i<messageLength-2; i++) {
                printf("%f, ", nextRoundMethods[((int)source-1)]->getValue()->at(i));
            }
            printf("\n");
#endif
		}
      
      		//free(incomingMessage);
		/*sort the converge and timing performance */
		quickSortProperties(convergenceAR, timeTaken, index, 0, (int)sizeWorld-2 );

		/*store the biggest convergence value == fastest rate of convergence */
		overallConvergence = convergenceAR[sizeWorld-2];

		if (overallConvergence < 1/10000) break; //exit while loop if convergence has reached the standard.
		
		/* Now process the results of this round, prepare for crossing of methods */
      
		for	(int i = 0; i < (int)(sizeWorld * MEMETICFREQUENCY) ;  i++) {
			//mixed strategy of the better strategies and the less effective strategies
            //randomly select a strategy from the faster half
            srand(time(NULL)+i);
            //int index1 = rand()%((int)(sizeWorld-1)/2);
            int index1 = rand()%((int)(sizeWorld-1));
            //randomly select a strategy from the slower half
            //int index2 = sizeWorld-2 - rand()%((sizeWorld-1)/2);
            int index2 = rand()%((int)(sizeWorld-1));
			mixedStrategy(nextRoundMethods[index1]->getValue(), nextRoundMethods[index2]->getValue());
		}
      


		for(int sourceI = 1; sourceI < sizeWorld; sourceI++) {
			/* Send a new round : (double)NumberInMethod, method array, and send the method iteration array */
    		/* Send the slave a new work unit */
            
            int messageSize = nextRoundMethods[sourceI-1]->getValue()->size();
            
            vector<double>* sendMessage = nextRoundMethods[sourceI-1]->getValue();
            double outMessage[(messageSize+1)];
            
            outMessage[0] = (double)messageSize/2;
            for (int j= 0; j<messageSize; j++) {
                outMessage[j+1] = sendMessage->at(j);
            }
            
    		MPI_Isend(&outMessage,             /* message buffer */
             		messageSize+1,                 /* one data item */
             		MPI_DOUBLE,           /* data item is an integer */
             		sourceI, /* to who we just received from */
             		WORKTAG,           /* user chosen message tag */
             		MPI_COMM_WORLD,
                      &request);   /* default communicator */
            
            //clean up
            delete nextRoundMethods[sourceI-1]->getValue();
  		}
      
      
}

  t_end = MPI_Wtime();
  
  	/* Tell all the slaves to exit by sending an empty message with the DIETAG. */
  	for (int rank = 1; rank < sizeWorld; ++rank) {
    	MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  	}
    //%%%%%%%%%%%%%%%%%%%%%%print out the history of methods that have been used
  //  historyOfCommands->printMap();
	cout<<"Time Spent: "<<t_end-t_begin<<endl;
    return overallConvergence;
   
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
void extractStrategy( vector<double> *incomingMessage,vector<double> *output) {
	for (int i = 2; i < incomingMessage->size(); i++) {
		output->at(i-2) = (incomingMessage->at(i));
		//only remains in the strategy: (double)numberOfMethods, Method Array, method frequency array
	}
}


/* mixes the strategy from s1 and s2 together */
/** preset at 1:1 mixing */
void mixedStrategy(vector<double>* s1, vector<double>* s2) {
    
	double count1, count2;
	count1 = (double)s1->size()/2;
    count2 = (double)s2-> size()/2;
    
	/* First check that if s1 and s2 are no longer divisible*/
	for(int i=count1; i < s1->size(); i++) {
		if (s1->at(i) <= 1)
			return;
	}
	for(int i=count2; i < s2->size(); i++) {
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


//Slave works to iteration and send the info back to the master
static doublylinkedlist* slave(
                  string filename,
                  vector<doublylinkedlist*>* groupGA,
                  std::vector<double> *edgeWeight,
                  std::vector<std::pair<int,int> > *coordinates,
                  std::vector<std::pair<int,int> > *vertexPair)
{
  	MPI_Status status;
    int rankWorld;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
    
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

  while (1) {
    /* Probe and Receive a message from the master */
	MPI_Probe(proc_root, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	/* Check the tag of the received message. */
    if (status.MPI_TAG == DIETAG) {
       
        return solutionDLL;
    }

    MPI_Get_count(&status, MPI_DOUBLE, &messageLength);
	incomingMessage = (double *)malloc(sizeof(double)*messageLength);
      

      
    /*receive the message from the root */
    MPI_Recv(incomingMessage,   //NOT SURE IF THIS IS RIGH T#($#######################
    		messageLength, 
    		MPI_DOUBLE, 
    		0,  //receive message from root
    		WORKTAG,
            MPI_COMM_WORLD, 
            &status);
      
    t1 = MPI_Wtime();
    
    retrieveStrategy(incomingMessage, &MethodSequence, &MethodIteration);

      
    double oldDist = solutionDLL-> getDistance();
      
    /* start the loop of work */
	for(int method=0; method < MethodSequence.size(); method++){
            //do the work
        int methodCode = MethodSequence.at(method);
        int methodIter =(int)MethodIteration.at(method);
            singleRoundImprovement(solutionDLL,
            					methodCode, 
            					filename,
            					groupGA,
            					edgeWeight, 
            					coordinates,
                                vertexPair,
                                methodIter );
	}
      
    newDist = solutionDLL-> getDistance();

	t2 = MPI_Wtime();
	double convergence = (oldDist - newDist )/oldDist;
    
    outgoingMessage.clear();
	outgoingMessage.push_back(t2-t1);
	outgoingMessage.push_back(convergence);
	outgoingMessage.insert(outgoingMessage.end(), MethodSequence.begin(), MethodSequence.end());
	outgoingMessage.insert(outgoingMessage.end(), MethodIteration.begin(), MethodIteration.end());

      
   
    double *outgoing = &outgoingMessage[0];
      
    /* Send the result back */
      MPI_Send(outgoing,
               outgoingMessage.size(),
               MPI_DOUBLE,
               0,
               WORKTAG,
               MPI_COMM_WORLD);
      
      //free(incomingMessage);
      outgoingMessage.clear();
      oldDist = newDist;
    
      MethodSequence.clear();
      MethodIteration.clear();
  }
  
}



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

//end of file
