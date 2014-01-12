/*
 *  MemeticMPI.c
 *
 *  Created by Yuan Shangguan
 *  Copyright 2013 
 *  University of Michigan Ann Arbor. All rights reserved.
 *
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iterator>
//#include <stdbool.h> // for boolean
#include "doublylinked.h"
#include "MST.hpp"
#include "GA.hpp"
#include "TSP_LK.h"


#define MessageTag 1
#define SumTag 2
#define ITERATION 10
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
    double t_begin, t_end;
    //get the rank of this processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
    //get the size of the world, find out how many processors are executing at the same time
    MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
    
    vector<int> MethodSequence;
    vector<int> MethodIteration;
    vector<doublylinkedlist*>* groupGA;
    //initialize parameters to store the problem
    std::vector<double> *edgeWeight;
	std::vector<std::pair<int,int> > *coordinates;
	std::vector<std::pair<int,int> > *vertexPair;
  	getEdgeWeight(edgeWeight, coordinates, vertexPair, filename); //function in MST.cpp, puts value into these variables
    
    
    /* Initializing the solution */
    string filename = "testDist.txt";
    //first we initialize the filename and get initial doublylinkedlist
    doublylinkedlist* solutionDLL = startingDll(string filename);
    
    /* method: 0: MST, 1:GA, 2:TwoOpt, 3:RayOpt, 4:KL, 5:StarOpt */
    MethodSequence->push_back(rankWorld%6); //initialize a single method
    MethodIteration->push_back(int ITERATION); //initialize the number of iterations for the initial method
    if (rankWorld == 0) {
        t_begin = MPI_Wtime();
    }

    
    /* start the loop of work */
	for(int method=0; method < MethodSequence->size(); method++){
        for (int iter = 0; iter < MethodIteration->at(method); iter++) {
            //do the work
            int methodCode = MethodSequence->at(method);
            double convergence = singleRoundImrovement(solutionDLL, 
            								methodCode, 
            								filename,
            								groupGA,
            								edgeWeight, 
            								coordinates, vertexPair);        
    	}		
	}
	
	
	
	
    //we then send communication between the processors
    //update method
    //update solution ? --> maybe only MST
    
    //spit out the best solution from the 0 processor, and the distance
    
    //repeat the loop till converge --> convergence criteria

    //spit out the best method and best solution from the 0 processor and the distance
    
    if (rankWorld == 0) {
        t_end = MPI_Wtime();
    }
    cout<<"Time Spent: "<<t_end-t_begin<<endl;
    
    /* Shut down MPI */
    MPI_Finalize();
    
    
    return 0;
}


//for each round, we use the method as presented in the method Sequence vector
//Each method loops through method iterations
double singleRoundImprovement(doublylinkedlist* solutionDLL, 
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
    return convergence;
}




//end of the file
