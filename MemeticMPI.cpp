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
#define ITERATION 1000
using namespace std;

void printArray (int arr[], int n);
void randomize ( int arr[], int n);
void swap (int *a, int *b);



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
	for(int method=0; method<MethodSequence->size(); method++){
        for (int iter = 0; iter < MethodIteration->at(method); iter++) {
            //do the work
            int methodCode = MethodSequence->at(method);
            singleRoundImrovement(solutionDLL, methodCode, filename,groupGA);
            
       
        
    		}		
	}

	
    
    //we then send communication between the processors
    //update method
    //update solution ? --> maybe only MST
    
    //spit out the best solution from the 0 processor, and the distance
    
    //repeat the loop till converge --> convergence criteria

    //spit out the best method and best solution from the 0 processor and the distance
    
    
    /* Shut down MPI */
    MPI_Finalize();
    
    
    
    
    
    
    
    return 0;
}


//for each round, we use the method as presented in the method Sequence vector
//Each method loops through method iterations
void singleRoundImprovement(doublylinkedlist* solutionDLL, int methodCode, string filename,vector<doublylinkedlist*>* groupGA)
{
    
    switch (methodCode) {
        case 0: //the MST method
            <#statements#>
            break;
        
        case 1: //the GA method
            
            break;
        case 2:
            
            break;
        case 3:
            
            break;
        case 4:
            
            break;
        case 5:
            
            break;
        default:
            break;
    }
    
    
    
    
}




//end of the file
