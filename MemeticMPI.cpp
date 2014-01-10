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
<<<<<<< HEAD
#include <iostream>
=======
>>>>>>> bba3fe97610c2f33050aad4636e75b8d959918f3
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
    //get the rank of this processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);
    //get the size of the world, find out how many processors are executing at the same time
    MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
    
    vector<int> MethodSequence;
    
    /* Initializing the solution */
    string filename = "testDist.txt";
    //first we initialize the filename and get initial doublylinkedlist
    doublylinkedlist* solutionDLL = startingDll(string filename);
    
    /* method: 0: MST, 1:GA, 2:TwoOpt, 3:RayOpt, 4:KL, 5:StarOpt */
    MethodSequence->push_back(rankWorld%6); //initialize a single method
    
<<<<<<< HEAD
    singleRoundImprovement(solutionDLL, MethodSequence);
    
    
    //do work
    //we then send communication between the processors
    //update method
    //update solution ? --> maybe only MST
    
    //repeat the loop till converge --> convergence criteria
=======
    /*
     * Generate initial Method selection and number of iterations
     */
    //Initial Sequence matrix
    int Sequence[6] = {0,1,2,3,4,5};
    Sequence = random_shuffle(Sequence.begin(), Sequence.end());
    //debug:
    printf("sequence");
    for (int i =0 ; i <6 ; i++) {
        printf("%d", Sequence[i]);
    }
    //frequency matrix
    int Frequency[6] = {0,0,0,0,0,0}
    Frequency(randWorld) = 1;
    printf("frequency");
    for (int i =0 ; i <6 ; i++) {
        printf("%d", Frequency[i]);
    }
>>>>>>> bba3fe97610c2f33050aad4636e75b8d959918f3
    
    
    
    
    /* Shut down MPI */
    MPI_Finalize();
    return 0;
}




//end of the file
