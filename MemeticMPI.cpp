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
#include <stdbool.h> // for boolean

#define MessageTag 1
#define SumTag 2
using namespace std;

int main(int argc, char** argv) {
	/* Initialize MPI environment */
	int ITERATION = 1000;//number of iterations
	int iteration;
	/* Assign size of the matrix */
    int world_size;		//equivalent to p
	int rankWorld;		//same as the id of the processor
	double t1,t2;      //start time and end time
    
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); //get the size of the world
	MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);//get the rank of this processor
	MPI_Status status;
	MPI_Request request;
	

    /*
     * Prepare initial condition of the problem
     */
    //First read in the problem
    vector<double> edgeWeight; //edgeWeight is coupled wth the vertexPair function
	vector<std::pair<int,int> > coordinates; //later expanded in getEdgeWeight function
	vector<std::pair<int,int> > vertexPair; //later expanded in getEdgeWeight function
    string filename = "testDist.txt"; //Give the problem a file name, the file has o be in the same folder as the code

    getEdgeWeight(&edgeWeight, &coordinates, &vertexPair, filename);
    
    int n = coordinates.size();
    int xPos[n], yPos[n],ind[n];
    int count = 0;
    for (vector<int>::iterator it = coordinates.begin(); it != coordinates.end(); it++) {
        ind[count] = count;
        xPos[count] = (*it).first;
        yPos[count] = (*it).second;
        count ++;
    }
    
    //finished initializing the element of doublylinkedlist
    doublylinkedlist* newDLL = new doublylinkedlist();
    newDLL->createList(ind, xPos, yPos, n);
    
    
    /*
     * Generate initial Method selection and number of iterations
     */
    //Initial Sequence matrix
    int[6] Sequence = {0,1,2,3,4,5};
    Sequence = random_shuffle(Sequence.begin(), Sequence.end());
    //debug:
    printf("sequence");
    for (int i =0 ; i <6 ; i++) {
        printf('%d', Sequence[i]);
    }
    //frequency matrix
    int[6] Frequency = {0,0,0,0,0,0}
    Frequency(randWorld) = 1;
    printf("frequency");
    for (int i =0 ; i <6 ; i++) {
        printf('%d', Frequency[i]);
    }
    
    
    
/*
//    
//    
//	//initialize elements
//	int iter;
//	for (iter = 0; iter < numOfElement; iter ++) {
//		elements[iter] = iter;
//		int globalCoords[2];
//		elementIndexToGlobalPosition(globalCoords, rankWorld, iter ,n, numOfCol, world_size);
//		if (globalCoords[0] == globalCoords[1]) {
//			elements[iter] = (double)globalCoords[0]*sin(sqrt((double)globalCoords[0]));
//		}
//		else {
//			elements[iter] = pow((double)globalCoords[0]+(double)globalCoords[1],1.1);
//		}
//	}
//	//first wait till initiation is done
//	MPI_Barrier(MPI_COMM_WORLD ) ;
//	//Print out the time
//	if (rankWorld == 0) {
//		t1 = MPI_Wtime();
//	}
//	
//    
//	//iteration of matrix update for 1000 times to get the sum
//	
//	for (iteration =0; iteration < ITERATION; iteration++) {
//		double bufferIn [n];
//		double bufferOut[n];
//		int iter;
//		
//		//send information to the adjacent ones
//		if (rankWorld != 0) {
//			//first prepare bufferout
//			for (iter = 0; iter < n; iter ++) {
//				bufferOut[iter] = elements[iter*numOfCol];
//			}
//			//send the first column to the left with immediate send
//			MPI_Isend(&bufferOut, n, MPI_DOUBLE, rankWorld-1, MessageTag, MPI_COMM_WORLD, &request);
//		}
//		//using blocking receive
//		if (rankWorld != world_size-1) {
//			MPI_Recv(&bufferIn, n, MPI_DOUBLE, rankWorld+1, MessageTag, MPI_COMM_WORLD, &status);
//		}
//		else {
//			for (iter=0;iter <n;iter ++) {
//				bufferIn[iter] = 0; //just for initiation, this buffer is not used...
//			}
//		}
//		
//		//wait till all information has been sent or received, now update teh matrix
//		if (rankWorld < world_size -1 && rankWorld > 0)
//			MPI_Wait(&request,&status);
//		
//		
//		//updates the matrix element
//		updateMatrixElement(numOfCol, bufferIn, elements, n, rankWorld, world_size);
//	}
//	
//	
//    
//	double sum = 0;
//	int i;
//	for (i = 0; i < numOfElement; i++) {
//		sum += fabs(elements[i]);
//	}
//	
//	double overallSum;
//	MPI_Reduce(&sum, &overallSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    
//	
//	if (rankWorld == 0) {
//		t2 = MPI_Wtime();
//		printf("Auto Detect: number of processors: %d; size of the matrix A = %d*%d;\n",world_size,n,n);
//		printf("Start Time: %lf; \n",t1);
//		printf("End Time: %lf; \nTime Elapse: %lf s = %lf min;\nverification sum: %lf \n",t2, t2-t1, (t2-t1)/60, overallSum);
//	}
//	
//	//free malloc memory
//	free(elements);
//	// Finalize the MPI environment. No more MPI calls can be made after this */
	MPI_Finalize();
	return 0;
    
}

