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
#include <stdbool.h> // for boolean

#define MessageTag 1
#define SumTag 2


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
	
	
	/* Check for input parameters */
	//test if the world_size is a perfect square
	if (!testSquare(world_size)){
		if (rankWorld ==0 ) {
			printf("Please Try Again With a Perfect Square Number of Processors\n");
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(EXIT_FAILURE);
	}
	//test if the matrix size (n*n) is smaller than the number of processors.
	else if(n < world_size) {
		if (rankWorld ==0) {
			printf("The processor number is too big (>n) and it is not considered for this assignment, aborted\n");
		}
		MPI_Abort(MPI_COMM_WORLD,1);
		exit(EXIT_FAILURE);
	}
	
    
	/*Find the number of elements to be stored in the processor the values in the processor is here */
	//step one, how many values are there in this processor
	numOfCol = ceil((rankWorld+1)*(double)n/(double)world_size)-1-ceil(rankWorld*(double)n/(double)world_size)+1;
	numOfElement = n*numOfCol;
	/*declare initial values*/
	double *elements;
	elements = (double *)malloc(sizeof(double)*(numOfElement));  // allocate n*numOfCol doubles
	
	
	//initialize elements
	int iter;
	for (iter = 0; iter < numOfElement; iter ++) {
		elements[iter] = iter;
		int globalCoords[2];
		elementIndexToGlobalPosition(globalCoords, rankWorld, iter ,n, numOfCol, world_size);
		if (globalCoords[0] == globalCoords[1]) {
			elements[iter] = (double)globalCoords[0]*sin(sqrt((double)globalCoords[0]));
		}
		else {
			elements[iter] = pow((double)globalCoords[0]+(double)globalCoords[1],1.1);
		}
	}
	//first wait till initiation is done
	MPI_Barrier(MPI_COMM_WORLD ) ;
	//Print out the time
	if (rankWorld == 0) {
		t1 = MPI_Wtime();
	}
	
    
	/*iteration of matrix update for 1000 times to get the sum */
	
	for (iteration =0; iteration < ITERATION; iteration++) {
		double bufferIn [n];
		double bufferOut[n];
		int iter;
		
		//send information to the adjacent ones
		if (rankWorld != 0) {
			//first prepare bufferout
			for (iter = 0; iter < n; iter ++) {
				bufferOut[iter] = elements[iter*numOfCol];
			}
			//send the first column to the left with immediate send
			MPI_Isend(&bufferOut, n, MPI_DOUBLE, rankWorld-1, MessageTag, MPI_COMM_WORLD, &request);
		}
		//using blocking receive
		if (rankWorld != world_size-1) {
			MPI_Recv(&bufferIn, n, MPI_DOUBLE, rankWorld+1, MessageTag, MPI_COMM_WORLD, &status);
		}
		else {
			for (iter=0;iter <n;iter ++) {
				bufferIn[iter] = 0; //just for initiation, this buffer is not used...
			}
		}
		
		//wait till all information has been sent or received, now update teh matrix
		if (rankWorld < world_size -1 && rankWorld > 0)
			MPI_Wait(&request,&status);
		
		
		//updates the matrix element
		updateMatrixElement(numOfCol, bufferIn, elements, n, rankWorld, world_size);
	}
	
	
    
	double sum = 0;
	int i;
	for (i = 0; i < numOfElement; i++) {
		sum += fabs(elements[i]);
	}
	
	double overallSum;
	MPI_Reduce(&sum, &overallSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
	
	if (rankWorld == 0) {
		t2 = MPI_Wtime();
		printf("Auto Detect: number of processors: %d; size of the matrix A = %d*%d;\n",world_size,n,n);
		printf("Start Time: %lf; \n",t1);
		printf("End Time: %lf; \nTime Elapse: %lf s = %lf min;\nverification sum: %lf \n",t2, t2-t1, (t2-t1)/60, overallSum);
	}
	
	//free malloc memory
	free(elements);
	// Finalize the MPI environment. No more MPI calls can be made after this
	MPI_Finalize();
	return 0;
    
}

