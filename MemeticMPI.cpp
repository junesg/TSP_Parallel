/*
 *  MemeticMPI.c
 *
 *
 *  Created by Yuan Shangguan
 *  Copyright 2013 University of Michigan Ann Arbor. All rights reserved.
 *
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> // for boolean


//#define DEBUG
#define MessageTag 1
#define SumTag 2


int main(int argc, char** argv) {
	/* Initialize MPI environment */
	int n = 1000;//the side length of the matrix
	int ITERATION = 1000;//number of iterations
	int iteration;
	//Assign size of the matrix
    int world_size;		//equivalent to p
	int rankWorld;		//same as the id of the processor
	int numOfCol;		//number of columns in processor rankWorld
	int numOfElement;	//number of elements in processor rankWorld
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

/*
 * This function takes in the ghost column and computes the updated column elements within the processor
 */
void updateMatrixElement(int numOfCol, double bufferIn[], double element[], int n, int rankWorld, int world_size){
	// declare a pointer variable to point to allocated heap space
	double *elementCopy;
	elementCopy = (double *)malloc(sizeof(double)*(n*numOfCol));  // allocate n*numOfCol doubles
	int i;
	int size;
    
	//update each term in the matrix, and store updated values into elementCopy
	for (i =0; i< n*numOfCol; i++) {
		int globalCoords[2];
		elementIndexToGlobalPosition(globalCoords, rankWorld, i ,n, numOfCol,world_size);
		//if its global index A(i,j), i=j or j =n-1
		if (globalCoords[0] == globalCoords[1] || globalCoords[1] == n-1) {
			elementCopy[i] = element[i];//original value is restored
		}
		
		//if nothing of the above is matching, then we use the matrix updating mechanism
		else {
			int k ;
			double up, down,right;
			//first, if the element does not belong to the first row or the last row or the last column,
			//then this element will be updated with internal elements
			if (globalCoords[0]<n-1) {
				down = fabs(element[i+numOfCol]);
			}
			if (globalCoords[0]>0) {
				up = fabs(element[i-numOfCol]);
			}
			if ((i+1)%numOfCol != 0){
				right = fabs(element[i+1]);
			}
			//if i belongs to the first row
			if (globalCoords[0] == 0) {
				up = fabs(element[i+(n-1)*numOfCol]);
			}
			//if i belongs to the last row
			if (globalCoords[0] == n-1){
				down = fabs(element[i-(n-1)*numOfCol]);
			}
			if ((i+1)%numOfCol == 0 && globalCoords[1]!=n-1) {
				int bufferId = floor ((double) i/(double)numOfCol);
				right = fabs(bufferIn[bufferId]);
			}
            
			elementCopy[i] = 0;
			//do the algorithm of summation with right, up and down values
			for (k = 1; k <11; k++) {
				elementCopy[i] += pow(down,1.0/(double)k) -pow(up,1.0/(double)k)*pow(right,1.0/(double)(k+1));
                
			}
			
			
			if (elementCopy[i] >= 10) {
				elementCopy[i] = 10;
			}
			if (elementCopy[i] <= -10) {
				elementCopy[i] =-10;
			}
			
		}
        
	}
	
	//after update, we store the values back to element
	for (i=0; i< n*numOfCol; i++) {
		element[i] = elementCopy[i];
	}
	free(elementCopy);
	
}




/*
 * This function transforms the elementID of the unit in the array of values stored in the processor to the global
 * matrix position (i,j) cooridnate.
 */
void elementIndexToGlobalPosition(int globalCoords[], int rankWorld, int elementId, int n, int numOfCol, int world_size) {
	
	int Col = ceil(rankWorld*(double)n/(double)world_size);//current very first column of this processor
	int addCol = elementId%numOfCol; //which column the current element is
	int Row = floor((double)elementId/(double)numOfCol);
	Col += addCol;
	
	//Check for system errors
	if (Row > n || Col > n) {
		printf("System Error\n");
		printf("**-------------\n");
		printf("Row is %d, Col is %d\n",Row, Col);
		printf("**-------------\n");
		exit(EXIT_FAILURE);
	}
	
	//assign value to global coordinates
	globalCoords[0] = Row;
	globalCoords[1] = Col;
}


/*
 * This function tests whether a given int is a perfect square.
 */
bool testSquare(int world_size) {
	double d_square = sqrt(world_size);
	return (world_size ==((int)d_square)*((int)d_square));
}
