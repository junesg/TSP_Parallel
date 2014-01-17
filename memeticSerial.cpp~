/*
 *  MemeticSerial.c
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
#include <time.h> //for wall clock time measurement
#include "doublylinked.h"
#include "MST.hpp"
#include "GA.hpp"
#include "TSP_LK.h"
#include "HashTable.hpp"
#include <time.h>
using namespace std;


#define ITERATION 1024 //each round of individual island development, we have this number of iterations
#define DEBUG



void singleRoundImprovement(doublylinkedlist* solutionDLL,
                            int methodCode,
                            string filename,vector<doublylinkedlist*>* groupGA,
                            vector<double> *edgeWeight, vector<std::pair<int,int> > *coordinates,
                            vector<std::pair<int,int> > *vertexPair,
                            int numberofIterations);
void quickSortProperties(
                         double *convergence,
                         double *timeTaken,
                         double *index,
                         int left,
                         int right);

void BreedingMethod (const int sizeWorld, LinkedHashEntry** nextRoundMethods, const double* index);
double conver_time_measure (double* converg, double* time, int pivot) ;
doublylinkedlist* startingDLL(string filename);

void runMethod(
               doublylinkedlist* solutionDLL,
                                   string filename,
                                   vector<doublylinkedlist*>* groupGA,
                                   std::vector<double> *edgeWeight,
                                   std::vector<std::pair<int,int> > *coordinates,
                                   std::vector<std::pair<int,int> > *vertexPair,
                                   double* incomingMessage);
void retrieveStrategy(double* incomingMessage, vector<double> *MethodSequence, vector<double> *MethodIteration);


int main(int argc, char** argv)
{
	/* Initialize serial environment */
    string filename = "pr1002.tsp";
    double timeLimit = 267.311;
    
    
    
    double convergence;
    double timeElapsed = 0;
    //initialize parameters to store the problem
    //initialize all the vectors for storing information of the problem.
    vector<doublylinkedlist*> groupGA;
    std::vector<double> edgeWeight;
    std::vector<std::pair<int,int> > coordinates;
    std::vector<std::pair<int,int> > vertexPair;
   
    getEdgeWeight(&edgeWeight, &coordinates, &vertexPair, filename);

    
    
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&\nSingle Method Loop Through: \n");
    for (int methodCode = 0 ; methodCode < 6; methodCode ++) {
     /*method code changes from 0 to 5 */
        /*variables */
		doublylinkedlist* solutionDLL;
        solutionDLL = startingDLL(filename);
        printf("Starting solution distance is : %f\n", solutionDLL->getDistance());
        
        do {
            time_t Tbegin;
            time(&Tbegin);
            //auto double Tbegin = clock();
            singleRoundImprovement(
                    solutionDLL,
                    methodCode,
                    filename,
                    &groupGA,
                    &edgeWeight,
                    &coordinates,
                    &vertexPair,
                    (int)ITERATION);
            //auto double Tend = clock();
            time_t Tend;
            time(&Tend);
            timeElapsed += double(std::difftime(Tend,Tbegin));
           // printf("time now %f\n",timeElapsed);
        }while (timeElapsed < timeLimit);
        
        printf("\n########################\n");
        solutionDLL->rearrangeList(0);
        solutionDLL->displayforward();
        printf("distance = %f\n", solutionDLL->getDistance());
        printf("########################\n");
        printf("Final Convergence is %f\n",convergence);
        printf("Time taken is %f",timeElapsed);
        printf(" filename is "); cout<<filename<<endl;
        printf("Method is %d \n", methodCode);
        printf("==========================\n\n");
		delete solutionDLL;
    }
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&\n End of Single Method Loop Through: \n");

    
//BELOW IS MULTIPLE LINE
    printf("\n\n\n&&&&&&&&&&&&&&&&&&&&&&&&&\n Multiple method combine Loop Through: \n");
    double incomingMessage[13] = {6,0,1,2,3,4,5, ITERATION/6, ITERATION/6,ITERATION/6,ITERATION/6,ITERATION/6,ITERATION/6};
        /*variables */
    convergence = 1;
    timeElapsed = 0;
    //initialize parameters to store the problem
    //initialize all the vectors for storing information of the problem.
    doublylinkedlist* solutionDLL;
    solutionDLL = startingDLL(filename);
    vector<double> MethodSequence, MethodIteration;
    retrieveStrategy(incomingMessage, &MethodSequence, &MethodIteration);
    printf("Starting solution distance is : %f\n", solutionDLL->getDistance());
    
    do {
        time_t Tbegin;
        time(&Tbegin);
        for(int method=0; method < MethodSequence.size(); method++){
            //do the work
            int methodCode = MethodSequence.at(method);
            //auto double Tbegin = clock();
            runMethod(
                                     solutionDLL,
                                     filename,
                                     &groupGA,
                                     &edgeWeight,
                                     &coordinates,
                                     &vertexPair,
                                     incomingMessage);
            //auto double Tend = clock();
        }
        time_t Tend;
        time(&Tend);
        timeElapsed += double(std::difftime(Tend,Tbegin));
    }while (timeElapsed< timeLimit);
    
    printf("\n########################\n");
    solutionDLL->rearrangeList(0);
    solutionDLL->displayforward();
    printf("distance = %f\n", solutionDLL->getDistance());
    printf("########################\n");
    printf("Time taken is %f",timeElapsed);
    printf(" filename is "); cout<<filename<<endl;
    printf("Method is 0,1,2,3,4,5 rotate \n");
    printf("==========================\n\n");
    
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&\n End of multiple Method Loop Through: \n");


    return 0;
}


//for each round, we use the method as presented in the method Sequence vector
//Each method loops through method iterations
void singleRoundImprovement(
            doublylinkedlist* solutionDLL,
			int methodCode,
            string filename,
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


//Slave works to iteration and send the info back to the master
void runMethod(
               doublylinkedlist* solutionDLL,
                  string filename,
                  vector<doublylinkedlist*>* groupGA,
                  std::vector<double> *edgeWeight,
                  std::vector<std::pair<int,int> > *coordinates,
                  std::vector<std::pair<int,int> > *vertexPair,
                  double* incomingMessage)
{
  	double convergence;
  	vector<double> MethodSequence, MethodIteration;
	double oldDist, newDist;
    //first we initialize the filename and get initial doublylinkedlist
    oldDist = solutionDLL->getDistance();

    retrieveStrategy(incomingMessage, &MethodSequence, &MethodIteration);
    oldDist = solutionDLL-> getDistance();
      
    /* start the loop of work */
	for(int method=0; method < MethodSequence.size(); method++){
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
	convergence = (oldDist - newDist )/oldDist;
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
