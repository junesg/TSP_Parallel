#include "mpi.h"
#include "doublylinked.h"
#include "MST.hpp"
#include "GA.hpp"
#include "TSP_LK.h"
#include "HashTable.hpp"

#define proc_root 0
#define WORKTAG 1
#define DIETAG 2
#define ITERATION 100 //each round of individual island development, we have this number of iterations
#define MEMETICFREQUENCY  0.4 //best method propagates to these 



void singleRoundImprovement(doublylinkedlist* solutionDLL, 
			int methodCode, string filename,vector<doublylinkedlist*>* groupGA, 
			vector<double> *edgeWeight, vector<std::pair<int,int> > *coordinates, 
			vector<std::pair<int,int> > vertexPair);
			
static void master() ;

void quickSortProperties( 
			double *convergence,  
			double *timeTaken,  
			double *index, 
			int left, 
			int right) ;

double conver_time_measure (double* converg, double* time, int pivot) ;

vector<double>* extractStrategy( double *incomingMessage) ;

void mixedStrategy(vector<double>* s1, vector<double>* s2) ;

static void slave(string filename);

void retrieveStrategy(vector<double>* incomingMessage, 
			vector<double>* MethodSequence, 
			vector<double>* MethodIteration);

doublylinkedlist* startingDLL(string filename);