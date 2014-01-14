#include "mpi.h"
#include "doublylinked.h"
#include "MST.hpp"
#include "GA.hpp"
#include "TSP_LK.h"
#include "HashTable.hpp"



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

vector<double> extractStrategy( vector<double> *incomingMessage);

void mixedStrategy(vector<double>* s1, vector<double>* s2) ;

static void slave(string filename,
                  vector<doublylinkedlist*>* groupGA,
                  std::vector<double> *edgeWeight,
                  std::vector<std::pair<int,int> > *coordinates,
                  std::vector<std::pair<int,int> > *vertexPair);

void retrieveStrategy(
			double *incomingMessage, 
			vector<double>* MethodSequence, 
			vector<double>* MethodIteration);

doublylinkedlist* startingDLL(string filename);