#include "mpi.h"
#include "doublylinked.h"
#include "MST.hpp"
#include "GA.hpp"
#include "TSP_LK.h"
#include "HashTable.hpp"
#include <time.h>


void singleRoundImprovement(doublylinkedlist* solutionDLL, 
			int methodCode, string filename,vector<doublylinkedlist*>* groupGA, 
			vector<double> *edgeWeight, vector<std::pair<int,int> > *coordinates, 
			vector<std::pair<int,int> > vertexPair,
            int numberofIterations);
			
static double master() ;

void quickSortProperties( 
			double *convergence,  
			double *timeTaken,  
			double *index, 
			int left, 
			int right);

void BreedingMethod (const int sizeWorld, LinkedHashEntry** nextRoundMethods, const double* index);

double conver_time_measure (double* converg, double* time, int pivot) ;

void extractStrategy( vector<double> *incomingMessage,vector<double> *output) ;

vector<double> mixedStrategy(vector<double>* s1, vector<double>* s2);

void mutateStrategy(vector<double>* strategy);

static doublylinkedlist* slave(
                  string filename,
                  vector<doublylinkedlist*>* groupGA,
                  std::vector<double> *edgeWeight,
                  std::vector<std::pair<int,int> > *coordinates,
                  std::vector<std::pair<int,int> > *vertexPair);

void retrieveStrategy(
			double *incomingMessage, 
			vector<double>* MethodSequence, 
			vector<double>* MethodIteration);

doublylinkedlist* startingDLL(string filename);


//end of file