#ifndef GA_H
#define GA_H

#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <stdio.h>
#include "TSP_LK.h"
#include "doublylinked.h"

//#include <random>       // std::default_random_engine


#define CROSSK 0.40  //percentage at where we cross over
#define POPULATION 100	//the initial population generated
//criterial for population  <= all combinations (n-1)!
#define breedPop 30 //the size of the breeding population
#define MUTATION 10 //how many links we mutate
//#define LISTSIZE 10 //size of the list --> will be replaced in the future by automatic size detection
#define MAXBREEDITERATION 20



using namespace std;

void GA_produceGroup(vector<pair<int,int> > coordinates, vector<doublylinkedlist*>* groupGA);

void GA_function(vector<doublylinkedlist*>* group, int numberOfIteration);

doublylinkedlist* crossOver1(doublylinkedlist *p1,doublylinkedlist* p2);

doublylinkedlist* GenerateOneSpecies(std::vector<pair<int,int> > coordinates,int seed, int* ind);

void GenerateInitPopulation(std::vector< pair<int,int> > coordinates, int* ind,vector<doublylinkedlist*>* groupGA );


double sortPopDistance(vector< doublylinkedlist* >  *list, vector<float> *distances, int left, int right);
int myrandom (int i);
void PopulationBreeding(std::vector<doublylinkedlist*>* group, double fitDistance );
doublylinkedlist* mutate(doublylinkedlist *) ;

#endif

//end of file
