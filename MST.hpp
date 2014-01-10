#ifndef MST_H
#define MST_H

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <set>
#include <queue>
#include <math.h>
#include <string.h>

#include "graph.h"  //for constructing the graph
#include "doublylinked.h" //for constructing a dll as a solution
#include "DisjointSets.h"  //for MST execution


//#define DEBUG
using namespace std;

void getEdgeWeight(std::vector<double>*, std::vector<std::pair<int,int> >*,std::vector<std::pair<int,int> >*, string);
void quickSort(double arr[], std::vector<double>*, std::vector<std::pair<int,int> >*, int left, int right);
std::vector< std::pair<int,int> > kurskalsAlgo(std::vector<double>*, std::vector<std::pair<int,int> >*, std::vector<std::pair<int,int> >*,DisjointSets*);
doublylinkedlist* DLLFromMST(std::vector<double> edgeWeight,std::vector<std::pair<int,int> > coordinates, std::vector<std::pair<int,int> > vertexPair );

#endif