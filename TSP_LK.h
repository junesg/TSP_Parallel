#ifndef TSP_H
#define TSP_H

#include <set>
#include "doublylinked.h"


#define MAXDEPTH 2
//#define LISTSIZE 10

using namespace std;

doublylinkedlist* ImprovePath(doublylinkedlist*, int , vector<int>*);
float distanceBetweenNodes(node* n1, node* n2);
//methods in the doublylinked.cpp file
doublylinkedlist* rayOpt(doublylinkedlist* ,int);
doublylinkedlist* starOpt(doublylinkedlist* , int ,int);
doublylinkedlist* TSP_LK (doublylinkedlist* ,int) ;
doublylinkedlist* TwoOpt(doublylinkedlist*,int );

#endif
