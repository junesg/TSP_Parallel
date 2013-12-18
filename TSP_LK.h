#ifndef TSP_H
#define TSP_H

#include <set>
#include "doublylinked.h"


#define MAXDEPTH 2
//#define LISTSIZE 10

using namespace std;

doublylinkedlist ImprovePath(doublylinkedlist P, int depth, vector<int> *R);
float distanceBetweenNodes(node* n1, node* n2);
//methods in the doublylinked.cpp file
doublylinkedlist* copyList(doublylinkedlist P, int start, int end);
doublylinkedlist rayOpt(doublylinkedlist P,int);
doublylinkedlist starOpt(doublylinkedlist P, int K,int);
doublylinkedlist* TSP_LK (doublylinkedlist tour,int) ;
doublylinkedlist TwoOpt(doublylinkedlist,int );

#endif
