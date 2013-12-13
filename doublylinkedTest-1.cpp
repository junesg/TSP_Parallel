#include <iostream>
#include <stdio.h>

#include "doublylinked.h"

#define LISTSIZE 6

using namespace std;

int main(){
    doublylinkedlist<int> aList;
    int ind[LISTSIZE] = {  4,  2,   1,  5,  0,  3};
    //int x[LISTSIZE]   = { 12, 10,  10,  12,  4, 12};
    //int y[LISTSIZE]   = {  5,  1,   1,   5,  4,  5};

    int x[LISTSIZE]   = {0,0,1,1,2,2};
    int y[LISTSIZE]   = {10,0,10,0,10,0};
    
    aList.createList(ind,x,y,LISTSIZE);
	aList.displayforward();

	doublylinkedlist<int> testList = aList.copyList(1,4);
	testList.displayforward();


   
}

