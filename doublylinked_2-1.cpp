#include <iostream>
#include <stdio.h>

#include "doublylinked.h"

#define LISTSIZE 7

using namespace std;

int main(){
//    // rayOpt DEBUG
//    doublylinkedlist<int> aList;
//    int ind[LISTSIZE] = {0,1,2,3,4,5};
//    int x[LISTSIZE]   = {0,0,1,1,2,2};
//    int y[LISTSIZE]   = {10,0,10,0,10,0};
//    aList.createList(ind,x,y,LISTSIZE);
//    aList.rayOpt();
    
    // kOpt DEBUG
    doublylinkedlist<int> bList;
//    int ind[LISTSIZE] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
//    int x[LISTSIZE]   = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
//    int y[LISTSIZE]   = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    int ind[LISTSIZE] = {0,1,2,3,4,5,6};
    int x[LISTSIZE]   = {1,4,7,3,4,4,6};
    int y[LISTSIZE]   = {0,1,2,1,4,20,5};
    bList.createList(ind,x,y,LISTSIZE);
    bList.starOpt(3);
    
}
