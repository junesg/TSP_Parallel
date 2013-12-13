#include "doublylinkedlist.h"

#define DEBUG

int main(){
    doublylinkedlist<int> aList;
    int x[1];
    x[0]=0;
    
    aList.createList(x,1);
    aList.displayforward();

    
   // aList.insertAfter(4,0);
    aList.insertAfter(5,0);
    aList.displayforward();

    aList.insertAfter(6,0);
    
    aList.displayforward();

    
    
    
    
}


