#ifndef DLL_H
#define DLL_H
#include "doublylinked.h"
#endif
#incldue <stdlib>
#include <stdio>
#include <set>
#include <math.h>

#define MAXDEPTH 100;
using namespace std;

doublylinkedlist ImprovePath(doublylinkedlist P, int depth, set<int>*);

int main() {
    
    std::set<int> R; //stores the indices of the flipped items
    doublylinkedlist list = ImprovePath(doublylinkedlist P, 0, &R);
    
}


float distanceBetweenNodes(node* n1, node* n2){
    return sqrt(pow((n1.x - n2.x),2) +pow((n1.y - n2.y),2));
}





/*

doublylinkedlist ImprovePath(doublylinkedlist P, int depth, set<int>*){
    
    if (depth < MAXDEPTH) {
        //recursion through all edges
        int firstHead = true;
        node* end = P.head -> prev
        for (node* p = P.head; p != P.head && !firstHead; ) {
            if (distanceBetweenNodes(p, p->next) > distanceBetweenNodes(end,p)) {
                //distance of end to this point is smaller
                if (imrpovement of the whole) {
                    <#statements#>
                }
            }
            
            
            
            
            
            //loop updates
            p = p->next;
            if (p= head)
                firstHead = false;
        }
    }
    else{
        
        
    }
    
    
}
*/
