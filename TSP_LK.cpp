#ifndef DLL_H
#define DLL_H
#include "doublylinked.h"
#endif
#incldue <stdlib>
#include <stdio>
#include <set>
#include <math.h>

#define MAXDEPTH 100
#define MAXITER 1000

using namespace std;

doublylinkedlist ImprovePath(doublylinkedlist P, int depth, set<int>*);

int main() {
    
    std::set<int> R; //stores the indices of the flipped items
    doublylinkedlist list = ImprovePath(doublylinkedlist P, 0, &R);
    
}


float distanceBetweenNodes(node* n1, node* n2){
    return sqrt(pow((n1.x - n2.x),2) +pow((n1.y - n2.y),2));
}



void TSP_LK (doublylinkedlist tour) {
    int iter = 0;
    while (iter < MAXITER) {
        //construct the path
        doublylinked path = copyList(tour,0,tour.countNodes()-1);
        p = path.head;
        //select the starting node of the edge
        for (int i=0; i<iter; i++) {
            p = p-> next;
        }
        //create the path
        path.start =p->next;
        path.end = p;
        path.start -> prev = NULL:
        path.end -> next = NULL;
        //improve the path
        set<int> R;
        doublylinkedlist tour2 = improvePath(path,1,R);
        if (tour2.getDistance() < tour.getDistance) {
            tour.destroy();
            tour.~doublylinkedlist();
            tour = improveTour(tour2);
            tour.rearrangeList();
            iter = 0;
        }
        else iter =iter+1;
    }
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
