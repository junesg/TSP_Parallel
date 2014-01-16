/*
 * All methods in this file creates a new dll and return the pointer that of the new dll
 */

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "TSP_LK.h"

using namespace std;



/*
int main() {
    doublylinkedlist* aList;
    aList = new doublylinkedlist();
    int ind[LISTSIZE] = {0, 1, 2, 3, 4, 5,6,7,8,9};
    int x[LISTSIZE]   = {0, 0, 1, 1, 2, 2,4,5,6,7};
    int y[LISTSIZE]   = {10,0, 10,0, 10, 0,2,3,4,5};
    aList->createList(ind,x,y,LISTSIZE);
    aList->displayforward();
    cout<<endl;
    cout<<"original distance is "<<aList->getDistance()<<endl;
    
    doublylinkedlist *bList;
//    bList = new doublylinkedlist();
    bList = TSP_LK(aList,5);
    
    cout<<endl<<"final reslult"<<endl;
    bList->displayforward();
    printf("debug");

    cout<<"distance is "<<bList->getDistance()<<endl;
    delete aList;
    delete bList;
   // aList->doublylinkedlist::~doublylinkedlist();
   // bList->doublylinkedlist::~doublylinkedlist();
    
}
*/


float distanceBetweenNodes(node* n1, node* n2){
    return sqrt(pow((n1->x - n2->x),2) +pow((n1->y - n2->y),2));
}


/*
 *TSP_LK takes in a pointer to the class of doublylinkedlist, which is the tour,
 *it returns a pointer to a new dynamically allocated doublylinkedlist
 *MAXITER is the number of iterations that the program iterates through
 */
doublylinkedlist* TSP_LK(doublylinkedlist* thisTour, int MAXITER) {
    int iter = 0;
    int countNode = thisTour->countNodes();
    doublylinkedlist *tour;
    tour = copyList(thisTour, 0 , countNode-1);//created tour##1
    
    while (iter < MAXITER) {
        node* p = tour->head;

        for (int i=0; i<iter; i++ ){
            p=p->next;
        }
        //rearrange to form a path
        tour->rearrangeList(p->data);

        //cout<<"REARRANGE"<<endl;
        //now call improve path on this path
        vector<int> R;
        for (int i=0; i<countNode; i++) {
            R.push_back(0);
        }
        
        //construct the improved path ##
        doublylinkedlist *tour3 = ImprovePath(tour, 1, &R); //construct path ##2
        
        if (tour3->getDistance() < tour->getDistance()) {
        	delete tour;
           	iter = 0;
           	tour = tour3;
        }
       else {
           iter =iter+1;
           delete tour3;
        }
    }
    
    return tour; //return path##1
}

/*
 *Improvepath takes in a pointer to a doublylinkedlist that is put into the improve path
 *It outputs a pointer to the doublylinkedlist that is dynamically allocated
 *depth is the number of iterations that a single path has been improved.
 *R is a pointer to a vector that records which edges have been used R.at(edge-data)=1 or not used R.at(edge->data)=0
 */
doublylinkedlist* ImprovePath(doublylinkedlist* Thispath, int depth, vector<int> *R){
   
    doublylinkedlist* path; //= new doublylinkedlist();
    path = copyList(Thispath,0, Thispath->countNodes()-1);  //construct tour ##1
    
    //cout<<"Improving tour "; path->displayforward();cout<<endl;

    //if there is three nodes in the path, no need to improve
    if (path->countNodes() <=3) {
        return path; //return tour##1
    }
    
    //if the depth is smaller than maxdepth, keeps on improving till a better path is found
    if (depth < MAXDEPTH) {
        for (node* p = path->head->next; p->next!= path->end; p = p->next) {
            if (R->at(p->data)==0) {
                if (distanceBetweenNodes(p,p->next) > distanceBetweenNodes(p,path->end)) {
                    //if tour length is improved
                    //cout<<"Tour length is improved"<<endl;
                    if (distanceBetweenNodes(p,p->next)+
                        distanceBetweenNodes(path->head,path->end) >
                        distanceBetweenNodes(p, path->end) +
                        distanceBetweenNodes(path->head, p->next)) {
                        //path.displayforward(); cout<<endl;
                        doublylinkedlist* tour;
                       // tour = new doublylinkedlist();
                        tour = copyList(path, 0, path->countNodes()-1); //copy a new tour from path ##2
                        tour->flipTwoItems(p->data, path->end->data);  //flip the two items in tour
                        tour->end = tour->head->prev;

                        delete path;
                        path = tour;
                        return path; //return path ##2
                    }
                    else {
                        path->flipTwoItems(p->data,path->end->data);  //just flip the edge
                        int thisData = p->data;
                        R->at(thisData)=1;
                        doublylinkedlist* result = ImprovePath(path, depth+1, R); //construct result ##3;
                        delete path;
                        //path->~doublylinkedlist();//destroy path ##1
                        return result; //return ##3
                    }
                }
            }
        }
    }
    else{
        float maxDist=0;
        node* maxNode;
        node* p;
        if (path->countNodes() <=3) {
            //cout<<"No Need improve path"<<endl;
            return path; //return ##1
        }
        path->end = path->head->prev;
        
        //get the node that will end up giving the biggest gain over the end to head
        for (p = path->head->next; p->next!=path->end; p = p->next) {
             if (distanceBetweenNodes(p,p->next) - distanceBetweenNodes(p,path->end) > maxDist) {
                 maxDist =distanceBetweenNodes(p,p->next) - distanceBetweenNodes(p,path->end);
                 maxNode = p;
             }
        }
        if (maxDist>0) {
            if (distanceBetweenNodes(maxNode,maxNode->next)+
                distanceBetweenNodes(path->head,path->end) >
                distanceBetweenNodes(maxNode, path->end) +
                distanceBetweenNodes(path->head, maxNode->next)) {
                doublylinkedlist *tour = copyList(path, 0, path->countNodes()-1); //create tour from path ##2
                path->end = path->head -> prev;
                tour->flipTwoItems(maxNode->data,path->end->data);
                //cout<<"IM: destroy original path: path in second if"<<endl;
                delete path;
                //path->doublylinkedlist::~doublylinkedlist(); //thisPath.~doublylinkedlist();
                return tour;//return the path of ##2
            }
        }
    }
    return path;
}


//RayOpt takes in a list, ouputs swapped two nodes (4 edges, 2 adjacent pairs)
doublylinkedlist* rayOpt(doublylinkedlist* Path,int NUMITERATIONS) //number of iteration for two Opt)
{
	int num_nodes = Path->countNodes();
	doublylinkedlist* P = copyList(Path, 0, num_nodes-1);
    node *p, *p1, *p3, *temp0, *temp1, *temp2, *temp3;
    int flag,m,n = 0;
    int temp[4];
    float current_distance = 0;
    float best_distance = P->getDistance();
    doublylinkedlist* tempList;
    
    while (n < NUMITERATIONS) {
        // Get pairs
        flag = 0;
        temp[0] = rand() % num_nodes ;
        //cout<<"temp0 is "<<temp[0];
        //cout<<"number of nodes: "<<num_nodes<<endl;
        
        temp[1] = P->getNextIndex(temp[0]); //PROBLEM
        temp[2] = rand() % num_nodes ;
        temp[3] = P->getNextIndex(temp[2]);
        
        //check if the pair are adjacent or if the nodes are the same
        if (temp[0]==temp[2] || temp[0]==temp[3] || temp[1]==temp[2]) flag = 1;

        
        // Save pairs in history
        if (flag==0) {
            n++;
            tempList = copyList(P, 0, num_nodes-1); //create templist ##1
            p = tempList->head;
            while (1) {
                if (p->data==temp[1]) p1 = p;
                else if (p->data==temp[3]) p3 = p;
                p=p->next;
                if (p==tempList->head) break;
            }
            
            temp0 = p1->prev;
            temp1 = p1->next;
            temp2 = p3->prev;
            temp3 = p3->next;
            
            p1->prev->next = p3;
            p1->next->prev = p3;
            p1->prev = temp2;
            p1->next = temp3;
            
            p3->prev->next = p1;
            p3->next->prev = p1;
            p3->prev = temp0;
            p3->next = temp1;
            tempList->end = tempList->head -> prev;
            
            current_distance = tempList->getDistance();
            //delete tempList;
            
            if (current_distance<best_distance) {
                best_distance = current_distance; // Update best distance

                delete P;
                //P->doublylinkedlist::~doublylinkedlist();
                P = copyList(tempList,0,num_nodes-1);

                n = 0;
            }
            delete tempList;
        }
    }
    P->end = P->head ->prev;
    return P;
}

//star opt takes in a dllist, outputs a new dll list that has better distance
doublylinkedlist* starOpt(doublylinkedlist* Path, int K,int NUMITERATIONS)
{
   // srand(3);
   	int num_nodes = Path->countNodes();
	doublylinkedlist* P = copyList(Path, 0, num_nodes-1);
    node *p;
    node *root[K];
    node *temp_root[K];
    int i,k,m,n = 0;
    int temp[2*K];
    int flag,matches;
    float current_distance = 0; // starting distance
    float best_distance = P->getDistance(); // starting best distance
   // vector < vector <int> > pairs;
    doublylinkedlist* tempList;
    //tempList= new doublylinkedlist();
    
    while (n<NUMITERATIONS) {
        // Find K pairs
        k = 0;
        while (k<K) {
            flag = 0;
            temp[2*k] = rand() % num_nodes;
            temp[2*k+1] = P->getNextIndex(temp[2*k]);
           // cout<<"DEBUG: "<<temp[2*k]<<endl;
            
            for (i=0;i<k;i++) { // Make sure pairs have not been seen in current set
                //if (temp[2*i]==temp[2*k] || temp[2*i+1]==temp[2*k] || temp[2*i]==temp[2*k+1] || temp[2*i+1]==temp[2*k+1]) flag = 1;
                if (temp[2*i]==temp[2*k]) flag = 1;
                
            }
            
            if (flag==0) k++;
        }
        
        
        // flag indicates that the apir is legal to move on
        if (flag==0) {

            n++;
            
            tempList = copyList(P,0,num_nodes-1);
            
            // Flip K items
            p = tempList->head;
            i = 0;
            while (1) {
                for (k=0;k<K;k++) {
                    if (p->data==temp[2*k]) {
                        root[i] = p;
                        temp_root[i] = p->next;
                        i++;
                        break;
                    }
                }
                p=p->next;
                if (p==tempList->head) break;
            }
            
            root[0]->next = temp_root[K-2];
            root[1]->next = temp_root[K-1];
            
            for (k=0;k<K;k++) {
                if (k>=2) root[k]->next = temp_root[k-2];
                root[k]->next->prev = root[k];
            }
            
            tempList->end = tempList->head->prev;
            
            // Check distance
            current_distance = tempList->getDistance();

            
            if (current_distance<best_distance) {
                best_distance = current_distance;
                delete P;
                P = copyList(tempList,0,num_nodes-1);
                n = 0;
            }
            delete tempList;
        }
    }

    return P;
}


//TwoOpt takes in a list, ouputs swapped two nodes (4 edges, 2 adjacent pairs)
//Memory Warning: This function always creates a new list, and it will return the newly created one
//Original implementation has history check, now we don't
doublylinkedlist* TwoOpt(doublylinkedlist* Path, int NUMITERATIONS)
{
	
    doublylinkedlist* P;
    int num_nodes = Path->countNodes();
    P = copyList(Path, 0, num_nodes-1);
    
    //  cout<<" P has number of nodes = "<<num_nodes<<endl;
    //P->displayforward(); cout<<endl;
    
    int flag,m,n = 0;
    int temp[4];
    float current_distance = 0;
    float best_distance = P->getDistance();
	//start a history records of which pair 
   // vector<vector<int> > pairs;
    doublylinkedlist* tempList;
    //tempList= new doublylinkedlist();
    
    while (n < NUMITERATIONS){// && n < num_nodes*(num_nodes-3)/2) { //
        // Get pairs
        flag = 0;
        
        temp[0] = rand() % num_nodes ;
    /*  cout<<"temp0 is "<<temp[0];
        cout<<"number of nodes: "<<num_nodes<<endl; */
        
        temp[1] = P->getNextIndex(temp[0]); //PROBLEM
        temp[2] = rand() % num_nodes ;
        temp[3] = P->getNextIndex(temp[2]);
        
        //check if the pair are adjacent or if the nodes are the same
        if (temp[0]==temp[2] || temp[0]==temp[3] || temp[1]==temp[2]) flag = 1;

        
        // Save pairs in history
        if (flag==0) {
            n++;

            
            tempList = copyList(P, 0, num_nodes-1); //create tempList ##1
            
            tempList->flipTwoItems(temp[0], temp[2]);
            tempList->end = tempList->head -> prev;
            
            current_distance = tempList->getDistance();

            
            if (current_distance<best_distance) {
                best_distance = current_distance; // Update best distance
                delete P;
                P = copyList(tempList,0,num_nodes-1);
                n = 0;
            } 
            delete tempList;
        }
        
    }
    
    //printf("debugBEST = %f\n",best_distance);
    P->end = P->head ->prev;

    return P;
}

//end of file
