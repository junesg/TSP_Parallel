#include "doublylinked.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <set>
#include <math.h>

#define MAXDEPTH 2
#define MAXITER 1
#define LISTSIZE 6

using namespace std;

#define LISTSIZE 6
#define NUMITERATIONS 1500

doublylinkedlist ImprovePath(doublylinkedlist P, int depth, set<int>*);
float distanceBetweenNodes(node* n1, node* n2);
//methods in the doublylinked.cpp file
doublylinkedlist copyList(doublylinkedlist P, int start, int end);
doublylinkedlist rayOpt(doublylinkedlist P);
doublylinkedlist starOpt(doublylinkedlist P, int K);
void TSP_LK (doublylinkedlist tour) ;
doublylinkedlist TwoOpt(doublylinkedlist );

int main() {
    doublylinkedlist aList;
    int ind[LISTSIZE] = {0, 1, 2, 3, 4, 5};
    int x[LISTSIZE]   = {0, 0, 1, 1, 2, 2};
    int y[LISTSIZE]   = {10,0, 10,0, 10, 0};
    aList.createList(ind,x,y,6);
    aList.displayforward();
    cout<<endl;
    cout<<"original distance is "<<aList.getDistance()<<endl;

    //TSP_LK (aList);
    cout<<"Number of nodes for alist: "<<aList.countNodes()<<endl;
    cout<<endl<<"final reslult"<<endl;
    //doublylinkedlist bList = starOpt(aList,3);
    doublylinkedlist bList = TwoOpt(aList);
    printf("debug");
    //aList.destroy();
   // aList.~doublylinkedlist();

    cout<<"distance is "<<bList.getDistance()<<endl;

}


float distanceBetweenNodes(node* n1, node* n2){
    return sqrt(pow((n1->x - n2->x),2) +pow((n1->y - n2->y),2));
}


void TSP_LK (doublylinkedlist tour) {
    int iter = 0;
    while (iter < MAXITER) {
       //construct the path
        int countNode = tour.countNodes();
        doublylinkedlist path = copyList(tour,iter%countNode,(iter+countNode-1)%countNode);

 
        set<int> R;
        doublylinkedlist tour2 = ImprovePath(path,1,&R);
        cout<<iter<<" round of improve path we get"<<endl;
        tour2.displayforward();
        cout<<endl<<"with distance:"<<tour2.getDistance()<<endl;
       if (tour2.getDistance() < tour.getDistance()) {
            tour.destroy();
            tour.~doublylinkedlist();
            tour = copyList(tour2, 0, countNode-1);
          //  tour.rearrangeList();
            iter = 0;
        }
        else iter =iter+1;
    }
}


doublylinkedlist ImprovePath(doublylinkedlist path, int depth, set<int>* R){
    //cout<<"in the function"<<endl;
    
    if (depth < MAXDEPTH) {
        //cout<<"enter"<<endl;
        for (node* p = path.head; p != path.end; p = p->next) {
            if ((*R).find(p->data)== (*R).end()) {
               // cout<<"investigate "<<p->data<<" and "<<p->next->data<<endl;
                if (distanceBetweenNodes(p,p->next) > distanceBetweenNodes(p,path.end)) {
                    //if tour length is improved
                  //  cout<<"Tour length is improved"<<endl;
                   // cout<<"debug";
                    if (distanceBetweenNodes(p,p->next)+
                        distanceBetweenNodes(path.head,path.end) >
                        distanceBetweenNodes(p, path.end) +
                        distanceBetweenNodes(path.head, p->next)) {
                     //   cout<<" getting better !"<<endl;
                            doublylinkedlist tour = copyList(path, 0, path.countNodes()-1); //completed tour
                     //       cout<<" Original tour "; tour.displayforward(); cout<<endl;
                            tour.flipTwoItems(p->data, p->next->data, path.end->data, path.head->data);
                           // cout<<"Final tour"<<endl; tour.displayforward();cout<<endl;
                            tour.end = tour.head -> prev;
                        path.destroy();
                        path.~doublylinkedlist();
                            return tour;
                    }
                
                    else {
                        path.flipTwoItems(p->data, p->next->data, path.end->data, path.head->data);
                        path.end = path.head->prev;
                        (*R).insert(p->data);
                     //   cout<< " we are not satisfied : "<<endl;
                        path.displayforward(); cout<<endl;
                        return ImprovePath(path, depth+1, R);
                    }
                }
            }
        }
    }
    else{
        float maxDist=0;
        node* maxNode;
        node* p;
        for (p = path.head; p!= path.end; p = p->next) {
            
             if (distanceBetweenNodes(p,p->next) - distanceBetweenNodes(p,path.end) > maxDist) {
                 maxDist =distanceBetweenNodes(p,p->next) - distanceBetweenNodes(p,path.end);
                 maxNode = p;
             }
        }
        if (maxDist>0) {
           // cout<<"investigate "<<p->data<<" and "<<p->next->data<<"and we got best"<<endl;
            if (distanceBetweenNodes(p,p->next)+
                distanceBetweenNodes(path.head,path.end) >
                distanceBetweenNodes(p, path.end) +
                distanceBetweenNodes(path.head, p->next)) {
                doublylinkedlist tour = copyList(path, 0, path.countNodes()); //completed tour
               // cout<<"before flipping: "; tour.displayforward();cout<<endl;
                tour.flipTwoItems(p->data, p->next->data, path.end->data, path.head->data);
                tour.displayforward();
                return tour;
            }
        }
        else{
            doublylinkedlist tour = copyList(path, 0, path.countNodes()); //completed tour
            tour.flipTwoItems(p->data, p->next->data, path.end->data, path.head->data);
            tour.end = path.head->prev;
            (*R).insert(p->data);
            path.destroy();
            path.~doublylinkedlist();
            //cout<< " we are not satisfied : "<<endl;
            return ImprovePath(tour, depth+1, R);
        }
    }
}


//copies the whole list from the "start" position to the "end" position
doublylinkedlist copyList(doublylinkedlist P, int start, int end)
{
    cout<<"start from "<<start<<" end at "<<end<<endl;
	doublylinkedlist copied;
	if (start > end) {
		return copied;
    }
	if (start < 0 || end > P.countNodes())
		return copied;
	node *pStart, *pEnd;
	node* p = P.head;
	
	for(int count = 0; count <= end; count ++) {
		if (count == start)
			pStart = p;
		if (count == end)
			pEnd = p;
		p = p-> next;
	}
    
	vector<int> ind,x,y;
	int count = 0;
    
	for(p = pStart; p!= pEnd; p = p-> next) {
		//cout<<"data "<<p->data<<"("<<p->x<<","<<p->y<<")"<<endl;
		ind.push_back((int) p->data);
		x.push_back((int) p->x);
		y.push_back((int) p->y);
		count++;
	}
	
    //cout<<"data "<<p->data<<"("<<p->x<<","<<p->y<<")"<<endl;
    ind.push_back(p->data);
    x.push_back(p->x);
    y.push_back(p->y);
    count++;
    
	int* arr = &ind[0];
	int* xPos = &x[0];
	int* yPos = &y[0];
	//cout<<"array done! "<<endl;
    
    //	for(int i = 0; i < count; i++){
    //		cout<<"arr: "<<arr[i]<< " xPos: "<<xPos[i]<<" yPos: "<<yPos[i]<<endl;
    //	}
	copied.createList(arr,xPos,yPos,count);
	return copied;
}


//RayOpt takes in a list, ouputs swapped two nodes (4 edges, 2 adjacent pairs)
doublylinkedlist rayOpt(doublylinkedlist P)
{
    node *p, *p1, *p3, *temp0, *temp1, *temp2, *temp3;
    int num_nodes = P.countNodes();
    
  //  cout<<" P has number of nodes = "<<num_nodes<<endl;
    P.displayforward(); cout<<endl;
    
    int flag,m,n = 0;
    int temp[4];
    float current_distance = 0;
    float best_distance = P.getDistance();
    vector<vector<int> > pairs;
    doublylinkedlist tempList;
    
    while (n < NUMITERATIONS && n < num_nodes*(num_nodes-3)/2) {
        // Get pairs
        flag = 0;
        
        temp[0] = rand() % num_nodes ;
        cout<<"temp0 is "<<temp[0];
        cout<<"number of nodes: "<<num_nodes<<endl;
        
        temp[1] = P.getNextIndex(temp[0]); //PROBLEM
        temp[2] = rand() % num_nodes ;
        temp[3] = P.getNextIndex(temp[2]);
        
        //check if the pair are adjacent or if the nodes are the same
        if (temp[0]==temp[2] || temp[0]==temp[3] || temp[1]==temp[2]) flag = 1;
        //check if this pair is in the history or not
        for (m=0;m<n;m++) {
            if (temp[0]==pairs[m][0] && temp[2]==pairs[m][1]) flag = 1;
            if (temp[2]==pairs[m][0] && temp[0]==pairs[m][1]) flag = 1;
        }
        
        // Save pairs in history
        if (flag==0) {
            pairs.resize(n+1);
            for (m=0;m<=n;m++) {
                pairs.at(m).resize(2);
            }
            pairs[n][0] = temp[0];
            pairs[n][1] = temp[2];
            
            tempList = copyList(P, 0, num_nodes-1);
            
           // tempList.flipNodes(temp[1], temp[3]);
        
            
            p = tempList.head;
            while (1) {
                if (p->data==temp[1]) p1 = p;
                else if (p->data==temp[3]) p3 = p;
                p=p->next;
                if (p==tempList.head) break;
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
            tempList.end = tempList.head -> prev;
            
            current_distance = tempList.getDistance();
            tempList.~doublylinkedlist();
            
            printf("%d. Trying flip: (%d %d)(%d %d), distance = %f\n",n,temp[0],temp[1],temp[2],temp[3],current_distance);
            
            if (current_distance<best_distance) {
                best_distance = current_distance; // Update best distance
                printf("BEST = %f\n",best_distance);
                cout<<"DIsplay: "<<endl;
                tempList.displayforward();cout<<endl;
               // tempList.end = tempList.head -> prev;
                P.displayforward(); cout<<endl;
                P.destroy();
                P.~doublylinkedlist();
                P = copyList(tempList,0,num_nodes-1);
                cout<<"DIsplay templist: "<<endl;
                //P.displayforward(); cout<<endl;
                printf("debug");
                n = 0;
            } else n++;
        }
    }
    printf("debugBEST = %f\n",best_distance);
    P.end = P.head ->prev;
    P.displayforward();
    printf("\n");
    return P;
}


doublylinkedlist starOpt(doublylinkedlist P, int K)
{
   // srand(3);
    node *p;
    node *root[K];
    node *temp_root[K];
    int num_nodes = P.countNodes();
    int i,k,m,n = 0;
    int temp[2*K];
    int flag,matches;
    float current_distance = 0; // starting distance
    float best_distance = P.getDistance(); // starting best distance
    vector < vector <int> > pairs;
    doublylinkedlist tempList;
    
    while (n<NUMITERATIONS) {
        // Find K pairs
        k = 0;
        while (k<K) {
            flag = 0;
            temp[2*k] = rand() % num_nodes;
            temp[2*k+1] = P.getNextIndex(temp[2*k]);
           // cout<<"DEBUG: "<<temp[2*k]<<endl;
            
            for (i=0;i<k;i++) { // Make sure pairs have not been seen in current set
                //if (temp[2*i]==temp[2*k] || temp[2*i+1]==temp[2*k] || temp[2*i]==temp[2*k+1] || temp[2*i+1]==temp[2*k+1]) flag = 1;
                if (temp[2*i]==temp[2*k]) flag = 1;
                
            }
            
            if (flag==0) k++;
        }
        
        // Make sure pairs have not been seen in history
        flag = 0;
        for (m=0;m<n;m++) {
            matches = 0;
            for (k=0;k<K;k++) {
                for (i=0;i<K;i++) {
                    //matches += temp[2*i]==pairs[m][2*k];
                }
            }
            if (matches==K) {
                flag = 1;
                break;
            }
        }
        
        // Test pairs
        if (flag==0) {
            //            // Save pairs in history
            //            pairs.resize(n+1);
            //            for (m=0;m<=n;m++) {
            //                pairs.at(m).resize(2*K);
            //            }
            //            for (k=0;k<K;k++) {
            //                pairs[n][2*k] = temp[2*k];
            //                pairs[n][2*k+1] = temp[2*k+1];
            //            }
            n++;
            
            tempList = copyList(P,0,num_nodes-1);
            
            // Flip K items
            p = tempList.head;
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
                if (p==tempList.head) break;
            }
            
            root[0]->next = temp_root[K-2];
            root[1]->next = temp_root[K-1];
            
            for (k=0;k<K;k++) {
                if (k>=2) root[k]->next = temp_root[k-2];
                root[k]->next->prev = root[k];
            }
            
            tempList.end = tempList.head->prev;
            
            // Check distance
            current_distance = tempList.getDistance();
            printf("%d. Trying flip, distance = %f\n",n,current_distance);
            for (k=0;k<K;k++) printf("%d ",temp[2*k]);
            printf("\n");
            
            if (current_distance<best_distance) {
                best_distance = current_distance;
                printf("BEST = %f\n",best_distance);
                P.destroy();
                P.~doublylinkedlist();
                P = copyList(tempList,0,num_nodes-1);
                n = 0;
            }
            tempList.destroy();
            tempList.~doublylinkedlist();
        }
    }
    printf("BEST = %f\n",best_distance);
    P.displayforward();
    printf("\n");
    return P;
}




//TwoOpt takes in a list, ouputs swapped two nodes (4 edges, 2 adjacent pairs)
doublylinkedlist TwoOpt(doublylinkedlist P)
{
    //node *p, *p1, *p3, *temp0, *temp1, *temp2, *temp3;
    int num_nodes = P.countNodes();
    
    //  cout<<" P has number of nodes = "<<num_nodes<<endl;
    P.displayforward(); cout<<endl;
    
    int flag,m,n = 0;
    int temp[4];
    float current_distance = 0;
    float best_distance = P.getDistance();
    vector<vector<int> > pairs;
    doublylinkedlist tempList;
    
    while (n < NUMITERATIONS && n < num_nodes*(num_nodes-3)/2) {
        // Get pairs
        flag = 0;
        
        temp[0] = rand() % num_nodes ;
        cout<<"temp0 is "<<temp[0];
        cout<<"number of nodes: "<<num_nodes<<endl;
        
        temp[1] = P.getNextIndex(temp[0]); //PROBLEM
        temp[2] = rand() % num_nodes ;
        temp[3] = P.getNextIndex(temp[2]);
        
        //check if the pair are adjacent or if the nodes are the same
        if (temp[0]==temp[2] || temp[0]==temp[3] || temp[1]==temp[2]) flag = 1;
        //check if this pair is in the history or not
        for (m=0;m<n;m++) {
            if (temp[0]==pairs[m][0] && temp[2]==pairs[m][1]) flag = 1;
            if (temp[2]==pairs[m][0] && temp[0]==pairs[m][1]) flag = 1;
        }
        
        // Save pairs in history
        if (flag==0) {
            pairs.resize(n+1);
            for (m=0;m<=n;m++) {
                pairs.at(m).resize(2);
            }
            pairs[n][0] = temp[0];
            pairs[n][1] = temp[2];
            
            tempList = copyList(P, 0, num_nodes-1);
            
           tempList.flipTwoItems(temp[0], temp[1], temp[2], temp[3]);
            tempList.end = tempList.head -> prev;
            
            current_distance = tempList.getDistance();
            tempList.~doublylinkedlist();

            printf("%d. Trying flip: (%d %d)(%d %d), distance = %f\n",n,temp[0],temp[1],temp[2],temp[3],current_distance);
            
            if (current_distance<best_distance) {
                best_distance = current_distance; // Update best distance
                printf("BEST = %f\n",best_distance);
                cout<<"DIsplay: "<<endl;
                tempList.displayforward();cout<<endl;
                // tempList.end = tempList.head -> prev;
                P.displayforward(); cout<<endl;
                P.destroy();
                P.~doublylinkedlist();
                P = copyList(tempList,0,num_nodes-1);
                cout<<"DIsplay templist: "<<endl;
                //P.displayforward(); cout<<endl;
                printf("debug");
                n = 0;
            } else n++;
        }
    }
    printf("debugBEST = %f\n",best_distance);
    P.end = P.head ->prev;
    P.displayforward();
    printf("\n");
    return P;
}