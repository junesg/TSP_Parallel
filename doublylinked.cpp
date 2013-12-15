#include <iostream>
#include <stdio.h>

#include "doublylinked.h"

using namespace std;
#define LISTSIZE 6
#define NUMITERATIONS 20

doublylinkedlist copyList(doublylinkedlist P, int start, int end);
void rayOpt(doublylinkedlist P);
void starOpt(doublylinkedlist P, int K);

//for debugging
int main(){
    doublylinkedlist aList;doublylinkedlist bList;
    int ind[LISTSIZE] = {0,1,2,3,4,5};
    int x[LISTSIZE]   = {0,0,1,1,2,2};
    int y[LISTSIZE]   = {10,0,10,0,10,0};
    aList.createList(ind,x,y,LISTSIZE);
    bList = copyList(aList, 2, 5);
    
    printf("RAY OPT\n");
    aList.displayforward();
    printf("\n");
    rayOpt(aList);
    
    printf("STAR OPT\n");
    aList.displayforward();
    printf("\n");
    starOpt(aList,3);
}



//copies the whole list from the "start" position to the "end" position
doublylinkedlist copyList(doublylinkedlist P, int start, int end)
{
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
void rayOpt(doublylinkedlist P)
{
    node *p, *p1, *p3, *temp0, *temp1, *temp2, *temp3;
    int num_nodes = P.countNodes();
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
        temp[1] = P.getNextIndex(temp[0]);
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
            
            tempList = copyList(P, 0,num_nodes-1);
            

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
            
            current_distance = tempList.getDistance();
            tempList.~doublylinkedlist();
            
            printf("%d. Trying flip: (%d %d)(%d %d), distance = %f\n",n,temp[0],temp[1],temp[2],temp[3],current_distance);
            
            if (current_distance<best_distance) {
                best_distance = current_distance; // Update best distance
                printf("BEST = %f\n",best_distance);
                P.~doublylinkedlist();
                P = copyList(tempList,0,num_nodes-1);
                n = 0;
            } else n++;
        }
    }
    printf("BEST = %f\n",best_distance);
    P.displayforward();
    printf("\n");
}

void starOpt(doublylinkedlist P, int K)
{
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
            
            // Check distance
            current_distance = tempList.getDistance();
            printf("%d. Trying flip, distance = %f\n",n,current_distance);
            
            if (current_distance<best_distance) {
                best_distance = current_distance;
                printf("BEST = %f\n",best_distance);
                P.~doublylinkedlist();
                P = copyList(tempList,0,num_nodes-1);
                n = 0;
            }
            tempList.~doublylinkedlist();
        }
    }
    printf("BEST = %f\n",best_distance);
    P.displayforward();
    printf("\n");
}