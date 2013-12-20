#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "doublylinked.h"

using namespace std;
#define LISTSIZE 6
#define NUMITERATIONS 20


//void rayOpt(doublylinkedlist P);
//void starOpt(doublylinkedlist P, int K);

//for debugging

/*
int main(){
    doublylinkedlist aList;doublylinkedlist bList;
    int ind[LISTSIZE] = {0,1,2,3,4,5};
    int x[LISTSIZE]   = {0,0,1,1,2,2};
    int y[LISTSIZE]   = {10,0,10,0,10,0};
    aList.createList(ind,x,y,LISTSIZE);
   // bList = copyList(aList, 2, 5);
    
    printf("RAY OPT\n");
    aList.displayforward();
    printf("\n");
    //rayOpt(aList);
    
   // bList.~doublylinkedlist();
    aList.destroy();
    aList.~doublylinkedlist();
    //printf("STAR OPT\n");
    //aList.displayforward();
    //printf("\n");
    //starOpt(aList,3);
}
*/


/*
 *Copies the whole list from the "start" position to the "end" position
 *for example to copy from 0 to num_nodes-1, start = 0, num_nodes-1=end
 *returns a different doublylinkedlist
 */
doublylinkedlist* copyList(doublylinkedlist* P, int start, int end)
{
    doublylinkedlist* copied;
    copied = new doublylinkedlist();
    
    cout<<"Start COPYLIST"<<endl;
	if (start > end) { //start has to come before end
		return copied; //return empty list
    }
    if (start < 0 || end > P->countNodes()) //check bounds of start and end
		return copied; //return empty list
    
	int count = end-start +1; //number of nodes to be added

	node *pStart, *pEnd;
	node* p = P->head;

	//int count = 0;
	for(int icount = 0; icount <= end; icount ++) {
		if (icount == start)
			pStart = p;
		if (icount == end)
			pEnd = p;
		p = p-> next;
	}

	int arr[count], xPos[count], yPos[count];
    //for loop to store the data
	int iter = 0;
	for(p = pStart; p!= pEnd; p = p-> next) {
		arr[iter] = (p->data);
		xPos[iter] = ((int) p->x);
		yPos[iter]=((int) p->y);
		iter ++;
	}
	
		arr[iter] = (p->data);
		xPos[iter] = ((int) p->x);
		yPos[iter]=((int) p->y);
	copied->createList(arr,xPos,yPos,count);
    cout<<"getting out of COPYLIST"<<endl;
	return copied;
}




/*create list inputs an array of indices, these indices will be stored as data
 *it also stores and the x,y positions of the points via input x and y arrays.
 *the n defines the number of nodes in the arry
 * NOTE that the list is not a cycle.
 */
void doublylinkedlist::createList(int ind[], int xPos[], int yPos[], int n)
{
	node *q;  //just a pointer
	//node *p = (node*)malloc(sizeof(node)); //the first node
    node *p = new node(ind[0],(float)xPos[0],(float)yPos[0]);
	//p->data = ind[0];
    //p->x = xPos[0];
    //p->y = yPos[0];
	head = p;
	p-> next = NULL;
	p-> prev = NULL;
    
	for (int i=1; i<n; i++){
		q = p;
		p = new node(ind[i], (float)(xPos[i]),(float)(yPos[i]));
        q->next = p;
		//p->data = ind[i];
        //p->x = xPos[i];
        //p->y = yPos[i];
		p->next = NULL;
		p->prev = q;
	}
    end = p;
    end -> next = head;
    head -> prev = end;
}


doublylinkedlist::~doublylinkedlist(){
    if (head!=NULL) {
        node* current = head->next;
        while( current != head ) {
            node* next = current->next;
            delete current;
            // delete current;
            current = next;
        }
        delete(current);
        cout<<"finished"<<endl;
    }
}

void doublylinkedlist::rearrangeList(int start){
	//Now loop through to make sure head is always the 0_th item

	node* p=head;
	if (p->data != start) {
	    p = p-> next;
        while (p->data!=start && p!= head) {
 	       p = p->next;
 	   }
	    head = p;
	}
    
    //also change the prev pointers
    node* q = head;
    for (p=head->next; p!=head; ) {
        p->prev = q;
        q = p;
		p = p->next;
    }
    end = q;
	head ->prev = q;
		
}

//takes in the data of a node, get the next node's data
int doublylinkedlist::getNextIndex(int n)
{
    node *p = head;
    int m=0;
    while (1) {
        if (p->data==n) {
            m = p->next->data;
            break;
        } else {
            p = p->next;
        }
        if (p==head) break;
    }
    return m;
}

//This function counts the number of nodes in the linkedlist
int doublylinkedlist::countNodes() {
    if (head == NULL) {
        return 0;
    }
    
    node *p = head;
	//end = head-> prev;
    p = p->next;
    int num_nodes = 1;

    while(p!=head) {
        p = p->next;
		num_nodes++;
    }
 //   cout<<"Number of nodes; "<<num_nodes<<endl;
    return num_nodes;
}

//This function returns the distance of a tour of the linkedlist
float doublylinkedlist::getDistance() {
    node *p;
    float distance = 0;
    p = head;
    while (1) {
        distance += sqrt(pow(p->x - p->next->x,2)+pow(p->y - p->next->y,2));
        p=p->next;
        if (p==head) break;
    }
   // distance +=sqrt(pow(head->x - end->x,2)+pow(head->y - end->y,2));
    return distance;
}

//Insert item after key
void doublylinkedlist::insertAfter(int item, int x, int y, int key) {
	node *q=head;
	int head2 = 0;
	//find q first
	while(head2<1) {
		if(q->data==key) break;
		q=q->next;
		if(q == head) head2++;
	}
    
	//if the q is not found.
	if(q==head && head2>1) {
		std::cout<<"key "<< key <<" not found"<<std::endl;
		return;
	}
    
	node *p=new node(item, (float)x, (float)y);
	p->next=q->next;
	p->prev=q;
	q->next=p;
}

//displaying list nodes in forward direction
void doublylinkedlist::displayforward() {
    node *p=head;
    if (head == NULL) {
        printf("Can not display\n");
        return;
    }
	while(1) {
        printf("%d ",p->data);
        //printf("%d (%d,%d), ",p->data,p->x,p->y);
        p=p->next;
        if (p==NULL || p==head) break;
	}
	cout<<endl;
}

//displaying list nodes in reverse direction
void doublylinkedlist::displaybackward() {
	node *p=head->prev;
	while(1) {
		printf("%d ",p->data);
		p=p->prev;
        if (p==head->prev) break;
	}
}

//fliping two edges according to 2opt
//This function takes in (pair11) and (pair21)
//It flips (pair11,pair12) and (pair21,pair22) 
void doublylinkedlist:: flipTwoItems(const int pair11, const int pair21){
	cout<<"Before flipping "<<pair11<< " and "<<pair21<<endl; displayforward();
    
	node *p11, *p12, *p21, *p22;
	p11 = p12 = p21 = p22 = NULL;
	int pair12 = getNextIndex(pair11);
	int pair22 = getNextIndex(pair21);
    
    if (pair12 == pair11 || pair12 == pair21 || pair12 == pair22 || pair11 == pair22 || pair22==pair21) {
        cout<<"Command error:"; printf(" pair11=%d, pair12=%d, pair21=%d, pair22=%d\n", pair11, pair12,pair21, pair22);
        return;
    }
    
    
	node *p; //pointer for iteration
	end = head->prev;
	p = head;
	//find the corresponding vectors
	while (p11 == NULL || p12 == NULL || p21 == NULL || p22 == NULL|| p!= end){
		if (p->data == pair11)
			p11 = p;
		if (p-> data == pair12)
			p12 = p;
		if (p-> data == pair21)
			p21 = p;
		if (p-> data == pair22)
			p22 = p;

		p = p->next;
	}
	if(end-> data == pair11)
		p11 = end;
	if(end-> data == pair12)
		p12 = end;
	if(end-> data == pair21)
		p21 = end;
	if(end-> data == pair22)
		p22 = end;
	
	//check if the pointers are all assigned
	if(p11 == NULL || p12 == NULL || p21 == NULL || p22 == NULL) {
		cout<<"Inpute parameter error, parameter out of bound."<<endl;
		return;
	}
	
	//Now we swap the edges
	if(p11->next != p12 && p11->prev != p12){
		cout<<"Input parameter error: parameters 1 and 2 do not form an edge"<<endl;
		return;
	}
	if(p21->next != p22 && p21 -> prev != p22) {
		cout<<"Input parameter error: parameters 3 and 4 do not form an edge"<<endl;
		return;
	}
	
	//first, we want to know the order of the pointers
	//this guarantees that p11->p12; p21 -> p22
	if(p11-> prev== p12){
		node *temp = p11;
		p11 = p12;
		p12 = temp;
	}
	if(p21-> prev== p22){
		node *temp = p21;
		p21 = p22;
		p22 = temp;
	}
	
	//Now swap the edges and reverse the prev/next pointers of in half of the paths
	node* temp = p12->next;
	node* temp2 = p21 -> prev;
	p11-> next = p21;
	p21-> prev = p11;
	p12 -> next = p22;
	p22 -> prev = p12;
	p21 -> next = temp2;
	p12 -> prev= temp;
	p = temp;
	while(p!= p21){
        
		node *ptp, *ptn;
		ptp = p->prev;
		ptn = p ->next;
		p -> prev = ptn;
		p -> next = ptp;
		p = p-> prev;
        
	}
    end = head->prev;

	cout<<"After flipping"<<endl; displayforward(); cout<<"END"<<endl;

}


 
//delete the ith node of the list 
void doublylinkedlist::deleteNode(int i){
	node *p;
	p = head;
    if (p->data == i) { //we delete the head
        head = p->next;
    }
    else {
		p= p->next;
		while(p!=head){
			if(p-> data !=i)
				p = p->next;
			else break;			
		}
	}
    p -> prev -> next = p -> next;
    p -> next -> prev = p -> prev;
    delete p;

}

//returns an integer array of the array
void doublylinkedlist::extractIndices(int* arr){
    arr[0] = ((int)head->data);
    node* p = head-> next;
    int count = 1;
    while (p!= head) {
        arr[count]=((int)p->data);
        p = p->next;
        count ++;
    }
}

//append list L2 to this function
void doublylinkedlist::appendList(doublylinkedlist l2){
    if (head != NULL ) {
        end = head-> prev;
        int num = l2.countNodes();
        //cout<<num<<" == number of l2"<<endl;
        for (int i=0; i< num; i++) {
            int arr[3];
            l2.getContentAt(i, arr);
            int lastKey = end ->data;
            //cout<<"Got data " << arr[0]<<","<<arr[1]<<","<<arr[2]<<endl;
            insertAfter(arr[0], arr[1], arr[2], lastKey);
			head = end->next->next;
        }
    }
    else {
        int num = l2.countNodes();
        int arr[3];
        int i = 0;
        int index[1]; int x[1]; int y[1];
        
        l2.getContentAt(i, arr);
        
        index[0] = arr[0]; x[0] = arr[1]; y[0] = arr[2];
        createList(index, x, y,1 ); //now the list is created
        for (int i=1; i< num; i++) {
            l2.getContentAt(i, arr);
            int lastKey = head-> prev->data;
            //cout<<"Got data " << arr[0]<<","<<arr[1]<<","<<arr[2]<<endl;
            insertAfter(arr[0], arr[1], arr[2], lastKey);
        }
    }
}

//get the content of key, x and y at the nodeCount (0 to numNodes-1)
void doublylinkedlist::getContentAt(int nodeCount, int arr[3]){
    node* p = head;
    for (int i = 0; i < nodeCount; i++) {
        p= p-> next;
    }
    arr[0] = p-> data;
    arr[1] = p-> x;
    arr[2] = p-> y;
}



//This function flips four edges with center n1, n3
// for instance, flipping 3,7 from 0-1-2-...-7
// then we end up changing edges (2,3)(3,4) and edges (6,7)(7,0)
//Done by TR.

void doublylinkedlist::flipNodes(int n1, int n3) {
    
    node *p, *p1, *p3, *temp0, *temp1, *temp2, *temp3;
    p = head;
    while (1) {
        if (p->data==n1) p1 = p;
        else if (p->data==n3) p3 = p;
        p=p->next;
        if (p==head) break;
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
}




void doublylinkedlist::flipKItems(int K, int *temp) {
    node *p = head;
    node *root[K];
    node *temp_root[K];
    int j,k;
    j = 0;
    while (1) {
        for (k=0;k<K;k++) {
            if (p->data==temp[2*k]) {
                root[j] = p;
                temp_root[j] = p->next;
                j++;
                break;
            }
        }
        p=p->next;
        if (p==head) break;
    }
    
    root[0]->next = temp_root[K-2];
    root[1]->next = temp_root[K-1];
    
    for (k=0;k<K;k++) {
        if (k>=2) root[k]->next = temp_root[k-2];
        root[k]->next->prev = root[k];
    }
    end = head -> prev;
}
    
//compares this current list with another input list aList
//starts from the head and examines every data afterwards
bool doublylinkedlist::compareList(doublylinkedlist aList){
	//aList.rearrangeList(0);
	//rearrangeList(0);
	node* thisp= head;
	
    aList.rearrangeList(head->data);
    node* thatp = aList.head;
    
	if(thisp->data != thatp->data)
		return false;
	thisp = thisp->next;
	thatp = thatp->next;
	while(thisp != head && thatp != aList.head) {
		//cout<<" "<<thisp->data<<" compared to "<<thatp->data<<endl;
		if(thisp->data != thatp->data)
			return false;
		thisp = thisp -> next; thatp = thatp -> next;
    }
    return true;
}
