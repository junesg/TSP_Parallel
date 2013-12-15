#ifndef DLL_H
#define DLL_H
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

struct node {   //node struct, has data in it and has pointers to previous and next nodes
    int data;
    float x;
    float y;
    node *prev;
    node *next;
};

class doublylinkedlist
{
public:
    node *head;
    node *end;
    doublylinkedlist() //constructor
    {
        head = NULL;
        end = NULL;
    }
    
   	void createList(int*,int*,int*,int);
    void destroy();

    doublylinkedlist removeParts(doublylinkedlist t1, doublylinkedlist p2);
	void insertAfter(int,int, int, int);
	void displayforward();
	void displaybackward();
    int countNodes();
    int getNextIndex(int);
    void findKPairs(int, vector<vector<int> > *);
    void findPairs(vector<vector<int> > *);
    float getDistance();
	void flipTwoItems(const int,const int,const int,const int);
    void flipNodes(int,int);
	void showDistance();
	void deleteNode(int);
	void rearrangeList();
    void rayOpt(int);
    void starOpt(int,int);
    void flipKItems(int, int*);
    void extractIndices(int*);
    void appendList(doublylinkedlist);
    void getContentAt(int, int arr[3]);
};




//create list inputs an array of indices, these indices will be stored as data
//it also stores and the x,y positions of the points via input x and y arrays.
//the n defines the number of nodes in the arry
//*** NOTE that the list is not a cycle.
void doublylinkedlist::createList(int* ind, int* x, int* y, int n)
{
	node *q;  //just a pointer
	node *p = new node; //the first node
	p->data = ind[0];
    p->x = x[0];
    p->y = y[0];
	head = p;
	p-> next = NULL;
	p-> prev = NULL;
	for (int i=1; i<n; i++){
		q = p;
		p = p->next = new node;
		p->data = ind[i];
        p->x = x[i];
        p->y = y[i];
		p->next = NULL;
		p->prev = q;
	}
    end = p;
    p -> next  = head;
    head -> prev = end;

}



void doublylinkedlist::destroy(){
    node *p, *q;
    p = head;
    q = head;
    for (p = head -> next; p != end; p = p->next) {
        q = p->prev;
        delete q;
    }
    delete p;
}





void doublylinkedlist::rearrangeList(){
	//Now loop through to make sure head is always the 0_th item
	node* p=head;
	if (p->data != 0) {
	    p = p-> next;
	    while (p->data!=0) {
 	       p = p->next;
 	   }
	    head = p;
	}
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
                //printf("(%2d,%2d)",temp[2*k],temp[2*k+1]);
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
    
}



 
void doublylinkedlist::findKPairs(int K, vector<vector<int> > *pairs) {
    int num_nodes = countNodes();

    int i,k,m,n = 0;
    int temp[2*K];
    int flag;
    int matches;

    while (n < 10) {
        k = 0;
        while (k<K) {
            flag = 0;
            
            temp[2*k] = rand() % num_nodes;
            temp[2*k+1] = getNextIndex(temp[2*k]);
        
            // Make sure pairs have not been seen in current set
            for (i=0;i<k;i++) {
                if (temp[2*i]==temp[2*k] || temp[2*i+1]==temp[2*k] || temp[2*i]==temp[2*k+1] || temp[2*i+1]==temp[2*k+1]) flag = 1;
            }
        
            if (flag==0) k++;
        }
        
        // Make sure valid pair set has not been seen in history - NOTE: DEBUG THIS!!!!!!!!!!!!!!!
        flag = 0;
        for (int m = 0;m < n;m++) {
            matches = 0;
            for (int k=0;k < K;k++) {
                for (int i=0;i < K;i++) {
                    matches += temp[2*i]==(*pairs)[m][2*k];
                }
            }
            if (matches==K) {
                flag = 1;
                break;
            }
        }
        
        // Record pair
        if (flag==0) {
            (*pairs).resize(n+1);
            for (int m=0;m<=n;m++) {
                (*pairs).at(m).resize(2*K);
            }
            for (k=0;k<K;k++) {
                (*pairs)[n][2*k] = temp[2*k];
                (*pairs)[n][2*k+1] = temp[2*k+1];
                printf("(%2d,%2d) %d",(*pairs)[n][2*k],(*pairs)[n][2*k+1],n);
            }
            printf("\n");
            n++;
        }
    }
}


 
void doublylinkedlist::findPairs(vector<vector<int> > *pairs) {
    int num_nodes = countNodes();
    int i,j,k,l,m,n = 0;
    int temp[4];
    int flag;
    for (i=0;i<num_nodes;i++) {
        for (j=0;j<num_nodes;j++) {
            flag = 0;
            temp[0] = i;
            temp[1] = getNextIndex(temp[0]);
            temp[2] = j;
            temp[3] = getNextIndex(temp[2]);
            
            // Make sure no two nodes are the same
            for (k=0;k<4;k++) {
                for (l=0;l<4;l++) {
                    if (k!=l && temp[k]==temp[l]) flag = 1;
                }
            }
            
            // Make sure pair hasn't been seen before
            for (m=0;m<n;m++) {
                if (temp[0]==(*pairs)[m][0] && temp[1]==(*pairs)[m][1] && temp[2]==(*pairs)[m][2] && temp[3]==(*pairs)[m][3]) flag = 1;
                if (temp[0]==(*pairs)[m][2] && temp[1]==(*pairs)[m][3] && temp[2]==(*pairs)[m][0] && temp[3]==(*pairs)[m][1]) flag = 1;
            }
            
            // Record pair
            if (flag==0) {
                (*pairs)[n][0] = temp[0];
                (*pairs)[n][1] = temp[1];
                (*pairs)[n][2] = temp[2];
                (*pairs)[n][3] = temp[3];
                n++;
            }
        }
    }
}



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
    node *p = head;
    p = p->next;
    int num_nodes = 1;
    while(p!=head) {
        num_nodes++;
        p = p->next;
    }
    return num_nodes;
}


 
float doublylinkedlist::getDistance() {
    node *p;
    float distance = 0;
    p = head;
    while (1) {
        distance += sqrt(pow(p->x - p->next->x,2)+pow(p->y - p->next->y,2));
        p=p->next;
        if (p==head) break;
    }
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
    
	node *p=new node;
	p->data=item;
    p->x=x;
    p->y=y;
	p->next=q->next;
	p->prev=q;
	q->next=p;
}

//diaplaying list nodes in forward direction
 
void doublylinkedlist::displayforward() {
	node *p=head;
	while(1) {
        printf("%d ",p->data);
        //printf("%d (%d,%d), ",p->data,p->x,p->y);
        p=p->next;
        if (p==NULL || p==head) break;
	}
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
//This function takes in (pair12,pair22) and (pair21,pair22) and flips them
 
void doublylinkedlist:: flipTwoItems(const int pair11, const int pair12, const int pair21, const int pair22){
	node *p11, *p12, *p21, *p22;
	p11 = p12 = p21 = p22 = NULL;
	node *p; //pointer for iteration
	
	node *tail = head->prev;
	//check to see if input vectors are repetitive
	if (pair11 == pair12 || pair21 == pair22 || pair11==pair21|| pair12==pair22){
		cout<<"Input parameter error, there is a duplicate in input value."<<endl;
		return;
	}
	p = head;
	//find the corresponding vectors
	while (p11 == NULL || p12 == NULL || p21 == NULL || p22 == NULL|| p!= tail){
		if (p->data == pair11){
			p11 = p;
		}
		if (p-> data == pair12){
			p12 = p;
		}
		if (p-> data == pair21)
			p21 = p;
		if (p-> data == pair22)
			p22 = p;
		p = p->next;
	}
	if(tail-> data == pair11){
		p11 = tail;
	}
	if(tail-> data == pair12){
		p12 = tail;
	}
	if(tail-> data == pair21){
		p21 = tail;
	}
	if(tail-> data == pair22){
		p22 = tail;
	}
	
	
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
}


 
void doublylinkedlist::showDistance() {
    node *p;
    float distance = 0;
    p = head;
    while (1) {
        distance += sqrt(pow(p->x - p->next->x,2)+pow(p->y - p->next->y,2));
		cout<<"From "<<p->data<<" <"<<p->x<<","<<p->y<<"> to "<<p->next->data<<" <"<<p->next->x<<","<<p->next->y<<"> distance= "<<distance<<endl;
        p=p->next;
        if (p==head) break;
    }
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
        node* end = head-> prev;
        int num = l2.countNodes();
        //cout<<num<<" == number of l2"<<endl;
        for (int i=0; i< num; i++) {
            int arr[3];
            l2.getContentAt(i, arr);
            int lastKey = head-> prev->data;
            //cout<<"Got data " << arr[0]<<","<<arr[1]<<","<<arr[2]<<endl;
            insertAfter(arr[0], arr[1], arr[2], lastKey);
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




#endif
