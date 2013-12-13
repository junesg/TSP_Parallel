#ifndef DLL_H
#define DLL_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;


template <class T>
class doublylinkedlist
{
private:
    struct node {   //node struct, has data in it and has pointers to previous and next nodes
        T data;
        float x;
        float y;
        node *prev;
        node *next;
    };
    node *head;
    
public:
    doublylinkedlist() //constructor
    {
        head = NULL;
    }
    
	void createList(T*,T*,T*,int);
	doublylinkedlist<T> copyList(int start, int end);
    doublylinkedlist<int> removeParts(doublylinkedlist<int> t1, doublylinkedlist<int> p2);
	void insertAfter(T,T,T,T);
	void displayforward();
	void displaybackward();
    void rayOpt();
    int countNodes();
    int getNextIndex(int);
    void findKPairs(int, vector<vector<int> > *);
    void findPairs(vector<vector<int> > *);
    float getDistance();
	void flipTwoItems(const T,const T,const T,const T);
    void flipNodes(int,int);
	void showDistance();
	void deleteNode(T);
	void rearrangeList();
    void starOpt(int K,int LISTSIZE);
    void flipKItems(int, int*);
    void extractIndices(int*);
    void appendList(doublylinkedlist<int>);
    void getContentAt(T, int arr[3]);
};


//create list inputs an array of indices, these indices will be stored as data
//it also stores and the x,y positions of the points via input x and y arrays.
//the n defines the number of nodes in the arry
template <class T>
void doublylinkedlist<T>::createList(T* ind, T* x, T* y, int n)
{
	node *q;  //just a pointer
	node *p = new node;  //the first node
	p->data = ind[0];
    p->x = x[0];
    p->y = y[0];
	head = p;
	p-> next = NULL;
	p-> prev = NULL;
	for (int i=1; i<n; i++){
		q = p;
		p = p->next=new node;
		p->data = ind[i];
        p->x = x[i];
        p->y = y[i];
		p->next = NULL;
		p->prev = q;
	}
	p -> next  = head; //TSP is a circle.
	head -> prev = p;

}

template <class T>
void doublylinkedlist<T>::rearrangeList(){
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


//copies the whole list from the "start" position to the "end" position
template <class T>
doublylinkedlist<T> doublylinkedlist<T>::copyList(int start, int end)
{
	doublylinkedlist<int> copied;
	if (start == end)
		return copied;
	if (start < 0 || end > countNodes())
		return copied;
	node *pStart, *pEnd;
	node* p = head;
	
	for(int count = 0; count <= end; count ++) {
		if (count == start)
			pStart = p;
		if (count == end)
			pEnd = p;
		p = p-> next;
	}

/*
	if (p -> data == start)
		pStart = p;
	if (p -> data == end)
		pEnd = p;
	p = p->next;
	while(p != head) {
		if (p -> data == start){
			pStart = p;
			//cout<<"Start "<<start<<","<<pStart->data<<" found"<<endl;
		}
		if (p -> data == end){
			pEnd = p;
			//cout<<"End "<<pEnd<<" found"<<endl;
		}
		p = p->next;
	}
*/	
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

	T* arr = &ind[0];
	T* xPos = &x[0];
	T* yPos = &y[0];
//	cout<<"array done! "<<endl;
	
	for(int i = 0; i < count; i++){
		cout<<"arr: "<<arr[i]<< " xPos: "<<xPos[i]<<" yPos: "<<yPos[i]<<endl;		
	}
//	cout<< "finish"<<endl;
	copied.createList(arr,xPos,yPos,count);
	return copied;
}


template<class T>
void doublylinkedlist<T>::rayOpt()
{
    int num_nodes = countNodes();
    
    // Find number of possible pairs
    int num_pairs = num_nodes*(num_nodes-3)/2;

    int i; // counters
    int K = 2;
    vector < vector <int> > pairs;
    pairs.resize(num_pairs);
    for (i=0;i<num_pairs;i++) pairs[i].resize(K);

    float current_distance;
    float best_distance;
    float first_best_distance;
    int index_best;
    
    int trial = 0;
    while (1) {
        // Find possible pairs
        
        findPairs(&pairs);
        
        //findKPairs(2,&pairs);
        
        // Display linked list
        displayforward();
        printf("\n");
        
        // Display pairs
        for (i=0;i<num_pairs;i++) {
            printf("(%d %d),(%d %d)\n",pairs[i][0],pairs[i][1],pairs[i][2],pairs[i][3]);
        }
        
        // Get current total distance
        best_distance = getDistance();
        
        // Find distance from switching each node pair
        i = 0;
        index_best = -1;
        for (i=0;i<num_pairs;i++) {
            printf("FLIPPING: %d -> %d, ",pairs[i][1],pairs[i][3]);
            flipNodes(pairs[i][1], pairs[i][3]);
            displayforward();
            printf(", ");
            
            // Check if total distance is improved
            current_distance = getDistance();
            
            if (current_distance >= best_distance) {
                printf("(d = %f)\n",current_distance);
                flipNodes(pairs[i][3], pairs[i][1]);
            } else {
                best_distance = current_distance;
                printf("(d = %f) (BEST)\n",best_distance);
                if (trial == 0) first_best_distance = best_distance;
                index_best = i;
                flipNodes(pairs[i][3], pairs[i][1]);
            }
        }
        
        // If no better distance was found, end search
        if (index_best==-1) {
            break;
        } else {
            flipNodes(pairs[index_best][1], pairs[index_best][3]);
        }
        
        // Update overall best distance
        if (best_distance<first_best_distance) {
            first_best_distance = best_distance;
        }
        trial++;
        printf("----------\n");
    }
    
    printf("...Done. Best distance = %f: (",first_best_distance);
    displayforward();
    printf(")\n");
}

template<class T>
void doublylinkedlist<T>::starOpt(int K,int LISTSIZE)
{
    int num_nodes = countNodes();
    int i,k,m,n = 0;
    int temp[2*K];
    int flag,trial,matches;
    float current_distance = 0; // starting distance
    float best_distance = getDistance(); // starting best distance
    vector < vector <int> > pairs;
    doublylinkedlist<int> tempList;
    
    displayforward();
    printf(", Distance = %f\n",best_distance);
    for (trial=0;trial<20;trial++) {
        // Find K pairs
        k = 0;
        while (k<K) {
            flag = 0;
            temp[2*k] = rand() % num_nodes;
            temp[2*k+1] = getNextIndex(temp[2*k]);
            for (i=0;i<k;i++) { // Make sure pairs have not been seen in current set
                if (temp[2*i]==temp[2*k] || temp[2*i+1]==temp[2*k] || temp[2*i]==temp[2*k+1] || temp[2*i+1]==temp[2*k+1]) flag = 1;
            }
            if (flag==0) k++;
        }
        
        // Make sure pairs have not been seen in history
        flag = 0;
        for (m=0;m<n;m++) {
            matches = 0;
            for (k=0;k<K;k++) {
                for (i=0;i<K;i++) {
                    matches += temp[2*i]==pairs[m][2*k];
                }
            }
            if (matches==K) {
                flag = 1;
                break;
            }
        }
        
        // Test pairs
        if (flag==0) {
            printf("BEST LIST = ");
            displayforward();
            tempList = copyList(0,LISTSIZE-1);
            tempList.flipKItems(K,temp);
            current_distance = tempList.getDistance();
            printf("CURRENT LIST = ");
            tempList.displayforward();
            printf(", Distance = %f",current_distance);
            
            // Keep flip and store in history if distance is improved, else revert
            if (current_distance<best_distance) {
                best_distance = current_distance; // Update best distance
                pairs.resize(n+1);
                for (m=0;m<=n;m++) {
                    pairs.at(m).resize(2*K);
                }
                for (k=0;k<K;k++) {
                    pairs[n][2*k] = temp[2*k];
                    pairs[n][2*k+1] = temp[2*k+1];
                }
                n++;
                flipKItems(K,temp);
                printf(" (BEST)");
            }
            printf("\n");
        }
    }
}



template<class T>
void doublylinkedlist<T>::flipKItems(int K, int *temp) {
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



template<class T>
void doublylinkedlist<T>::findKPairs(int K, vector<vector<int> > *pairs) {
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


template<class T>
void doublylinkedlist<T>::findPairs(vector<vector<int> > *pairs) {
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

template<class T>
int doublylinkedlist<T>::getNextIndex(int n)
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

template<class T>
int doublylinkedlist<T>::countNodes() {
    node *p = head;
    p = p->next;
    int num_nodes = 1;
    while(p!=head) {
        num_nodes++;
        p = p->next;
    }
    return num_nodes;
}


template<class T>
float doublylinkedlist<T>::getDistance() {
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
//This function flips four edges with center n1, n3
// for instance, flipping 3,7 from 0-1-2-...-7
// then we end up changing edges (2,3)(3,4) and edges (6,7)(7,0)
//Done by TR.
template<class T>
void doublylinkedlist<T>::flipNodes(int n1, int n3) {
    
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

//Insert item after key
template<class T>
void doublylinkedlist<T>::insertAfter(T item, T x, T y, T key) {
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
template<class T>
void doublylinkedlist<T>::displayforward() {
	node *p=head;
	while(1) {
        printf("%d ",p->data);
        //printf("%d (%d,%d), ",p->data,p->x,p->y);
        p=p->next;
        if (p==head) break;
	}
}

//displaying list nodes in reverse direction
template<class T>
void doublylinkedlist<T>::displaybackward() {
	node *p=head->prev;
	while(1) {
		printf("%d ",p->data);
		p=p->prev;
        if (p==head->prev) break;
	}
}

//fliping two edges according to 2opt
//This function takes in (pair12,pair22) and (pair21,pair22) and flips them
template<class T>
void doublylinkedlist<T>:: flipTwoItems(const T pair11, const T pair12, const T pair21, const T pair22){
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


template<class T>
void doublylinkedlist<T>::showDistance() {
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
template<class T>
void doublylinkedlist<T>::deleteNode(T i){ 
	node *p;
	p = head;
	if(p->data != i){
		p= p->next;
		while(p!=head){
			if(p-> data !=i)
				p = p->next;
			else break;			
		}
	}
	p -> prev -> next = p -> next;
	p -> next -> prev = p -> prev;
    if(p == head)
        head = head-> next;
}


//returns an integer array of the array
template<class T>
void doublylinkedlist<T>::extractIndices(int* arr){
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
template<class T>
void doublylinkedlist<T>::appendList(doublylinkedlist<int> l2){
    node* end = head-> prev;
    int num = l2.countNodes();
    cout<<num<<" == number of l2"<<endl;
    for (int i=0; i< num; i++) {
        int arr[3];
        l2.getContentAt(i, arr);
        int lastKey = head-> prev->data;
        cout<<"Got data " << arr[0]<<","<<arr[1]<<","<<arr[2]<<endl;
        insertAfter(arr[0], arr[1], arr[2], lastKey);
    }
}

//get the content of key, x and y at the nodeCount (0 to numNodes-1)
template<class T>
void doublylinkedlist<T>::getContentAt(T nodeCount, int arr[3]){
    node* p = head;
    for (int i = 0; i < nodeCount; i++) {
        p= p-> next;
    }
    arr[0] = p-> data;
    arr[1] = p-> x;
    arr[2] = p-> y;
}



#endif
