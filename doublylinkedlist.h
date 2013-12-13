#ifndef DLL_H
#define DLL_H

#include<iostream>

using namespace std;

template <class T>
class doublylinkedlist
{
	private:
		struct node {   //node struct, has data in it and has pointers to previous and next nodes
			T data;
			node *prev;
			node *next;
		};

		node *head;

	public:
		doublylinkedlist() //constructor
		{	
			head = NULL;
		}

	void createList(T*,int);
	void insertAfter(T,T);
	void deletenode(T);
	void displayforward();
	void displaybackward();
	void flipTwoItems(T,T,T,T);
};






template <class T>
void doublylinkedlist<T>::createList(T* x, int n)
{
	node *q;  //just a pointer
	node *p = new node;  //the first node
	p-> data = x[0];
	head = p;
	p-> next = NULL;
	p-> prev = NULL;
	for (int i=1; i<n; i++){
		q = p;
		p = p->next=new node;
		p-> data = x[i];
		p->next = NULL;
		p->prev = q;
	}
	p -> next  = head; //TSP is a circle.
	head -> prev = p;
}


//Insert item after key 
template<class T> 
void doublylinkedlist<T>::insertAfter(T item,T key) { 
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
		cout<<"key "<< key <<" not found"<<endl;
		return; 
	}

	node *p=new node; 
	p->data=item; 
	p->next=q->next; 
	p->prev=q; 
	q->next=p; 
} 

/*
//deleting item from double linked lisr 
template<class T> 
void doublylinkedlist<T>::deletenode(T item) { 
	if(head==NULL) { 
		cout<<"List is empty"<<endl; 
		return; 
	} 
	if(head->data==item) { 
		head=head->next; 
		head->prev=NULL; 
		return; 
	} 
	if(tail->data==item) { 
		tail=tail->prev; 
		tail->next=NULL; 
		return; 
	} 
	node *p=head->next; 
	while(p!=NULL) { 
		if(p->data==item)  break; 
		p=p->next; 
	} 

	if(p==NULL) { 
		cout<<item<<"not found "<<endl; 
		return; 
	} 
		
	(p->prev)->next=p->next; 
	(p->next)->prev=p->prev; 
	return; 
} 
*/


//diaplaying list elements in forward direction
template<class T> 
void doublylinkedlist<T>::displayforward() {
	node *p=head; 
	cout<<"\n Doubly linked list (Forward):"; 
	cout<<p->data;
	p = p->next;
	while(p!=head) { 
		cout<<","<<p->data; 
		p=p->next; 
	}	 
	cout<<endl;
} 

//displaying list elements in reverse direction 
template<class T> 
void doublylinkedlist<T>::displaybackward() { 
	node *p=head->prev; 
	cout<<"\n Doubly linked list (Backward)"; 
	while(p!=head) { 
		cout<<p->data<<""; 
		p=p->prev; 
	} 
} 

/*
template<class T>
void doublylinkedlist<T>:: compairLists(doublylinkedlist<T> list2){
}
*/
//fliping two edges according to 2opt
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

#endif


