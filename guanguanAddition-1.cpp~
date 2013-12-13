
void showDistance();


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



//Now loop through to make sure head is always the 0_th item
p=head;
if (p->data != 0) {
    p = p-> next;
    while (p->data!=0) {
        p = p->next;
    }
    head = p;
}



node* findNode(T);
void deleteNode(T);


//find the pointer to the node with data == i
node* doublylinkedlist<T>::findNode(T i) {
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
	if(p->data != i) {
		cout<<"Error: node "<<i<<" not found!"<<endl;
		return null;
	}
	return p;
}


//delete the ith node of the list
template<class T>
void doublylinkedlist<T>::deleteNode(T i){ 
	node *p = findNode(i);
	p -> prev -> next = p -> next;
	p -> next -> prev = p -> prev;	
	free(p);
}



