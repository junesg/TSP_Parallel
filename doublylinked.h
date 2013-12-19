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
    node(int data1, float x1, float y1): data(data1),x(x1), y(y1){}
};

class doublylinkedlist
{
public:
    node *head;
    node *end;
    doublylinkedlist() //constructor
    {
        node* head = NULL;
        node* end = NULL;
    }
    
   	void createList(int*,int*,int*,int);
    ~doublylinkedlist();
	void insertAfter(int,int, int, int);
	void displayforward();
	void displaybackward();
    int countNodes();
    int getNextIndex(int);
    float getDistance();
	void flipTwoItems(const int,const int);
	//void showDistance();
	void deleteNode(int);
	void rearrangeList(int start);
    void extractIndices(int*);
    void appendList(doublylinkedlist);
    void getContentAt(int, int arr[3]);
    void flipNodes(int,int);
    void flipKItems(int,int*);
	bool compareList(doublylinkedlist);
};

//copyList is the only function that dynamically allocates
//memory
doublylinkedlist* copyList(doublylinkedlist* , int start, int end);

#endif
