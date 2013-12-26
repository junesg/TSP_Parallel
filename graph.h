#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>


using namespace std;


typedef struct Vertices //the vertices in the graph
{
    int index;
    int xPos;
    int yPos;
    bool visited;
    vector <int> neighbors;
}Vertices;

typedef struct Edges //the edges as given by the MST
{
    int fromIndex;
    int toIndex;
    //double length;
    bool revised;
}Edges;


class Graph
{
public:
    ~Graph();
    Graph(){
        std::vector<Vertices*> allVert;
        std::vector<Edges*> allEdge;
    }
    
    void createGraph(std::vector< pair<int, int> > initialsolution, std::vector<std::pair<int,int> > coordinates);
    bool isConnected(int,int);
    void printGraphContent();
    std::vector<Vertices*> allVert;
    std::vector<Edges*> allEdge;
    Vertices* findVertex(int ind);
    void DFS(int, vector<int>*);
};


struct node{
    int info;
    struct node *next;
};


class stack{
    struct node *top;
public:
    stack(); //constructor of stack
    void push(int); //push a node index
    int pop();
    bool isEmpty();
    void display();
};

stack::stack(){
    top = NULL;
}

//push a new node onto the stack
void stack::push(int data){
    node *p;
    p = new node;
    p->info = data;
    p->next = NULL;
    if(top!=NULL){
        p->next = top;
    }
    top = p;
}

//take a node out of the stack
int stack::pop(){
    struct node *temp;
    int value;
    if(top==NULL){
        cout<<"\nThe stack is Empty"<<endl;
    }else{
        temp = top;
        top = top->next;
        value = temp->info;
        delete temp;
    }
    return value;
}

bool stack::isEmpty(){
    return (top == NULL);
}

void stack::display(){
    struct node *p = top;
    if(top==NULL){
        cout<<"\nNothing to Display\n";
    }else{
        cout<<"\nThe contents of Stack\n";
        while(p!=NULL){
            cout<<p->info<<endl;
            p = p->next;
        }
    }
}


