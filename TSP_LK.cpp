#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "TSP_LK.h"
//#include "doublylinked.h"

using namespace std;

//debug
#define LISTSIZE 10





int main() {
    doublylinkedlist aList;
    int ind[LISTSIZE] = {0, 1, 2, 3, 4, 5,6,7,8,9};
    int x[LISTSIZE]   = {0, 0, 1, 1, 2, 2,4,5,6,7};
    int y[LISTSIZE]   = {10,0, 10,0, 10, 0,2,3,4,5};
    aList.createList(ind,x,y,LISTSIZE);
    aList.displayforward();
    cout<<endl;
    cout<<"original distance is "<<aList.getDistance()<<endl;

    aList = TSP_LK(aList,5);
    
    cout<<endl<<"final reslult"<<endl;
    aList.displayforward();
    printf("debug");

    cout<<"distance is "<<aList.getDistance()<<endl;
    aList.destroy();
    aList.~doublylinkedlist();
}



float distanceBetweenNodes(node* n1, node* n2){
    return sqrt(pow((n1->x - n2->x),2) +pow((n1->y - n2->y),2));
}



doublylinkedlist TSP_LK (doublylinkedlist thisTour, int MAXITER) {
    int iter = 0;
    int countNode = thisTour.countNodes();
    doublylinkedlist tour = copyList(thisTour, 0 , countNode-1);//created tour##1
    
    while (iter < MAXITER) {
        doublylinkedlist path = copyList(tour,0,countNode-1);//construct the path##2
        node* p = path.head;
        for (int i=0; i<iter; i++ ){
            p=p->next;
        }
        //rearrange to form a path
        path.rearrangeList(p->data);
        //now call improve path on this path
        vector<int> R;
        for (int i=0; i<countNode; i++) {
            R.push_back(0);
        }
        
        //construct the improved path ##
        doublylinkedlist tour2 = ImprovePath(path, 1, &R); //construct the path##3
        //path is also destroyed ##2
        
        if (tour2.getDistance() < tour.getDistance()) {
           tour.destroy();
           tour.~doublylinkedlist(); //destroys the tour##1
           tour = tour2; //tour now points to tour2
           iter = 0;
        }
       else {
           iter =iter+1;
           tour2.destroy();
           tour2.~doublylinkedlist(); //destorys the tour2##3
        }

      //  cout<<"TSP_LK destroy  path"<<endl;
       // path.destroy();
       // path.~doublylinkedlist(); //destroy the current path##2
       // cout<<"destroy vector"<<enl;
       // R.~vector();
    }
    cout<<"TSP_LK destroy original path"<<endl;
    thisTour.destroy();
    thisTour.~doublylinkedlist();
    return tour;
}

//doublylinkedlist that is put into the improve path "thisPath" will be destroyed;
doublylinkedlist ImprovePath(doublylinkedlist thisPath, int depth, vector<int> *R){
    doublylinkedlist path = copyList(thisPath,0,thisPath.countNodes()-1); //create a new path inside this function ##1
    cout<<"Improving tour "; path.displayforward();cout<<endl;
    
    //if there is three nodes in the path, no need to improve
    if (path.countNodes() <=3) {
        return path;
    }
    
    //if the depth is smaller than maxdepth, keeps on improving till a better path is found
    if (depth < MAXDEPTH) {
        for (node* p = path.head->next; p->next!= path.end; p = p->next) {
            if ((*R)[p->data]==0) {
                if (distanceBetweenNodes(p,p->next) > distanceBetweenNodes(p,path.end)) {
                    //if tour length is improved
                    cout<<"Tour length is improved"<<endl;
                   // cout<<"debug";
                    if (distanceBetweenNodes(p,p->next)+
                        distanceBetweenNodes(path.head,path.end) >
                        distanceBetweenNodes(p, path.end) +
                        distanceBetweenNodes(path.head, p->next)) {
                        //path.displayforward(); cout<<endl;
                        doublylinkedlist tour = copyList(path, 0, path.countNodes()-1); //copy a new tour from path ##2
                        cout<<" Original tour "; tour.displayforward(); cout<<endl;
                        tour.flipTwoItems(p->data, path.end->data);                     //flip the two items in tour
                        tour.end = tour.head -> prev;
                        cout<<"IM: destroy original path: path in first if"<<endl;
                        path.destroy();
                        path.~doublylinkedlist(); //destroy path ##1'
                        cout<<"IM: destroy original path: thisPath in first if"<<endl;
                        thisPath.destroy();
                        thisPath.~doublylinkedlist();//destroy original path
                        return tour; //return the path ##2
                    }
                    else {
                        path.flipTwoItems(p->data,path.end->data);  //just flip the edge
                        int thisData = p->data;
                        (*R)[thisData]=1;
                        cout<<"IM: destroy original path: thisPath in first else"<<endl;

                        thisPath.destroy();
                        thisPath.~doublylinkedlist();
                        return ImprovePath(path, depth+1, R); //return improvement of ##1
                    }
                }
            }
        }
    }
    
    else{
        float maxDist=0;
        node* maxNode;
        node* p;
        if (path.countNodes() <=3) {
            cout<<"No Need improve path"<<endl;
            return path;
        }
        path.end = path.head->prev;
       // printf("path starts at %d, ends at %d\n", path.head->data, path.end->data);
        
        //get the node that will end up giving the biggest gain over the end to head
        for (p = path.head->next; p->next!=path.end; p = p->next) {
             if (distanceBetweenNodes(p,p->next) - distanceBetweenNodes(p,path.end) > maxDist) {
                 maxDist =distanceBetweenNodes(p,p->next) - distanceBetweenNodes(p,path.end);
                 maxNode = p;
             }
        }
        if (maxDist>0) {
            if (distanceBetweenNodes(maxNode,maxNode->next)+
                distanceBetweenNodes(path.head,path.end) >
                distanceBetweenNodes(maxNode, path.end) +
                distanceBetweenNodes(path.head, maxNode->next)) {
                doublylinkedlist tour = copyList(path, 0, path.countNodes()-1); //create tour from path ##2
                path.end = path.head -> prev;
                tour.flipTwoItems(maxNode->data,path.end->data);
                cout<<"IM: destroy original path: path in second if"<<endl;
                path.destroy(); thisPath.destroy();
                path.~doublylinkedlist(); thisPath.~doublylinkedlist(); //destroy ##1 and the original thisPath
                return tour;
            }
        }
        else{
            thisPath.destroy(); thisPath.~doublylinkedlist(); //destroy original thisPath
            }
    }
    return path;

}

/*
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
*/


//RayOpt takes in a list, ouputs swapped two nodes (4 edges, 2 adjacent pairs)
doublylinkedlist rayOpt(doublylinkedlist P,int NUMITERATIONS) //number of iteration for two Opt)
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


doublylinkedlist starOpt(doublylinkedlist P, int K,int NUMITERATIONS)
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
                cout<<"COPY IN DEBUG: "<<num_nodes<<" nodes"<<endl;
                tempList.displayforward();cout<<endl<<endl;
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
doublylinkedlist TwoOpt(doublylinkedlist P, int NUMITERATIONS)
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
            
            tempList.flipTwoItems(temp[0], temp[2]);
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
