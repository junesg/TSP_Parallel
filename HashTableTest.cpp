#include <stdio.h>
#include <stdlib.h>#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "HashTable.hpp"
using namespace std;

int main() {

	int TABLE_SIZE = 6;

	HashMap* ABC= new HashMap(TABLE_SIZE);
	std::vector<int> thisVector;
	for(int i = 0; i< 10; i++)
	 	thisVector.push_back(i);
	ABC->put(0 , &thisVector);
	std::vector<int> thatVector;
	for(int i = 10; i< 20; i++)
	 	thatVector.push_back(i);
	 ABC->put(0 , &thatVector);
		
	cout<< "roundnumber is "<< ABC->getNumberOfRound(0)<<endl;
	vector<int> *avect = ABC->get(0);
	cout<<"The string at 0 position is "<<endl;
	for(int i = 0; i< avect->size(); i++)
	{
		cout<<avect->at(i)<<",";
	
	}
	cout<<endl;
	return 0;



}
#include <vector>
#include <iostream>
#include "HashTable.hpp"
using namespace std;

int main() {

	int TABLE_SIZE = 6;
	
	HashMap* ABC= new HashMap(TABLE_SIZE);
	std::vector<int> thisVector;
	for(int i = 0; i< 10; i++)
	 	thisVector.push_back(i);
	ABC->put(0 , &thisVector);
	std::vector<int> thatVector;
	for(int i = 10; i< 20; i++)
	 	thatVector.push_back(i);
	 ABC->put(0 , &thatVector);
		
	cout<< "roundnumber is "<< ABC->getNumberOfRound(0)<<endl;
	vector<int> *avect = ABC->get(0);
	cout<<"The string at 0 position is "<<endl;
	for(int i = 0; i< avect->size(); i++)
	{
		cout<<avect->at(i)<<",";
	
	}
	cout<<endl;
	return 0;



}