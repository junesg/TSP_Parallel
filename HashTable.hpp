#ifndef HASH_TABLE_HPP
#define HASH_TABLE_HPP

#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

/* LinkedHashEntry is the class to be stored in the Hash Table */
/* The values in the storage is method sequence and method iterations history */
class LinkedHashEntry {
private:
    int key; //the key is the rank of the processor 0 to #proc-1
    std::vector<double>* value; //this value includes the the following: 
    					// n: the count of terms in the method sequence
    					// the array of method sequence
    					// m: the count of terms in the method iteration sequence
    					// the array of iteration sequence
    LinkedHashEntry *next;
public:
    LinkedHashEntry(int key, std::vector<double>* value) {
        this->key = key;
        this->value = value;
        this->next = NULL;
    }
    
    LinkedHashEntry(int key) {
        this->key = key;
        this->value = NULL;
        this->next = NULL;
    }
    //function to delete this entry
    ~LinkedHashEntry() {
    	delete value;   
    }
    
    int getKey() {
        return key;
    }
    
    vector<double>* getValue() {
        return value;
    }
    
    void setValue(vector<double>* value) {
        this->value = value;
    }
    
    LinkedHashEntry *getNext() {
        return next;
    }
    
    void setNext(LinkedHashEntry *next) {
        this->next = next;
    }
};











/* HashMap is the implementation of HashMap class used to store past commands */
class HashMap {
private:
    LinkedHashEntry **table;
    int TABLE_SIZE;
public:
    HashMap(int TABLE_SIZE) {
    	this->TABLE_SIZE = TABLE_SIZE;
        table = new LinkedHashEntry*[TABLE_SIZE];
        for (int i = 0; i < TABLE_SIZE; i++)
            table[i] = NULL;
    }
    
    vector<double>* get(int key) { //key starts from zero to tablesize-1
    	//first create a null pointer to be returned
    	vector<double>* nullVect1;
        //int hash = (key % TABLE_SIZE);
        if (table[key] == NULL) {
            return NULL;
        }
        else {
            LinkedHashEntry *entry = table[key];
            while (entry != NULL && entry->getKey() != key)
                entry = entry->getNext();
            if (entry == NULL)
                return NULL;
            else
                return entry->getValue();
        }
    }
    
    int getNumberOfRound(int key) {
        vector<double>* nullVect;
        int rounds = 0;
	    if (table[key] == NULL) {
            return rounds;
        }
        else {
            LinkedHashEntry *entry = table[key];
            while (entry != NULL && entry->getKey() != key) {
                rounds ++;
            }
            rounds ++;
            if (entry == NULL)
                return rounds;
            else
                return rounds;
        }
	}
    
    void put(int key, vector<double>* value) {
        if (table[key] == NULL)
            table[key] = new LinkedHashEntry(key, value);
        else {
            LinkedHashEntry *entry = table[key];
            while (entry->getNext() != NULL)
                entry = entry->getNext();
            if (entry->getKey() == key)
                entry->setValue(value);
            else
                entry->setNext(new LinkedHashEntry(key, value));
        }
    }
    
    void remove(int key) {
        if (table[key] != NULL) {
            LinkedHashEntry *prevEntry = NULL;
            LinkedHashEntry *entry = table[key];
            while (entry->getNext() != NULL && entry->getKey() != key) {
                prevEntry = entry;
                entry = entry->getNext();
            }
            if (entry->getKey() == key) {
                if (prevEntry == NULL) {
                    LinkedHashEntry *nextEntry = entry->getNext();
                    delete entry;
                    table[key] = nextEntry;
                } else {
                    LinkedHashEntry *next = entry->getNext();
                    delete entry;
                    prevEntry->setNext(next);
                }
            }
        }
    }
//    
//    void printMap() {
//        printf("Printing HashMap hsitory:\n");
//        for (int key=0; key< TABLE_SIZE; key++) {
//
//            if (table[key] != NULL) {
//                printf("\nTable Entry %d\n",key);
//                printf("==========================\n",key);
//                vector<double>* values = this->get(key);
//                while (!values->empty() && values != NULL) {
//                        for (int i = 0; i < values->size(); i++) {
//                            printf("%f, ", values->at(i));
//                        }
//                        printf("\n ----------------\n");
//                    this->remove(key);
//                    values = this->get(key);
//                }
//            }
//            printf("=====================\n");
//        }
//    }
    
    ~HashMap() {
        for (int i = 0; i < TABLE_SIZE; i++)
            if (table[i] != NULL) {
                LinkedHashEntry *prevEntry = NULL;
                LinkedHashEntry *entry = table[i];
                while (entry != NULL) {
                    prevEntry = entry;
                    entry = entry->getNext();
                    delete prevEntry;
                }
            }
        delete[] table;
    }
};

#endif

//end of file