# Traveling Salesman Problem MPI Implementation

	*Examples of how to use each section of the various algorithms are all enclosed in the source code.
	*The algorithms used in this whole program:
	
	1. Minimum Spanning Tree (kurskals algorithm)
	2. Depth First Search 
	3. Genetic Algorithm
	4. Memetic Algorithm
	5. Two-Opt
	6. LK (for TSP)
	7. Star-Opt(novel)
	8. Ray-Opt(novel)
	9. master-slave MPI
	10. Memetic algorithm with sub-heuristic trading
	
	
	*Data Structures actually employed:
	1. doubly linked list
	2. hash table
	3. disjoint sets
	4. graph
	5. lots and lots of vectors.....


## Method I have coded and developed:

### MST (create initial solution from a file)
#### Usage: MST.cpp, graph.cpp, doublylinked.cpp, and DisjointSets.cpp
	*	calling files
	*	call DLlFromMST(string file name)
	*	the MST method reads from the file, performs minimum spanning tree, and uses a depth first search to reach a final solution
	
	
### GA 
#### Usage: TSP_LK.cpp, and doublylinkedlist.cpp
	*	Genetic algorithm
	*	the final method should take in a doubly linked list and return the best or a batch of best solutions to the problem
	
	
### Two-Opt
#### Usage: randonly genearte a pair of edges and do a crossing between the four vertices


### Ray opt 
#### Usage: randonly genearte a pair of edges and do a direct swap of their precedents and followers. 



### KL
#### Usage: recursively finds the path that can be improved and improve upon the tour


### Star opt
#### Usage: choose an arbitrary number of opt and randomly select the edges to be swapped. 




###MemeticMPI.cpp (overarching mpi file regulating the parallel behavior)
*** compiler
mpic++ MemeticMPI.cpp doublylinked.cpp MST.cpp GA.cpp graph.cpp TSP_LK.cpp DisjointSets.cpp
*** running
mpirun -np <number of processes> <program name and arguments>


###Compiler: I haven't the chance to write a make file yet, so this will do:
To Compile:

####parallel
 *  To Compile: mpic++ MemeticMPI.cpp doublylinked.cpp MST.cpp GA.cpp graph.cpp TSP_LK.cpp DisjointSets.cpp -o output

 * To Run: mpirun -np 6 ./a.out

####serial
 * To Compile: g++ memeticSerial.cpp doublylinked.cpp MST.cpp GA.cpp graph.cpp TSP_LK.cpp DisjointSets.cpp -g
 * To run: ./a.out


When running in flux (MPI parallel super computer machine): 
1. change file name 
2. change iteration size
3. In the serial: give a fixed amount of time
4. adjust GA parameters, depth for LK and for two-opt , ray-opt and star opt
