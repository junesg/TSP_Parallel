# Traveling Salesman Problem MPI Implementation

## Method I have coded and developed:


### MST (create initial solution from a file)
#### Usage: MST.cpp, graph.cpp, doublylinked.cpp, and DisjointSets.cpp
	*	calling files
	*	call DLlFromMST(string file name)
	*	the MST method reads from the file, performs minimum spanning tree, and uses a depth first search to reach a final solution
	
	
### GA 
#### Usage: TSP_LK.cpp, and doublylinkedlist.cpp
	*	Genetic algorithm
	*	Currently this is still an int-main implemented file 
	*	the final method should take in a doubly linked list and return the best or a batch of best solutions to the problem
	
	
### Two-Opt


### Ray opt 


### KL

### Star opt




###MemeticMPI.cpp (overarching mpi file regulating the parallel behavior)
*** compiler
mpic++ MemeticMPI.cpp doublylinked.cpp MST.cpp GA.cpp graph.cpp TSP_LK.cpp DisjointSets.cpp
*** running
mpirun -np <number of processes> <program name and arguments>



TODOS:
when the np is smaller than 6, we need to find a way to loop through all methods before starting to let the parallel system migrate



MPI: to prevent usage of too much RAM we don't want to send too many messages before posting receives. we also don't want the slowest processor to get jammed. 


ONE OF THE PROBLEM Of current:
 you can't declare double arr[some var representing int] = {...}
 
 
  mpic++ MemeticMPI.cpp doublylinked.cpp MST.cpp GA.cpp graph.cpp TSP_LK.cpp DisjointSets.cpp -g
 

To Compile:
mpic++ MemeticMPI.cpp doublylinked.cpp MST.cpp GA.cpp graph.cpp TSP_LK.cpp DisjointSets.cpp -g
mpirun -np 6 ./a.out

g++ memeticSerial.cpp doublylinked.cpp MST.cpp GA.cpp graph.cpp TSP_LK.cpp DisjointSets.cpp -g
./a.out


When running in flux (MPI parallel super computer machine): 
1. change file name 
2. change iteration size
3. In the serial: give a fixed amount of time
4. adjust GA parameters, depth for LK and for two-opt , ray-opt and star opt