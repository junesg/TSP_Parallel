
/*
 * The master program will allocate tasks to the slaves
 */
static void master(void)
{
    int sizeWorld, rank;
    MPI_Status status;
    MPI_Request request;
    
    
    /*
     * Prepare initial condition of the problem
     */
    //First read in the problem
    //Give the problem a file name, the file has o be in the same folder as the code
    string filename = "testDist.txt";
    doublylinkedlist* newDLL = startingDLL(filename);
    
    
    for (rank = 1; rank < sizeWorld; ++rank) {
        
        /* Find the next item of work to do */
        work = get_next_work_item();
        
        /* Send it to each rank */
        
        MPI_Send(&work,             /* message buffer */
                 1,                 /* one data item */
                 MPI_INT,           /* data item is an integer */
                 rank,              /* destination process rank */
                 WORKTAG,           /* user chosen message tag */
                 MPI_COMM_WORLD);   /* default communicator */
    }
    
    
    
    
    /* Loop over getting new work requests until there is no more work
     to be done */
    
    work = get_next_work_item();
    while (work != NULL) {
        
        /* Receive results from a slave */
        
        MPI_Recv(&result,           /* message buffer */
                 1,                 /* one data item */
                 MPI_DOUBLE,        /* of type double real */
                 MPI_ANY_SOURCE,    /* receive from any sender */
                 MPI_ANY_TAG,       /* any type of message */
                 MPI_COMM_WORLD,    /* default communicator */
                 &status);          /* info about the received message */
        
        /* Send the slave a new work unit */
        
        MPI_Send(&work,             /* message buffer */
                 1,                 /* one data item */
                 MPI_INT,           /* data item is an integer */
                 status.MPI_SOURCE, /* to who we just received from */
                 WORKTAG,           /* user chosen message tag */
                 MPI_COMM_WORLD);   /* default communicator */
        
        /* Get the next unit of work to be done */
        
        work = get_next_work_item();
    }
    
    /* There's no more work to be done, so receive all the outstanding
     results from the slaves. */
    
    for (rank = 1; rank < ntasks; ++rank) {
        MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }
    
    /* Tell all the slaves to exit by sending an empty message with the
     DIETAG. */
    
    for (rank = 1; rank < ntasks; ++rank) {
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }
    
    
    
    
    /*
     * Generate initial Method selection and number of iterations
     */
    //Initial Sequence matrix
    int Sequence[6]={0,1,2,3,4,5};
    randomize(Sequence, 6);
    //    printArray(Sequence,6);
    int Frequency[6] = {0,0,0,0,0,0};
    Frequency[rankWorld] = 1;
    //    printArray(Frequency,6);
    
	MPI_Finalize();
	return 0;
    
}


static void
slave(void)
{
    int work; //how many unit of tasks to do
    int results; //
    MPI_Status status;
    
    while (1) {
        
        // Receive a message from the master
        //blocking receive
        MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
        
        //If the tag is to terminate, do nothing.
        if (status.MPI_TAG == TERMINATE_TAG) {
            return;
        }
        
        //otherwise do the work
        result = do_work(work);
        
        //return the result to the master
        MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}



/*
 * Functions below are helper functions for initializing the master
 */
//find the initial doublylinkedlist from the problem file
doublylinkedlist* startingDll(string filename)
{
    vector<double> edgeWeight; //edgeWeight is coupled wth the vertexPair function
    vector<std::pair<int,int> > coordinates; //later expanded in getEdgeWeight function
    vector<std::pair<int,int> > vertexPair; //later expanded in getEdgeWeight function
    getEdgeWeight(&edgeWeight, &coordinates, &vertexPair, filename);
    
    int n = coordinates.size();
    int xPos[n], yPos[n],ind[n];
    int count = 0;
    for (vector<std::pair<int,int> >::iterator it = coordinates.begin(); it != coordinates.end(); it++) {
        ind[count] = count;
        xPos[count] = (*it).first;
        yPos[count] = (*it).second;
        count ++;
    }
    //finished initializing the element of doublylinkedlist
    doublylinkedlist* newDLL = new doublylinkedlist();
    newDLL->createList(ind, xPos, yPos, n);
    return newDll;
}


// A utility function to print an array
void printArray (int arr[], int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

// A function to generate a random permutation of arr[]
void randomize ( int arr[], int n )
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    srand (time(NULL));
    
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
        
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
    }
}
//function used to swap the content of two positions
void swap (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}




/*
 * Functions below are helper functions for doing work
 */

static unit_result_t
do_work(unit_of_work_t work)
{
    /* Fill in with whatever is necessary to process the work and
     generate a result */
}



/*
 * Functions below are helper functions for processing the results
 */
static void
process_results(unit_result_t result)
{
    /* Fill in with whatever is relevant to process the results returned
     by the slave */
}




