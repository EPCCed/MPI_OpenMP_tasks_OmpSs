// #undef _USE_EXTRAE_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "TaskQueue.h"
#ifdef _USE_EXTRAE_
#include "ExtraeConstants.h"
#endif

/*************************************************************************************/

// Definition of the task queue, receiving information required for the functions
// * sort is used to order the nodes in the routine EnQueue
// * seek is used to decide if a node is selected in the routine DeQueue
// * data includes information required for the previous functions
void CreateTaskQueue (ptr_TaskQueue tsk_queue, int num_nodes, void *data, 
											SortType_func sort, Seek_func seek) {
	// Parameters validation
	if ((tsk_queue == NULL) || (num_nodes < 1)) {
		printf ("Incorrect parameters in CreateTaskQueue\n"); PrintTrace (); exit (-1);
	}
	// To work correctly, the size of the queue has to be greater than one
	if (num_nodes <= 1) num_nodes = 2;
	// Malloc memory to avoid overwriting (when queuing) the first num_nodes positions of q
	CreateInts (&(tsk_queue->q), num_nodes);
	// Initialize the contents of the local vector to -1 and 
	// the rest of components to correspondent values
	InitInts (tsk_queue->q, num_nodes, -1, 0);
	tsk_queue->output = tsk_queue->input  = 0;
	tsk_queue->size = num_nodes;
	tsk_queue->reset_value = -1;
#ifdef _MULTI_THREADS_
	// Create the lock to control the multithread computation
 	omp_init_lock( &(tsk_queue->lock) );
#endif
#ifdef _BLK_COND_VARS_
	// Create the required structures to manage condition variables
	pthread_mutex_init (&(tsk_queue->mut), NULL);
	pthread_cond_init (&(tsk_queue->cond), NULL);
 	tsk_queue->root = 1;
//  #pragma omp flush(tsk_queue->root)
#endif
	// Initialization of the sort and the seek structures
	tsk_queue->data = data;
	tsk_queue->sort = sort;
	tsk_queue->seek = seek;
}

// Liberate the structures included in the queue tsk_queue.
// The routine returns the data, which have to be liberated 
// by the calling routine.
void *RemoveTaskQueue (ptr_TaskQueue tsk_queue) {
	void *data = NULL;

	if (tsk_queue != NULL) {
		// If a non null pointer has arrived to the routine
#ifdef _BLK_COND_VARS_
		// Liberation of the structures to manage condition variables
		pthread_cond_destroy (&(tsk_queue->cond));
		pthread_mutex_destroy (&(tsk_queue->mut));
#endif
#ifdef _MULTI_THREADS_
		// Liberation of the lock for multithread computation
 		omp_destroy_lock( &(tsk_queue->lock) );
#endif
		// Liberation of the local vector on the queue is implemented
		RemoveInts (&(tsk_queue->q));
		data = tsk_queue->data;
	}

	// Return the data to be liberated for the calling routine
	return data;
}

// Print all the nodes included in the queue tsk_queue
void PrintQueue (TaskQueue tsk_queue) {
	int i = tsk_queue.output, n = tsk_queue.size;
	int end = tsk_queue.input, *pi = tsk_queue.q;

	// Parameters validation
	if ((pi == NULL) || (i < 0) || (end < 0) || 
			(n <= 1) || (i >= n) || (end >= n)) {
		printf ("Incorrect parameters in PrintQueue\n"); PrintTrace (); exit (-1);
	} else {
		// Printing the data included in the queue
		while (i != end) {
			printf ("%d ", pi[i++]); if (i == n) i = 0;
		}
		printf ("\n");
	}
}

// Print all the nodes included in the queue tsk_queue, where
// * f1 is the number of positions occupied by each node
// * f2 is the number of significant positions occupied by each node
void PrintFQueue (TaskQueue tsk_queue, int f1, int f2) {
	int i = tsk_queue.output, n = tsk_queue.size;
	int end = tsk_queue.input, *pi = tsk_queue.q;
	char formato[10];

	// Parameters validation
	if ((pi == NULL) || (i < 0) || (end < 0) || 
			(n <= 1) || (i >= n) || (end >= n)) {
		printf ("Incorrect parameters in PrintFQueue\n"); PrintTrace (); exit (-1);
	} else {
		// Define the format used to print the data
		if (f1 <= 0) 
			sprintf (formato, "%%d ");
		else if (f2 <= 0)
			sprintf (formato, "%%%dd", f1);
		else
			sprintf (formato, "%%%d.%dd", f1, f2);
		// Printing the data included in the queue
		while (i != end) {
			printf (formato, pi[i++]); if (i == n) i = 0;
		}
		printf ("\n");
	}
}

// Print all the nodes included in the queue tsk_queue
// on the file defined by f
void FPrintQueue (FILE * f, TaskQueue tsk_queue) {
	int i = tsk_queue.output, n = tsk_queue.size;
	int end = tsk_queue.input, *pi = tsk_queue.q;

	// Parameters validation
	if ((f == NULL) || (pi == NULL) || (i < 0) || (end < 0) || 
			(n <= 1) || (i >= n) || (end >= n)) {
		printf ("Incorrect parameters in FPrintQueue\n"); PrintTrace (); exit (-1);
	} else {
		// Printing on f the data included in the queue
		while (i != end) {
			fprintf (f, "%d ", pi[i++]); if (i == n) i = 0;
		}
		fprintf (f, "\n");
	}
}

// Print all the nodes included in the queue tsk_queue, where
// * f1 is the number of positions occupied by each node
// * f2 is the number of significant positions occupied by each node
// on the file defined by f
void FPrintFQueue (FILE *f, TaskQueue tsk_queue, int f1, int f2) {
	int i = tsk_queue.output, n = tsk_queue.size;
	int end = tsk_queue.input, *pi = tsk_queue.q;
	char formato[10];

	// Parameters validation
	if ((f == NULL) || (pi == NULL) || (i < 0) || (end < 0) || 
			(n <= 1) || (i >= n) || (end >= n)) {
		printf ("Incorrect parameters in FPrintFQueue\n"); PrintTrace (); exit (-1);
	} else {
		// Define the format used to print the data
		if (f1 <= 0) 
			sprintf (formato, "%%d ");
		else if (f2 <= 0)
			sprintf (formato, "%%%dd", f1);
		else
			sprintf (formato, "%%%d.%dd", f1, f2);
		// Printing the data included in the queue
		while (i != end) {
			fprintf (f, formato, pi[i++]); if (i == n) i = 0;
		}
		fprintf (f, "\n");
	}
}

// This routine reads the contents of the file, whose name is in filename,
// and these values are put on the queue tsk_queue
void ReadQueue (char *filename, ptr_TaskQueue tsk_queue) {
	int beg = tsk_queue->output, n = tsk_queue->size, j = 0;
	int i = tsk_queue->input, *pi = tsk_queue->q;
	int dimTsk = 0, *tasks = NULL;

	// Parameters validation
	if ((filename == NULL) || (pi == NULL) || (beg < 0) || (i < 0) || 
			(n <= 1) || (beg >= n) || (i >= n)) {
		printf ("Incorrect parameters in ReadQueue\n"); PrintTrace (); exit (-1);
	} else {
		dimTsk = ReadInts (filename, &tasks);
		if (dimTsk == -1) {
			printf ("%s does not exist!!\n", filename); exit (-1);
		} else if (dimTsk >= n) {
			printf ("The size of %s is incorrect (%d,%d)!!\n", 
								filename, dimTsk, NumNodesQueue (*tsk_queue)); 
			exit (-1);
		}
		for (j=0; j<dimTsk; j++) {
			pi[i++] = tasks[j]; if (i == n) i = 0;
			if (i == beg) {
				printf ("The queue has been filled by the file %s!!\n", filename);
				exit (-1);
			}
		}
		RemoveInts (&tasks);
		tsk_queue->input = i;
	}
}

// This routine reads the contents of the file, whose name is in filename,
// and these values changes the current values of the queue tsk_queue
void ReadOnQueue (char *filename, TaskQueue tsk_queue) {
	int i = tsk_queue.output, n = tsk_queue.size, j = 0;
	int end = tsk_queue.input, *pi = tsk_queue.q;
	int dimTsk = 0, *tasks = NULL;

	// Parameters validation
	if ((filename == NULL) || (pi == NULL) || (i < 0) || (end < 0) || 
			(n <= 1) || (i >= n) || (end >= n)) {
		printf ("Incorrect parameters in ReadOnQueue\n"); PrintTrace (); exit (-1);
	} else {
		dimTsk = ReadInts (filename, &tasks);
		if (dimTsk == -1) {
			printf ("%s does not exist!!\n", filename); exit (-1);
		} else if (dimTsk != NumNodesQueue (tsk_queue)) {
			printf ("The size of %s is incorrect (%d,%d)!!\n", 
								filename, dimTsk, NumNodesQueue (tsk_queue)); 
			exit (-1);
		}
		for (j=0; j<dimTsk; j++) {
			pi[i++] = tasks[j]; if (i == n) i = 0;
		}
		RemoveInts (&tasks);
	}
}

// Write all the nodes included in the queue tsk_queue
// on the file whose name is filename
void WriteQueue (char *filename, TaskQueue tsk_queue) {
	FILE *f = NULL;

	f = fopen (filename, "w");
	if (f == NULL) {
		printf ("Error opening %s\n", filename); exit (-1);
	} else {
		printf ("Print %s containing a task queue\n", filename);
		FPrintQueue (f, tsk_queue);
		fclose (f);
	}
}

// Write all the nodes included in the queue tsk_queue, where
// * f1 is the number of positions occupied by each node
// * f2 is the number of significant positions occupied by each node
// on the file whose name is filename
void WriteFQueue (char *filename, TaskQueue tsk_queue, int f1, int f2) {
	FILE *f = NULL;

	f = fopen (filename, "w");
	if (f == NULL) {
		printf ("Error opening %s\n", filename); exit (-1);
	} else {
		printf ("Print %s containing a task queue\n", filename);
		FPrintFQueue (f, tsk_queue, f1, f2);
		fclose (f);
	}
}

/*************************************************************************************/

// Return 1 if the queue tsk_queue is empty
int EmptyQueue (TaskQueue tsk_queue) {
	// Parameter validation
	if ((tsk_queue.output < 0) || (tsk_queue.input < 0)) {
		printf ("Incorrect parameters in EmptyQueue\n"); PrintTrace (); exit (-1);
	} 
	return (tsk_queue.output ==  tsk_queue.input);
}

// Return the number of nodes in the queue tsk_queue
int NumNodesQueue (TaskQueue tsk_queue) {
	// Parameter validation
	if ((tsk_queue.output < 0) || (tsk_queue.input < 0) || (tsk_queue.size <= 1)) {
		printf ("Incorrect parameters in NumNodesQueue\n"); PrintTrace (); exit (-1);
	} 
	return ((tsk_queue.input>=tsk_queue.output)?0:tsk_queue.size)+
																					tsk_queue.input-tsk_queue.output;
}

// Set the value of the reset_value from the input
void SetResetQueue (ptr_TaskQueue tsk_queue) {
	// Parameter validation
	if ((tsk_queue == NULL) || (tsk_queue->input < 0)) {
		printf ("Incorrect parameters in SetResetQueue\n"); PrintTrace (); exit (-1);
	}
#ifdef _MULTI_THREADS_
	// Close the lock, useful for multithread codes
	omp_set_lock( &(tsk_queue->lock) );
#endif
	// Initialize the reset_value
	tsk_queue->reset_value = tsk_queue->input;
#ifdef _BLK_COND_VARS_
	// Initialize the root
 	tsk_queue->root = 1;
//  #pragma omp flush(tsk_queue->root)
#endif
#ifdef _MULTI_THREADS_
	// Open the lock, useful for multithread codes
	omp_unset_lock( &(tsk_queue->lock) );
#endif
}

// Set the value of output to 0, while the value of the input will be 
// * reset_value, if value is negative.
// * value, if it is nonnegative.
void ResetQueue (ptr_TaskQueue tsk_queue, int value) {
	// Parameters validation
	if ((tsk_queue == NULL) || ((value < 0) && (tsk_queue->reset_value < 0))) {
		printf ("Incorrect parameters in SetResetQueue\n"); PrintTrace (); exit (-1);
	}
#ifdef _MULTI_THREADS_
	// Close the lock, useful for multithread codes
	omp_set_lock( &(tsk_queue->lock) );
#endif
	// Initialize the input and output
 	tsk_queue->output = 0 ;
	if (value < 0) 
		tsk_queue->input = tsk_queue->reset_value ;
	else
		tsk_queue->input = value ;
#ifdef _BLK_COND_VARS_
	// Initialize the root
 	tsk_queue->root = 1;
//  #pragma omp flush(tsk_queue->root)
#endif
#ifdef _MULTI_THREADS_
	// Open the lock, useful for multithread codes
	omp_unset_lock( &(tsk_queue->lock) );
#endif
}

// The threads can be stopped trying to DeQueue a node.
// This routine wakes all the threads
void ClearQueue (ptr_TaskQueue tsk_queue) {
#ifdef _USE_EXTRAE_
	// Initialization of the event related to ClearQueue
	Extrae_event (961422222, EV_CLQUEUE);
#endif
	// Parameter validation
	if (tsk_queue == NULL) {
		printf ("Incorrect parameters in ClearQueue\n"); PrintTrace (); exit (-1);
	}
#if _BLK_COND_VARS_
	// Only useful for energy-saving techniques
//  #pragma omp flush(tsk_queue->root)
  if (tsk_queue->root) {
    tsk_queue->root = 0;
//    #pragma omp flush(tsk_queue->root)
	#ifdef VERBOSE                     
		printf("Waking up all threads Thread: %d\n", thread_num);
	#endif
		pthread_mutex_lock (&(tsk_queue->mut));
		pthread_cond_broadcast (&(tsk_queue->cond));
		pthread_mutex_unlock (&(tsk_queue->mut));
	}
#endif
#ifdef _USE_EXTRAE_
	// Finalization of the event related to ClearQueue
	Extrae_event (961422222, EV_NULL);
#endif
}

/*************************************************************************************/

// Insert the node tsk in the queue tsk_queue
void EnQueue (ptr_TaskQueue tsk_queue, int tsk) {
	int i, j;

#ifdef _USE_EXTRAE_
	// Initialization of the event related to EnQueue
	Extrae_event (961422222, EV_ENQUEUE);
#endif
	// Parameters validation
	if ((tsk_queue == NULL) || (tsk_queue->q  == NULL) || (tsk_queue->output < 0) || 
			(tsk_queue->input < 0) || (tsk_queue->size <= 1) || 
			((NumNodesQueue (*tsk_queue)+1) >= tsk_queue->size)) {
		printf ("Incorrect parameters in EnQueue (%d,%d,%d)\n", 
						tsk_queue->output, tsk_queue->input, tsk_queue->size); 
		PrintTrace (); exit (-1);
	}
#ifdef _MULTI_THREADS_
	// Close the lock, useful for multithread codes
 	omp_set_lock( &(tsk_queue->lock) );
#endif
	// Read the actual value of input, and update its value
	i = tsk_queue->input;
	#pragma omp atomic
		tsk_queue->input++; 
	if (tsk_queue->input == tsk_queue->size) tsk_queue->input = 0;
	// Change the order of the nodes, if the parameter sort exists
	if (tsk_queue->sort != NULL) { 
		j = ((i == 0)?tsk_queue->size:i)-1;
		// Look for the correct position in the queue, using the function sort 
		while ((i != tsk_queue->output) && (tsk_queue->sort(tsk_queue->data, tsk, tsk_queue->q[j]))) {
			tsk_queue->q[i] = tsk_queue->q[j];
			i = j; j = ((i == 0)?tsk_queue->size:i)-1;
		}
	}
	// Insert the node tsk in its final position
  tsk_queue->q[i] = tsk;
#ifdef _MULTI_THREADS_
	// Open the lock, useful for multithread codes
   omp_unset_lock( &(tsk_queue->lock) );	  
#endif
#ifdef _BLK_COND_VARS_
	// Close the mutex lock, useful to manage the condition variables
	pthread_mutex_lock (&(tsk_queue->mut));
	// Usually, only one thread has to be woken up
	if (tsk_queue->seek == (SeekOwner)) //	#ifdef _OWNR_FILE_
		// In this case, a specific thread has to be woken up, but it is impossible
		// to specify a thread using cond_signal. Therefore, the broadccast is used.
		pthread_cond_broadcast (&(tsk_queue->cond));
	else
		// Otherwise, only a thread is woken up
		pthread_cond_signal (&(tsk_queue->cond));
	// Open the mutex lock, useful to manage the condition variables
 	pthread_mutex_unlock (&(tsk_queue->mut));
#endif
#ifdef _USE_EXTRAE_
	// Finalization of the event related to EnQueue
	Extrae_event (961422222, EV_NULL);
#endif
}

// Select a node of the queue tsk_queue.
// The parameter tid is only used if the function seek exists
int DeQueue (ptr_TaskQueue tsk_queue, int tid) {
	int tsk = -1;
	int i, j;

#ifdef _USE_EXTRAE_
	// Initialization of the event related to DeQueue
	Extrae_event (961422222, EV_DEQUEUE);
#endif
	// Parameters validation
	if ((tsk_queue == NULL) || (tsk_queue->q  == NULL) || (tsk_queue->output < 0) || 
			(tsk_queue->input < 0) || (tsk_queue->size <= 1) || 
			((tsk_queue->seek != NULL) && (tid < 0))) {
		printf ("Incorrect parameters in DeQueue\n"); PrintTrace (); exit (-1);
	}
	// If the queue tsk_queue is not empty
	if (tsk_queue->input != tsk_queue->output)  {
#ifdef _MULTI_THREADS_
		// Close the lock, useful for multithread codes
		omp_set_lock( &(tsk_queue->lock) );
#endif
		// Verify again if the queue tsk_queue is not empty
		if (tsk_queue->input != tsk_queue->output) {
			// If the function seek doesn't exist
			if (tsk_queue->seek == NULL) {
				// Return the first element of the queue
				// and modify the actual value of output
				tsk = tsk_queue->q[tsk_queue->output++];
				if (tsk_queue->output == tsk_queue->size) tsk_queue->output = 0;
			} else { // If the function seek exists
				// Look for the first node, which satisfies the condition defined in the function seek.
				// Go from the output to the input
				i = tsk_queue->output; tsk = tsk_queue->q[i];
				while ((i != tsk_queue->input) && (tsk_queue->seek (tsk_queue->data, tsk, tid) == 0)) {
					i++; if (i == tsk_queue->size) i = 0;
					tsk = tsk_queue->q[i];
				}
				// If no node satisfies the condition
				if (i == tsk_queue->input) tsk = -1;
				else {
					// The node in the position i is tsk, and the queue has to be compressed
					j = ((i == 0)?tsk_queue->size:i)-1;
					while (i != tsk_queue->output) {
						tsk_queue->q[i] = tsk_queue->q[j];
						i = j; j = ((i == 0)?tsk_queue->size:i)-1;
					}
					// Update the value of output
					#pragma omp atomic
						tsk_queue->output++;
					if (tsk_queue->output == tsk_queue->size) tsk_queue->output = 0;
				}
			}
		}
#ifdef _MULTI_THREADS_
		// Open the lock, useful for multithread codes
		omp_unset_lock( &(tsk_queue->lock) );
#endif
	}
#ifdef _USE_EXTRAE_
	// Finalization of the event related to DeQueue
	Extrae_event (961422222, EV_NULL);
#endif

	// Return the selected node.
	return tsk;
}

// Select a node of the queue tsk_queue.
// The parameter tid is only used if the function seek exists
// The execution of the thread is stopped if the queue tsk_queue is empty.
int DeQueueBl (ptr_TaskQueue tsk_queue, int tid) {
	int tsk = -1, cond = 1;
	int i, j;

#ifdef _USE_EXTRAE_
	// Initialization of the event related to DeQueueBl
	Extrae_event (961422222, EV_DEQUEBL);
#endif
	// Parameters validation
	if ((tsk_queue == NULL) || (tsk_queue->q  == NULL) || (tsk_queue->output < 0) || 
			(tsk_queue->input < 0) || (tsk_queue->size <= 1) || 
			((tsk_queue->seek != NULL) && (tid < 0))) {
		printf ("Incorrect parameters in DeQueueBl\n"); PrintTrace (); exit (-1);
	}
	// If the queue tsk_queue is not empty
	if (tsk_queue->input != tsk_queue->output)  {
#ifdef _MULTI_THREADS_
		// Close the lock, useful for multithread codes
		omp_set_lock( &(tsk_queue->lock) );
#endif
		// Verify again if the queue tsk_queue is not empty
		if (tsk_queue->input != tsk_queue->output) {
			// If the function seek doesn't exist
			if (tsk_queue->seek == NULL) {
				// Return the first element of the element
				// and modify the actual value of output
				tsk = tsk_queue->q[tsk_queue->output++];
				if (tsk_queue->output == tsk_queue->size) tsk_queue->output = 0;
			} else { // If the function seek exists
				// Look for the first node, which satisfy the condition defined in the function seek.
				// Go from the output to the input
				i = tsk_queue->output; tsk = tsk_queue->q[i];
				while ((i != tsk_queue->input) && (tsk_queue->seek (tsk_queue->data, tsk, tid) == 0)) {
					i++; if (i == tsk_queue->size) i = 0;
					tsk = tsk_queue->q[i];
				}
				// If no node satisfies the condition
				if (i == tsk_queue->input) tsk = -1;
				else {
					// The node in the position i is tsk, and the queue has to be compressed
					j = ((i == 0)?tsk_queue->size:i)-1;
					while (i != tsk_queue->output) {
						tsk_queue->q[i] = tsk_queue->q[j];
						i = j; j = ((i == 0)?tsk_queue->size:i)-1;
					}
					// Update the value of output
					#pragma omp atomic
						tsk_queue->output++;
					if (tsk_queue->output == tsk_queue->size) tsk_queue->output = 0;
				}
			}
		}
		// Verify if a node has been found
		cond = (tsk == -1);
#ifdef _MULTI_THREADS_
		// Open the lock, useful for multithread codes
		omp_unset_lock( &(tsk_queue->lock) );
#endif
	}
#ifdef _BLK_COND_VARS_
	// If no node satisfies the condition the execution of the thread is stopped
	if (cond) {
		pthread_mutex_lock (&(tsk_queue->mut));
		pthread_cond_wait (&(tsk_queue->cond), &(tsk_queue->mut));
		pthread_mutex_unlock (&(tsk_queue->mut));
	}
#endif
#ifdef _USE_EXTRAE_
	// Finalization of the event related to DeQueue
	Extrae_event (961422222, EV_NULL);
#endif

	// Return the selected node.
	return tsk;
}

/*************************************************************************************/

// Seek for the node tsk in the queue tsk_queue
int SeekTaskInQueue (ptr_TaskQueue tsk_queue, int tsk) {
	int fnd = 0;
	int i;

	// Parameters validation
	if ((tsk_queue == NULL) || (tsk_queue->q  == NULL) || (tsk_queue->output < 0) || 
			(tsk_queue->input < 0) || (tsk_queue->size <= 1)) {
		printf ("Incorrect parameters in SeekTaskInQueue\n"); PrintTrace (); exit (-1);
	}
	if (tsk_queue->input != tsk_queue->output)  {
#ifdef _MULTI_THREADS_
		// Close the lock, useful for multithread codes
		omp_set_lock( &(tsk_queue->lock) );
#endif
		i = tsk_queue->output;
		// Begin from the output to the input
		while ((i != tsk_queue->input) && (tsk != tsk_queue->q[i])) {
			i++; if (i == tsk_queue->size) i = 0;
		}
		fnd = (tsk == tsk_queue->q[i]);
#ifdef _MULTI_THREADS_
		// Open the lock, useful for multithread codes
		omp_unset_lock( &(tsk_queue->lock) );
#endif
	}

	// Return if the node tsk exists in the task queue
	return fnd;
}

/*************************************************************************************/

