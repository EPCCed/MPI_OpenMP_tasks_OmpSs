#ifndef TaskType

#define TaskType 1

#ifdef _BLK_COND_VARS_
	#include <pthread.h>
#endif
#ifdef _MULTI_THREADS_
	#include "omp.h"
#endif
#include "SortFunction.h"

typedef struct TaskQueue {
  int input;                 // position from where the node is EnQueued
  int output;                // position from where the node is DeQueued
  int reset_value;           // default value of the input and output value
  int *q;                    // vector where the values are stored
  int size;                  // size of the vector q
	void * data;               // includes information required for the functions
	SortType_func sort;        // function uses to order the nodes in the routine EnQueue
	Seek_func seek;            // function uses to decide if a node is selected in the routine DeQueue
#ifdef _BLK_COND_VARS_
	pthread_cond_t cond;       // condition variable used to block the thread
	pthread_mutex_t mut;       // mutex required for the condition variable
	int root;                  // to verify if some thread arrives to some code section
#endif
#ifdef _MULTI_THREADS_
  omp_lock_t lock;           // lock uses in multithread access to the queue
#endif
} TaskQueue, *ptr_TaskQueue;

/*************************************************************************************/

// Definition of the task queue, receiving information required for the functions
// * sort is used to order the nodes in the routine EnQueue
// * seek is used to decide if a node is selected in the routine DeQueue
// * data includes information required for the previous functions
extern void CreateTaskQueue (ptr_TaskQueue tsk_queue, int num_nodes, void *data, 
											SortType_func sort, Seek_func seek);

// Liberate the structures included in the queue tsk_queue.
// The routine returns the data, which have to be liberated 
// by the calling routine.
extern void *RemoveTaskQueue (ptr_TaskQueue tsk_queue);

// Print all the nodes included in the queue tsk_queue
extern void PrintQueue (TaskQueue tsk_queue);

// Print all the nodes included in the queue tsk_queue, where
// * f1 is the number of positions occupied by each node
// * f2 is the number of significant positions occupied by each node
extern void PrintFQueue (TaskQueue tsk_queue, int f1, int f2);

// Print all the nodes included in the queue tsk_queue
// on the file defined by f
extern void FPrintQueue (FILE * f, TaskQueue tsk_queue);

// Print all the nodes included in the queue tsk_queue, where
// * f1 is the number of positions occupied by each node
// * f2 is the number of significant positions occupied by each node
// on the file defined by f
extern void FPrintFQueue (FILE *f, TaskQueue tsk_queue, int f1, int f2);

// This routine reads the contents of the file, whose name is in filename,
// and these values are put on the queue tsk_queue
extern void ReadQueue (char *filename, ptr_TaskQueue tsk_queue);

// This routine reads the contents of the file, whose name is in filename,
// and these values changes the current values of the queue tsk_queue
extern void ReadOnQueue (char *filename, TaskQueue tsk_queue);

// Write all the nodes included in the queue tsk_queue
// on the file whose name is filename
extern void WriteQueue (char *filename, TaskQueue tsk_queue);

// Write all the nodes included in the queue tsk_queue, where
// * f1 is the number of positions occupied by each node
// * f2 is the number of significant positions occupied by each node
// on the file whose name is filename
extern void WriteFQueue (char *filename, TaskQueue tsk_queue, int f1, int f2);

/*************************************************************************************/

// Return 1 if the queue tsk_queue is empty
extern int EmptyQueue (TaskQueue tsk_queue);

// Return the number of nodes in the queue tsk_queue
extern int NumNodesQueue (TaskQueue tsk_queue);

// Set the value of the reset_value from the input
extern void SetResetQueue (ptr_TaskQueue tsk_queue);

// Set the value of output to 0, while the value of the input will be 
// * reset_value, if value is negative.
// * value, if it is nonnegative.
extern void ResetQueue (ptr_TaskQueue tsk_queue, int value);

// The threads can be stopped trying to DeQueue a node.
// This routine wakes all the threads
extern void ClearQueue (ptr_TaskQueue tsk_queue);

/*************************************************************************************/

// Insert the node tsk in the queue tsk_queue
extern void EnQueue (ptr_TaskQueue tsk_queue, int tsk);

// Select a node of the queue tsk_queue.
// The parameter tid is only used if the function seek exists
extern int DeQueue (ptr_TaskQueue tsk_queue, int tid);

// Select a node of the queue tsk_queue.
// The parameter tid is only used if the function seek exists
// The execution of the thread is stopped if the queue tsk_queue is empty.
extern int DeQueueBl (ptr_TaskQueue tsk_queue, int tid);

/*************************************************************************************/

// Seek for the node tsk in the queue tsk_queue
extern int SeekTaskInQueue (ptr_TaskQueue tsk_queue, int tsk);

/*************************************************************************************/

#endif

