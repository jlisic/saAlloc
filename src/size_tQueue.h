/* header file for intQueue */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* create a queue ADT */
struct queue { 
  size_t item;
  struct queue * head;
  struct queue * tail;
};

typedef struct queue queue;
typedef struct queue * queuePtr;

/* function to create a new Q */
queuePtr createQueue( size_t item );

/* function to delete a Q and set to NULL */
void deleteQueue( queuePtr q);

/* function to remove a node in a queue */
queuePtr removeQueue( queuePtr q, queuePtr qi); 

/* function to add a node in a queue with a particular item */
queuePtr addQueue(queuePtr q, queuePtr qi, size_t item);

/* function to search a Q for an item, this function returns the first occurence */
queuePtr searchQueue( queuePtr q, size_t item);

/* make a priority Q assignment */
queuePtr insertMinValueQueue( size_t i, double a, double * A, queuePtr q ); 






