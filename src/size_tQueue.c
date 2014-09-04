/* intQueue src */
#include "size_tQueue.h"


/* function to create a new Q */
queuePtr createQueue( size_t item ) {

  queuePtr Q; /* declare our new pointer */

  Q = malloc(sizeof(queue));
  //assert( Q != NULL );

  Q->head = NULL;
  Q->tail = NULL;
  Q->item = item;
  return(Q);
} 


/* function to delete a Q and set to NULL */
void deleteQueue( queuePtr q) {

while( q != NULL )
 q = removeQueue(q,q);  

return;
}


/* function to remove a node in a queue */
queuePtr removeQueue( queuePtr q, queuePtr qi) { 

  //assert( q != NULL);
  //assert( qi != NULL);


  if( qi->head != NULL)  {
    (qi->head)->tail = qi->tail;
  }
  else { 
    q = qi->tail;
  }
  if( qi->tail != NULL)  {
    (qi->tail)->head = qi->head;
  } 

  free(qi);
  qi = NULL;

  return(q);
}


/* function to add a node before qi in a queue with a particular item */
queuePtr addQueue(queuePtr q, queuePtr qi, size_t item) {
  queuePtr qNew,qs;
 
  /* create new queue */
  qNew = createQueue( item );

  /* in case there isn't anything there */
  if( q == NULL) {
    q = qNew;
    return(q);
  } 
  /* mess around with the tail node */ 
  else if( qi == NULL ) {
    /* we need to find the last pointer */
    qs = q;
    while(qs->tail != NULL) qs = qs->tail; 
    qNew->head = qs;
    qs->tail = qNew;
  }   
  /* mess around with the head node */ 
  else if( qi == q ) {
    qNew->tail = q;  /* qNew -> q */
    q->head = qNew;   
    q = qNew;
  }
  /* for everyone else */ 
  else {
    qNew->tail = qi;
    qNew->head = qi->head;
    (qi->head)->tail = qNew;
    qi->head   = qNew;
  }

  return(q); 
}


/* function to search a Q for an item, this function returns the first occurence */
queuePtr searchQueue( queuePtr q, size_t item) {
  queuePtr qi;

  //assert(q != NULL );

  qi = q;

  while( qi != NULL ){
    if( qi->item == item ) break; 
    qi = (qi->tail);
  } 

  return(qi);
}





