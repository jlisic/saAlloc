#ifndef HEADER_SA
#define HEADER_SA
/*** headers  ***/


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "R.h"
#include "Rmath.h"
#include "dist.h"

/*** functional defines ***/

#define SA_GETINDEX(x) (size_t) ( runif(0.0,1.0) * (double) (x) ) 

#define SA_GETINDEX_DOUBLE(x) ( runif(0.0,1.0) * (double) (x) ) 



/*** typedefs  ***/
/* initialize sa function                                                     */
typedef void ( * initFunctionPtr )(
            size_t *, /* current state */
            double *,  /* current cost */
            double *, /* distance matrix */
            void *,   /* administrative data */
            size_t,   /* number of distance matricies */
            size_t    /* number of elements within a state */
            );

/* randomState function                                                        */
typedef size_t ( * randomStateFunctionPtr )(
            size_t *, /* current state */
            double *,  /* current cost */
            size_t *, /* new state */
            double *,  /* new cost */
            double *, /* distance matrix */
            void *,   /* administrative data */
            size_t,   /* number of distance matricies */
            size_t    /* number of elements within a state */
            );

/* costChange function                                                        */
typedef double ( * costChangeFunctionPtr )(
            size_t,   /* new Index */
            size_t *, /* current state */
            double *,  /* current cost */
            size_t *, /* new state */
            double *,  /* new cost */
            double *, /* distance matrix */
            void *,   /* administrative data */
            size_t,   /* number of distance matricies */
            size_t    /* number of elements within a state */
            );

/* update function                                                            */
typedef void ( * updateFunctionPtr )(
            size_t ,  /* 1 if accepted 0 if not */
            size_t,   /* new Index */
            size_t *, /* current state */
            double *, /* current cost */
            size_t *, /* new state, destructive */
            double *, /* new cost, destructive */
            double *, /* distance matrix */
            void *,   /* administrative data */
            size_t,   /* number of distance matricies */
            size_t    /* number of elements within a state */
            );
    
/* cooling schedule                                                           */
typedef double ( * coolFunctionPtr )( 
            size_t,   /* iteration */ 
            double,   /* costChange */ 
            size_t,   /* new Index */ 
            size_t *, /* current state */ 
            double *, /* current cost */
            double *, /* distance matrix */
            void *, /* administrative data */ 
            size_t,   /* number of distance matricies */
            size_t    /* number of elements within a state */
            );                      

/* diagnostic function                                                        */
typedef void ( * diagFunctionPtr )(
            size_t,   /* new Index */
            size_t *, /* current state */
            double *, /* current cost */
            void *,   /* administrative data */
            size_t,   /* number of distance matricies */
            size_t    /* number of elements within a state */
            );

/******************************************************************************/
/************************** FUNCTION PROTOTYPES *******************************/
/******************************************************************************/

/* the simulated annealing algorithm */
/* given an initial state s_0, i = 0, and objective function Q.
 * 1. select a potential move t_{i+1}
 * 2. if the move improves the objective function then s_{i+1} = t_{i+1}
 * 3. else if the Pr( U < T( i ) ), for a random variate U then s_{i+1} = t_{i+1}
 * 4. else s_{i+1} = s_{i}
 * 5. if T(i) < T (some pre-determined stopping point), then terminate
 * 6. else goto 1
 */ 
double *  sa( 
    size_t * I,              /* initial state                                 */
    double * Q,              /* initial costs                                 */
    size_t * J,              /* final state                                   */
    double * R,              /* final costs                                   */
    double * D,              /* distance matrix                               */
    double * costChange,     /* cost change                                   */
    void * A,                /* administrative data                           */
    size_t dN,               /* number of distance matricies                  */
    size_t N,                /* number of elements within a state             */
    size_t m,                /* max number of iterations                      */
    size_t auxFunctionIter,  /* how often to run the auxiliar update function */
    initFunctionPtr initSAFunction,  /* initialize function                           */
    randomStateFunctionPtr randomStateSAFunction, /* random state function */
    costChangeFunctionPtr costChangeFunction,
                             /* cost change function                          */
    updateFunctionPtr updateFunction,/* update function                               */
    coolFunctionPtr coolFunction,    /* cooling schedule                              */ 
    diagFunctionPtr diagFunction     /* diagnostic function                           */ 
    ); 


/************************** SUPPORT FUNCTIONS *********************************/

/* function that frees a multidimensional array of size n */
void freeMDA( void ** X, size_t depth ); 

/* function that prints out a matrix of type double */
void printMatrixFullDbl(double ** X , size_t row, size_t col ); 

/* function that prints out a matrix of type integer */
void printMatrixFullSize_t(size_t ** X , size_t row, size_t col ); 

#endif
