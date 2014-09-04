#ifndef HEADER_OPTCV
#define HEADER_OPTCV


/*** headers  ***/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "R.h"
#include "Rmath.h"
#include "dist.h"
#include "sa.h"

#ifdef CLI
#define MATHLIB_STANDALONE /* set for stand alone R.h, R from C  */
#endif


/******************************************************************************/
/************************** FUNCTION PROTOTYPES *******************************/
/******************************************************************************/


/************************** SUBSTRATA PROBLEM FUNCTIONS ***********************/

/* define an administrative data type */
typedef struct 
{
  size_t j;     /* index to share second randomly selected var */
  size_t H;     /* Number of strata */
  size_t k;       /* the size of an element of x                   */
  size_t NhMax; /* max stratum size*/
  size_t * Nh;  /* size of strata vector */
  size_t ** L;  /* label matrix */
  double *** C; /* contribution tensor */
  double * T;   /* target variance vector */
  double * W;   /* within variance vector */
  double ** V;  /* variance matrix */
  double * x;   /* A link to the Data Set */
  double ** Total; /* Total size */
  double temp;  /* temp denominator for cooling function*/
  double * acres;   /* not used */
  double * NhAcres; /* not used */
  size_t * size;
  size_t * NhSize;
  double acreDif;
} 
optCV_adminStruct;

typedef optCV_adminStruct * optCV_adminStructPtr;


/* packSubstrata is a function that build all required data sets
 * and converts the output to the adminSubstrataStructPtr to
 * a void pointer to pass to sa. 
 */
void * optCV_packSubstrata( 
    size_t * I,     /* initial assignment                            */
    double * D,     /* distance matrix                               */
    double * x,     /* data D by N by K                              */
    int * aInt,     /* integer admin data from file                  */
    double * aDbl,  /* double admin data from file                   */
    size_t dN,      /* number of distance matricies                  */
    size_t N,       /* number of elements within a state             */
    size_t k,       /* the size of an element of x                   */
    size_t NInt,    /* the number of items in the admin integer data */
    size_t NDbl     /* the number of items in the admin integer data */
); 


/* clean up for the substrata administrative data */
void optCV_deleteSubstrata( optCV_adminStructPtr  a, size_t dN, size_t N );

/* init function */
void optCV_init (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
            ); 

/* cost change function */
size_t optCV_randomState (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * J,            /* new state, destructive                  */
            double * R,            /* new cost, destructive                   */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
            );

/* cost change function */
double optCV_costChange (
            size_t i,              /* new Index                               */
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * J,            /* new state, destructive                  */
            double * R,            /* new cost, destructive                   */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
            );

/* update from change in i function */
void optCV_update (
            size_t accept, /* 1 if accepted, 0 if not */
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            size_t * J, /* new state, destructive */
            double * R, /* new cost, destructive */
            double * D, /* distance matrix */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
            );
    
/* cooling schedule                                                           */
double optCV_cool ( 
            size_t iter,   /* iteration */
            double diff, /* cost change */ 
            size_t i,   /* new Index */ 
            size_t * I, /* current state */ 
            double * Q, /* current cost */
            double *, /* distance matrix */
            void * A, /* administrative data */ 
            size_t dN,   /* number of distance matricies */
            size_t N    /* number of elements within a state */
            );                      

/* diagnostic print                                                           */
void optCV_diag( 
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
    ); 

/* count the number of lables */
size_t optCV_labelCount( size_t * label, size_t N ); 

/* this is a function that creates and retrurns an array of sizes for the labels */
size_t * optCV_labelTotalPSUs ( size_t * label, size_t N, size_t H ); 
size_t * optCV_labelTotalSegments ( size_t * label, size_t N, size_t H, size_t *  size ); 
double * optCV_labelTotalAcres ( size_t * label, size_t N, size_t H, double * acres ); 

/* quick function to determine the maximum int in an array */
size_t optCV_arrayMaxSize_t( size_t * a, size_t n ); 

/* quick function to determine the maximum index in a double array */
size_t optCV_arrayMaxIndexDbl( double * a, size_t n ); 

  /* this matrix identifies each item with all other items that can have labels assigned */
size_t ** optCV_labelCreateMaster( size_t * label, size_t N, size_t H, size_t NhMax );

/* creates a matrix of differences between item i, and all the points in a 
 * stratum h.  Rows are items, and columns are strata.
 */ 
double *** optCV_createContribMatrix( size_t * label, size_t dN, size_t N, size_t H, double * D, size_t ** L, size_t * Nh);

/* function that creates a variance MDA (matrix) [commodity][strata] */
double ** optCV_createVarMatrix( size_t * label, size_t dN, size_t N, size_t H, double * D, size_t ** L, size_t * Nh, size_t * NhSize);
 
/* function that creates a total MDA (matrix) [commodity][strata] */
double ** optCV_createTotalMatrix( size_t * label, size_t k, size_t dN, size_t N, size_t H, double * x);

#endif
