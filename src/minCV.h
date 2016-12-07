#ifndef HEADER_MINCV
#define HEADER_MINCV


/*** headers  ***/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "R.h"
#include "Rmath.h"
#include "dist.h"
#include "sa.h"
#include "cv.h"
#include "alloc.h"

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
  size_t Hi;                //index to share second randomly selected var 
  size_t Hj;                //index to share second randomly selected var 
  size_t H;                 //Number of strata 
  size_t R;                 //the number of observations of a PSU                    
  size_t J;                 //the number of domains 
  size_t K;                 //the number of variables                    
  size_t * Nh;              //Population size by strata 
  size_t * candidate_Nh;    //Population size by strata 
  double * nh;              //sample size
  double * candidate_nh;    //candidate sample size
  double *** var;           //contribution tensor 
  double *** candidate_var; //variance matrix for candidate 
  double *** mu;            //mean matrix for R
  double *** candidate_mu;  //mean matrix for andidate 
  double * T;               //target CV 
  double * locationAdj;
  double * scaleAdj;
  double *** locationAdj_mu;
  double *** scaleAdj_mu;
  double *** candidate_locationAdj_mu;
  double *** candidate_scaleAdj_mu;
  double * x;               //A link to the Data Set 
  double **  Total;         //Total size 
  double totalProb;         //total probability of all strata
  double * totalStrataProb; //probability of selecting a strata
  double temp;              //temp denominator for cooling function
  double * prob;            //vector of sampling weights 
  double * probMatrix;      //matrix of sampling weights , one weight for each stratum 
  double * penalty;           //penalty coefficient for objective function 
  double p;                 //exponent for penalty  
  double * cv;              // cv
  double * candidate_cv;              // cv
  size_t iterSampleSize;
  size_t preserveSatisfied; //    1 - do not go above a prior met constraint, 
                            //    0 - allow moving above a prior met constraint 
  size_t *** domain;         //domain MDA
  size_t fpc;              // 1 use a finite population correction factor on CV calculation
                           // 0 do not use a finite population correction factor on CV calculation
} 
minCV_adminStruct;

typedef minCV_adminStruct * minCV_adminStructPtr;


/* packSubstrata is a function that build all required data sets
 * and converts the output to the adminSubstrataStructPtr to
 * a void pointer to pass to sa. 
 */
void * minCV_pack( 
  double * x,               // A link to the Data Set 
  size_t * I,               // assignments
  size_t N,                 // Number of strata 
  size_t K,                 // Number of strata 
  size_t H,                 // Number of strata 
  size_t R,                 // The number of observations of a PSU                    
  size_t J,                 // The number of targeted characteristics                    
  double * nh,              // sample size
  double * T,               // target CV 
  double * locationAdj,
  double * scaleAdj,
  double temp,              // temp denominator for cooling function
  double * prob,            // vector of sampling weights 
  double * probMatrix,      // matrix of sampling weights , one weight for each stratum 
  double * penalty,          // penalty coefficients for objective function 
  double p,                 // exponent for penalty  
  size_t iterSampleSize,    // number of iterations to optimize the sample size
  size_t preserveSatisfied, //    1 - do not go above a prior met constraint, 
                            //    0 - allow moving above a prior met constraint 
  int * domain,             // domain HxKxJ data structure that will be turned into an MDA
  int * fpc                 // 1 - use fpc 
                            // 0 - do not use fpc
); 


/* clean up for the substrata administrative data */
void minCV_delete( minCV_adminStructPtr  a, size_t dN, size_t N );

    
/* cooling schedule                                                           */
double minCV_cool ( 
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
void minCV_diag( 
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
    ); 

/* init function */
void minCV_init (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
            ); 

/* get strata move, conditioned on a selected PSU */
size_t minCV_getMoveStrata( 
    size_t i, 
    size_t * I, 
    double * prob, 
    double * probMatrix,
    size_t N,
    size_t H
    );

/* select a PSU */
size_t minCV_getIndex( double * prob, double totalProbability, double * totalStrataProbability ); 

/* select a random state, psu and strata to move to */
size_t minCV_randomState (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * IFin,         /* new state                               */
            double * QFin,         /* new cost                                */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t K,              /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
    ); 


/* function to select an index */
size_t minCV_getIndex( double * prob, double totalProbability, double * totalStrataProbability ); 


// randomly select strata
size_t minCV_getMoveStrata( 
    size_t i, 
    size_t * I, 
    double * prob, 
    double * probMatrix,
    size_t N,
    size_t H
    );


/* This is a function to check the change in substrata cost                   */
double minCV_costChange (
            size_t i,              /* new Index                               */
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * IFin,         /* new state                               */
            double * QFin,         /* new cost                                */
            double * D,            /* distance matrix     (NOT USED)          */
            void * A,              /* administrative data                     */
            size_t K,              /* number of variables                     */
            size_t N               /* number of PSUs                          */
    );


/* update function */
void minCV_update (
            size_t accept, /* 1 if accepted 0 if not, 2 handles additional functions */
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            size_t * IFin, /* new state, destructive */
            double * QFin, /* new cost, destructive */
            double * D, /* distance matrix */
            void * A,   /* administrative data */
            size_t K,  /* number of distance matricies */
            size_t N,   /* number of elements within a state */
            double * costChange /* cost Change */
            ); 


/* cooling schedule                                                           */
double minCV_cool ( 
            size_t iter,   /* iteration */
            double diff, /* cost change */ 
            size_t i,   /* new Index */ 
            size_t * I, /* current state */ 
            double * Q, /* current cost */
            double * D, /* distance matrix */
            void * A, /* administrative data */ 
            size_t dN,   /* number of distance matricies */
            size_t N    /* number of elements within a state */
            ); 

#endif
