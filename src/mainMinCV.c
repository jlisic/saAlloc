/* Jonathan Lisic */
/* standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "R.h"

/* program headers */
#include "dist.h"
//#include "file.h"
#include "sa.h"

#include "minCV.h"

/* this needs to be updated for other platforms if not using gcc */
/* and also for running on 32 bit archs */
#define F_SIZET "%zu"

/**************************** MAIN R ******************************************/
/* read in input and output from R                                            */
/*  -- old --
1      double * x,       the column major listing of items      k*n*d  
2      int * kInt,       the length of an individual item       1     
3      int * dInt,       number of commodities                  1      
4      int * NInt,       number of individual items             1      
5      int * iterInt,    number of interations                  1      
6      int * IInt,       destructive: init/final allocation     n      
7      double * Q,       cost value of each commodity           d      
8      double * adminDataDbl,   / Administrative data double             ANDbl  
9      int * adminDataInt,      / Administrative data integer            ANInt  
10      int * ANIntInt,  / number of double admin recordds        1      
11      int * ANDblInt,  / number of integer admin recors         1      
12      int * dup,       / UNTESTED/UNUSED dup vector             n      
13      double * costChangeDbl,      / (costChange,T) per item          iter * 3 
14      double * doubleSampleSize,   / output sample size                        
15      int * auxFunctionIterInt,    / how often to run the auxiliary function 1 
16      int * costChangeSizeInt      / size of a cost change element           1 
*/
/**************************** MAIN R ******************************************/

void R_minCV (
    double * x,      /* the column major listing of items              K*N*R          1*/
    int * NInt,      /* number of PSUs                                                2*/
    int * JInt,      /* number of domains                                             3*/
    int * KInt,      /* number of variables                                           4*/
    int * HInt,      /* number of strata                                              5*/
    int * RInt,      /* number of observations of each PSU                            6*/
    int * iterSampleSizeInt,   /* number of interations of the sampling algorithm  1  7*/
    int * IInt,      /* destructive: init/final allocation               n            8*/
    int * domainInt, /*                                              H*K*J            9*/
    double * prob,       /*                                            N             10*/ 
    double * probMatrix, /*                                            N*H           11*/ 
    double * target,           /*                                                    12*/
    double * locationAdjDouble, /*                                                   13*/
    double * scaleAdjDouble,    /*                                                   14*/
    double * pDouble,           /*                                                   15*/
    double * penalty,           /*                                                   16*/
    double * nh,      /* sample size by stratum                             H        17*/
    double * costChangeDbl,      /* data saved at each iteration    5 * H + K        18*/
    double * temp,               /* temperature for SA                               19*/
    int * auxFunctionIterInt,    /* how often to run the auxiliary function 1        20*/
    int * costChangeSizeInt,     /* size of a cost change element           1        21*/
    int * iterInt,               /*                                                  22*/
    int * preservedSatisfiedInt
  ) { 
  
  /* note that we get ints from R and we want to work with size_t, which are not the same for most systems */ 
  size_t i,h,j;
  
  size_t N = *NInt;
  size_t J = *JInt;
  size_t K = *KInt;
  size_t H = *HInt;
  size_t R = *RInt;

  size_t preserveSatisfied = *preservedSatisfiedInt;

  /* make size_t transfers */
  size_t iter = (size_t) *iterInt; /* number of iterations fo the algorithm */
  size_t iterSampleSize = (size_t) *iterSampleSizeInt; /* number of iterations fo the algorithm */
  size_t auxFunctionIter = (size_t) *auxFunctionIterInt; /* run the auxiliary update program this often */
  size_t costChangeSize = (size_t) *costChangeSizeInt;
  
  double * D = NULL;

  /***** 0.2 ASSIGN OUR MEMORY *****/

  size_t * I = malloc(sizeof(size_t) * N );
  size_t * IFin = malloc(sizeof(size_t) * N );
  
  /* copy over some data */
  for( i=0; i < N; i++) {
    I[i] = (size_t) IInt[i];
    IFin[i] = I[i];
  }

  /* objective function value returned     */
  double * QFin  =  calloc(1, sizeof(double));
  double * Q  =  calloc(1, sizeof(double));
  



  /***** 0.3 CREATE DATA SET *****/
  
  /* app specific administrative data      */
  Rprintf("Creating Data Set\n"); 
  void * A =minCV_pack( 
    x,               // A link to the Data Set 
    I,               // assignments
    N,               // Number of strata 
    K,               // Number of strata 
    H,               // Number of strata 
    R,               // The number of observations of a PSU                    
    J,               // The number of domains 
    nh,              // sample size
    target,          // target CV 
    locationAdjDouble,
    scaleAdjDouble,
    *temp,           // temp denominator for cooling function
    prob,            // vector of sampling weights 
    probMatrix,      // matrix of sampling weights , one weight for each stratum 
    penalty,         // penalty coefficient for objective function 
    *pDouble,               // exponent for penalty  
    iterSampleSize,  // number of iterations to optimize the sample size
    preserveSatisfied, //    1 - do not go above a prior met constraint, 
                     //    0 - allow moving above a prior met constraint 
    domainInt        // domain HxKxJ data structure that will be turned into an MDA
   );


  /* record initial values into the output */
  minCV_adminStructPtr a = (minCV_adminStructPtr) A; 
  for( h = 0; h < 6; h++) costChangeDbl[h] = -1.0; 
 
  /* add on sample size */
  for( h = 0; h < H; h++) costChangeDbl[6+h] = a->candidate_nh[h]; 
  
  
  /* add on current cv */
  for( j = 0; j < J; j++) costChangeDbl[6+H+j] = a->candidate_cv[j]; 
   
  a = NULL;


  //minCV_diag( 1, I, Q, A, K, N); 


  /***** 1.0 RUN OUR ALGORITHM  *****/
  Rprintf("Starting Run\n");

   sa( 
     I,                       // initial state
     Q,                       // initial costs
     IFin,                    // final state
     QFin,                    // final costs
     D,                       // distance matrix
     costChangeDbl + costChangeSize,           // returned cost Change
     A,                       // administrative data
     K,                       // number of variables 
     N,                       // number of elements within a state
     iter,                    // max number of iterations
     auxFunctionIter,         // how often to update mu/var
     costChangeSize,
     minCV_init,              // init function
     minCV_randomState,       // function to generate random state
     minCV_costChange,        // cost change function 
     minCV_update,            // update function 
     minCV_cool,              // cooling schedule
     minCV_diag               // diagnostic printout function  
    ); 
  /***** 2.0 WRITE RESULTS *****/

  /* copy index back */
  for( i=0; i < N; i++) IInt[i] = (int) I[i];

  /* clean up */
  free(I);
  free(IFin);
  free(QFin);

//  deleteMDA( (void *) domain, MDA_SIZE_T, 3, H, K); 

  /* app specific clean up */ 
  minCV_delete( A, K, N );

  return;
}


