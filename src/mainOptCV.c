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

#include "optCV.h"

/* this needs to be updated for other platforms if not using gcc */
/* and also for running on 32 bit archs */
#define F_SIZET "%zu"

/**************************** MAIN R ******************************************/
/* read in input and output from R                                            */
/**************************** MAIN R ******************************************/
void R_optCV (
      double * x,      /* the column major listing of items      k*n*d  */
      int * kInt,      /* the length of an individual item       1      */
      int * dInt,      /* number of commodities                  1      */
      int * NInt,      /* number of individual items             1      */
      int * iterInt,   /* number of interations                  1      */
      int * IInt,      /* destructive: init/final allocation     n      */
      double * Q,      /* cost value of each commodity           d      */
      double * adminDataDbl,   /* Administrative data double             ANDbl  */
      int * adminDataInt,      /* Administrative data integer            ANInt  */
      int * ANIntInt,  /* number of double admin recordds        1      */
      int * ANDblInt,  /* number of integer admin recors         1      */
      int * dup,       /* UNTESTED/UNUSED dup vector             n      */
      double * costChangeDbl /* (costChange,T) per item          iter * 2 */
    ) { 
  
  /* note that we get ints from R and we want to work with size_t, which are not the same for most systems */ 
  size_t N;                          /* number of items                       */
  size_t dN;                         /* sets of number of items               */ 
  size_t k;                          /* the length of each vector in x        */ 
  size_t iter;                       /* the number of iterations              */ 
  size_t * I;                        /* initial assignment                    */
  size_t * J;                        /* final assignment                      */
  double * R;                        /* objective function value returned     */
  double * D;                        /* distance matrix                       */
 
  void * A;                          /* app specific administrative data      */
  size_t ANInt;
  size_t ANDbl;

  size_t i;

  size_t * size;

  /* make size_t transfers */
  k = (size_t) *kInt;
  dN = (size_t) *dInt;
  N = (size_t) *NInt;
  iter = (size_t) * iterInt;
  ANInt = (size_t) * ANIntInt;
  ANDbl = (size_t) * ANDblInt;
  /***** 0.2 ASSIGN OUR MEMORY *****/

  I          = (size_t * ) malloc(sizeof(size_t) * N );

  /* copy over some data */
  for( i=0; i < N; i++) 
    I[i] = (size_t) IInt[i];


  J          = (size_t * ) malloc(sizeof(size_t) * N );

  R          = (double * ) malloc(sizeof(double) * dN * 2 );

  /***** 0.3 MAKE INITIAL ASSIGNMENT *****/
  
  /* create a distance matrix */


  /* created to take care of the convert issue */ 
  size          = (size_t * ) malloc(sizeof(size_t) * N );

  for( i=0; i < N; i++) size[i] =  (size_t) adminDataInt[i];

  D = createDistMatrix( x, k, dN, N, squaredEuclidianDist, size ); 
 
  free(size);


#ifdef DEBUG
  i = 0;
  while( i < dN) 
  {
    printDistMatrix( D, i, N );
    i++;
  }
#endif

  A = optCV_packSubstrata(I,D,x,adminDataInt, adminDataDbl, dN,N,k, ANInt, ANDbl); 
  /***** 1.0 RUN OUR ALGORITHM  *****/

   sa( 
    I,                       // initial state                                 //
    Q,                       // initial costs                                 //
    J,                       // final state                                   //
    R,                       // final costs                                   //
    D,                       // distance matrix                               //
    costChangeDbl,              // returned cost Change                          //
    A,                       // administrative data                           //
    dN,                      // number of distance matricies                  //
    N,                       // number of elements within a state             //
    iter,                       // max number of iterations                      //
    optCV_init,          // init function                                 //
    optCV_randomState,
    optCV_costChange,    // cost change function                          //
    optCV_update,        // update function                               //
    optCV_cool,          // cooling schedule                              // 
    optCV_diag           // diagnostic printout function                  //
    ); 

  /***** 2.0 WRITE RESULTS *****/

  /* copy index back */
  for( i=0; i < N; i++) 
    IInt[i] = (int) I[i];
  
  /* clean up */
  free(I);
  free(J);
  free(D);
  free(R);



  /* app specific clean up */ 
  optCV_deleteSubstrata( A, dN, N );

  return;
}


