
#include "mainSubstrata2.h"


/**************************** MAIN R ******************************************/
/* read in input and output from R                                            */
/**************************** MAIN R ******************************************/
void R_substrata2 (
      double * x,      /* the column major listing of items      k*n*d  */
      int * kInt,      /* the length of an individual item       1      */
      int * dInt,      /* number of commodities                  1      */
      int * NInt,      /* number of individual items             1      */
      int * iterInt,   /* number of interations                  1      */
      int * IInt,      /* destructive: init/final stratification n      */
      double * Q,      /* cost value of each commodity           d      */
      double * adminDataDbl,   /* Administrative data double             ANDbl  */
      int * adminDataInt,      /* Administrative data integer            ANInt  */
      int * ANIntInt,  /* number of double admin recordds        1      */
      int * ANDblInt,  /* number of integer admin recors         1      */
      int * dup,       /* UNTESTED/UNUSED dup vector             n      */
      double * costChangeDbl, /* (costChange,T) per item          iter * 2 */
      int * auxFunctionIterInt, /* how often to run the auxiliary function 1 */
      int * costChangeSizeInt      /* size of a cost change element           1 */
    ) { 
  
  /* note that we get ints from R and we want to work with size_t, which are not the same for most systems */ 
  size_t i;

  /* make size_t transfers */
  size_t k = (size_t) *kInt;    /* vector size for each distance matrix */
  size_t dN = (size_t) *dInt;   /* number of distance matricies */
  size_t N = (size_t) *NInt;    /* number of vectors (obs) */
  size_t iter = (size_t) * iterInt; /* number of iterations fo the algorithm */
  size_t auxFunctionIter = (size_t) *auxFunctionIterInt; /* run the auxiliary update program this often */
  size_t ANInt = (size_t) * ANIntInt;
  size_t ANDbl = (size_t) * ANDblInt;
  size_t costChangeSize = (size_t) *costChangeSizeInt;

  /***** 0.2 ASSIGN OUR MEMORY *****/

  size_t * I = malloc(sizeof(size_t) * N );
  size_t * J = malloc(sizeof(size_t) * N );

  /* copy over some data */
  for( i=0; i < N; i++) I[i] = (size_t) IInt[i];



  double * R = malloc( sizeof(double) * dN );

  /***** 0.3 MAKE INITIAL ASSIGNMENT *****/
  

  /* created to take care of the convert issue */ 
  size_t * size = malloc(sizeof(size_t) * N );

  for( i=0; i < N; i++) size[i] =  (size_t) adminDataInt[i];

  /* create distance matrix */
  double * D = createDistMatrix( x, k, dN, N, squaredEuclidianDist, size ); 
 
  free(size);
  
  
  /* app specific administrative data      */
  void * A = substrata2_packSubstrata(I,D,adminDataInt, adminDataDbl, dN,N, ANInt, ANDbl); 


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
    auxFunctionIter,
    costChangeSize,
    substrata2_init,          // init function                                 //
    substrata2_randomState,
    substrata2_costChange,    // cost change function                          //
    substrata2_update,        // update function                               //
    substrata2_cool,          // cooling schedule                              // 
    substrata2_diag           // diagnostic printout function                  //
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
  substrata2_deleteSubstrata( A, dN, N );

  return;
}


