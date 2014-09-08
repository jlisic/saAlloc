#include "sa.h"


/* the simulated annealing algorithm */
/* given an initial state s_0, i = 0, and objective function Q.
 * 1. select a potential move t_{i+1}
 * 2. if the move improves the objective function then s_{i+1} = t_{i+1}
 * 3. else if the Pr( U < T( i ) ), for a random variate U then s_{i+1} = t_{i+1}
 * 4. else s_{i+1} = s_{i}
 * 5. if T(i) < T (some pre-determined stopping point), then terminate
 * 6. else goto 1
 */ 

double * sa( 
    size_t * I,              /* initial state                                 */
    double * Q,              /* initial costs                                 */
    size_t * J,              /* final state                                   */
    double * R,              /* final costs                                   */
    double * D,              /* distance matrix                               */
    double * costChange,     /* the cost change                               */
    void * A,                /* administrative data                           */
    size_t dN,               /* number of distance matricies                  */
    size_t N,                /* number of elements within a state             */
    size_t m,                /* max number of iterations                      */
    initFunctionPtr initSAFunction,  /* initialize function                           */
    randomStateFunctionPtr randomStateSAFunction,
    costChangeFunctionPtr costChangeFunction,
                             /* cost change function                          */
    updateFunctionPtr updateFunction,/* update function                               */
    coolFunctionPtr coolFunction,    /* cooling schedule                              */ 
    diagFunctionPtr diagFunction     /* diagnostic function                           */ 
    ) 
{

  size_t newI = 0;           /* potential index                               */
  size_t i;                  /* iteration index                               */
  double newCostChange=0;    /* change in total cost                          */
  double T;                  /* value of cooling schedule                     */
  double U;                  /* U(0,1) random variable                        */


  /* take care of R random number generator */
#ifdef CLI 
  set_seed(1,1); 
  costChange = (double *) malloc(sizeof(double) * m * 3);
#endif

#ifndef CLI
  GetRNGstate();
#endif
  i = 0;

  /* initialize the algorithm */
  initSAFunction( I, Q, D, A, dN, N);
 

  #ifndef CLI
    Rprintf("\n****************** init ****************\n");
  #endif 
  #ifdef CLI
    printf("\n****************** init ****************\n");
  #endif 
  diagFunction(i,I,Q,A,dN,N);
 
  while( i < m ) 
  {
    if( i % 100 == 0) Rprintf("iter: %d\n",(int) i);
    #ifndef CLI 
    if( i % 1000 == 0) Rprintf("%d ",(int) i);
    #endif 
    #ifdef CLI
    if( i % 1000 == 0) printf("%d ",(int) i);
    #endif

    // 1. select a potential move t_{i+1} 
    newI = randomStateSAFunction( I, Q, J, R, D, A, dN, N);


    /* if there are no possible movements we terminate */
    if( newI > N ) return( costChange ); 


    newCostChange = costChangeFunction( newI, I, Q, J, R, D, A, dN, N);
    costChange[i*3] = newCostChange;
    costChange[i*3+1] = -1;
    costChange[i*3+2] = 0;


//#ifdef DEBUG3 
    Rprintf("\ndelta=%f,", newCostChange);
//#endif

#ifdef DEBUG2 
    Rprintf("\n****************** post cost change  ****************\n");
    diagFunction(i,I,Q,A,dN,N);
#endif

    /* 2. if the move improves the objective function then s_{i+1} = t_{i+1} */
    if( newCostChange <= 0 )
    {
      costChange[i*3+2] = 1;
#ifdef DEBUG
      printf("Accepted\t");
#endif

      updateFunction(1, newI, I, Q, J, R, D, A, dN, N);   
#ifdef DEBUG2
      Rprintf("\n****************** post accepted update  ****************\n");
      diagFunction(i,I,Q,A,dN,N);
#endif
    }
    else {

      T = coolFunction( i, newCostChange, newI, I, Q, D, A, dN, N);   
#ifndef CLI
      U = runif(0.0,1.0); 
#endif
#ifdef CLI
      U = (double) rand() / (double) RAND_MAX; 
#endif
      costChange[i*3+1] = U;
      
#ifdef DEBUG
      printf("Non-optimal, T = %f, U = %f:  ", T, U);
#endif     
      /* 3. else if the Pr( U < T( i ) ), for a random variate U then s_{i+1} = t_{i+1} */
      if( U < T ) 
      {
#ifdef DEBUG
        printf("Accepted\t");
#endif
        costChange[i*3+2] = 1;
        updateFunction(1, newI, I, Q, J, R, D, A, dN, N);   
#ifdef DEBUG
        printf("\n****************** post accepted update non-opt  ****************\n");
        diagFunction(i,I,Q,A,dN,N);
#endif
      } 
      else
      {
#ifdef DEBUG
        printf("Not Accepted\t");
#endif
      }  
      /* 4. else s_{i+1} = s_{i} */
      updateFunction(0, newI, I, Q, J, R, D, A, dN, N);   
#ifdef DEBUG2
      Rprintf("\n****************** post accepted update reject  ****************\n");
      diagFunction(i,I,Q,A,dN,N);
#endif

    }

//#ifdef DEBUG  
  printf("\n************************* i = %d **************************\n",(int) i);
  for( size_t ii=0; ii < N; ii++) printf("\n I[%d]: %d", (int) ii, (int) I[ii] ); 
//#endif


    /* 6. ommitted, it is implemented as the while loop */ 
    /* 6. else goto 1 */
#ifdef DEBUG2
  Rprintf("\n");
#endif


    i++; 
  } 
  
#ifndef DEBUG
  diagFunction(i,I,Q,A,dN,N);
#endif

#ifndef CLI
  PutRNGstate();
#endif


  return( costChange );

}


/************************** SUPPORT FUNCTIONS *********************************/


/* function that frees a multidimensional array of size n */
void freeMDA( void ** X, size_t t ) {
  size_t i;


  for( i= 0; i < t; i++) free(X[i]);
  free(X);

  return;
}


/* function that prints out a matrix of type double */
void printMatrixFullDbl(double ** X , size_t row, size_t col ) {
  size_t i,j;

  for(i = 0; i < row; i++) {
    #ifndef CLI 
    Rprintf("%d:\t",(int) i);
    #endif 
    #ifdef CLI
    printf("%d:\t",(int) i);
    #endif
    for(j = 0; j < col; j++) {
      #ifndef CLI 
      Rprintf("%0.4f\t",X[i][j]);
      #endif 
      #ifdef CLI
      printf("%0.4f\t",X[i][j]);
      #endif
    }
    #ifndef CLI 
    Rprintf("\n");
    #endif 
    #ifdef CLI
    printf("\n");
    #endif
  }

}


/* function that prints out a matrix of type integer */
void printMatrixFullSize_t(size_t ** X , size_t row, size_t col ) {
  size_t i,j;

  for(i = 0; i < row; i++) {
    #ifndef CLI 
    Rprintf("%d:\t",(int) i);
    #endif 
    #ifdef CLI
    printf("%d:\t",(int) i);
    #endif
    for(j = 0; j < col; j++) {
      #ifndef CLI 
      Rprintf("%d\t",(int) X[i][j]);
      #endif 
      #ifdef CLI
      printf("%d\t",(int) X[i][j]);
      #endif
    }
    #ifndef CLI 
    Rprintf("\n");
    #endif 
    #ifdef CLI
    printf("\n");
    #endif
  }

}


