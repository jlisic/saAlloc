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
    size_t auxFunctionIter,  /* how often to run the auxiliar update function */
    size_t costChangeSize,
    initFunctionPtr initSAFunction,  /* initialize function                   */
    randomStateFunctionPtr randomStateSAFunction,
    costChangeFunctionPtr costChangeFunction,
                             /* cost change function                          */
    updateFunctionPtr updateFunction,/* update function                       */
    coolFunctionPtr coolFunction,    /* cooling schedule                      */ 
    diagFunctionPtr diagFunction     /* diagnostic function                   */ 
    )  {

  size_t newI = 0;           /* potential index                               */
  size_t i;                  /* iteration index                               */
  double newCostChange=0;    /* change in total cost                          */
  double T;                  /* value of cooling schedule                     */
  double U;                  /* U(0,1) random variable                        */

/* take care of R random number generator */
  GetRNGstate();
  i = 0;

/* initialize the algorithm */
//Rprintf("\n****************** init ****************\n");
  initSAFunction( I, Q, D, A, dN, N);
 
  while( i < m ) {
//Rprintf("\n****************** start = %d ****************\n", (int) i );
//diagFunction(i,I,Q,A,dN,N);


    // print status every 1000th 
    if( m > 1 ) { 
      if( m < 10000 ) {
        Rprintf("Percent Complete: %4.2f\r", (float) (100*i)/m);
      } else if( (1000 * i) % 65535 == 0 ) {
        Rprintf("Percent Complete: %4.2f\r", (float) (100*i)/m);
        //R_checkUserInterrupt();
      }
    }
    
    
// 1. select a potential move t_{i+1} 
    newI = randomStateSAFunction( I, Q, J, R, D, A, dN, N);

    // if there are no possible movements we terminate 
    if( newI > N ) {
      return( costChange ); 
    }

    costChange[i * costChangeSize + 1] = -1;
    costChange[i * costChangeSize + 2] = 0;        // accept 
    costChange[i * costChangeSize + 3] = newI + 1; // move from 
  
    newCostChange = costChangeFunction( newI, I, Q, J, R, D, A, dN, N);
//diagFunction(i,I,Q,A,dN,N);
//printf("costChangeFunction - %f\n", newCostChange);  
    costChange[i * costChangeSize] = newCostChange;

    // 2. if the move improves the objective function then s_{i+1} = t_{i+1} 
    if( newCostChange <= 0 ) {

      costChange[i * costChangeSize + 2] = 1;
      updateFunction(1, newI, I, Q, J, R, D, A, dN, N, NULL);   

    } else {

      T = coolFunction( i, newCostChange, newI, I, Q, D, A, dN, N);   
      U = runif(0.0,1.0); 
      
      costChange[i * costChangeSize + 1] = U;
      costChange[i * costChangeSize + 2] = T;
      
      // 3. else if the Pr( U < T( i ) ), for a random variate U then s_{i+1} = t_{i+1} 
      if( U < T ) {

        updateFunction(1, newI, I, Q, J, R, D, A, dN, N, NULL);   

      } 
// 4. else s_{i+1} = s_{i} 
// accept state 0:
// copy back to the candidate states the original values
      updateFunction(0, newI, I, Q, J, R, D, A, dN, N, NULL );   

    }
      
// update cost change 
// accept state 2:
// change back candidate states and copy back  
    updateFunction(2, newI, I, Q, J, R, D, A, dN, N, &(costChange[i * costChangeSize ]) );   

// 6. ommitted, it is implemented as the while loop 
// 6. else goto 1 
//

    i++; 
  } 
 
  PutRNGstate();
  
// output status   
  Rprintf("Percent Complete: %4.2f\n", (int) (i*100)/m  );


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
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < col; j++) {
      Rprintf("%0.4f\t",X[i][j]);
    }
    Rprintf("\n");
  }

}


/* function that prints out a matrix of type integer */
void printMatrixFullSize_t(size_t ** X , size_t row, size_t col ) {
  size_t i,j;

  for(i = 0; i < row; i++) {
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < col; j++) {
      Rprintf("%d\t",(int) X[i][j]);
    }
    Rprintf("\n");
  }

}


