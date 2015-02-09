/* this is a class of functions to efficiently calculate and re-calculate the cv */

#include "cv.h"

/* function that creates a variance MDA (matrix) [commodity][strata][obs] */
double *** cv_createMeanMatrix( 
    size_t * I, 
    size_t dN, 
    size_t N, 
    size_t kN, 
    size_t H, 
    double * x, 
    size_t * Nh
  ) {

  size_t d,h,k,i,j;

  /* allocate memory for a matrix of strata by observation for each characteristic*/ 
  double *** mu = cv_double3DMatrixCreate(dN,H,kN); 

  size_t xN = N * dN * kN;

  #pragma omp for private(d,k,j,h) 
  for(i = 0; i < xN; i++ ) { 

     d = i / (N*kN); /* commodity index */
     k = i % (kN);   /* panel index */
     j = i - d*N*kN; /* observation index */
     j = j /kN;

     //printf(" d = %d, k = %d, j = %d\n", (int) d, (int) k, (int) j);


     h = I[j];       /* stratum */

     #pragma omp atomic
     mu[d][h][k] += x[i] / (double) Nh[h];  
  }
    
  return(mu);
}



/* function that creates a variance MDA (matrix) [commodity][strata] */
double *** cv_createVarMatrix( 
    size_t * I, 
    size_t dN, 
    size_t N, 
    size_t kN, 
    size_t H, 
    double * x, 
    double *** mu, 
    size_t * Nh 
  ) {

  size_t i; /* element in data structure index */

  size_t d; /* characteristic index */
  size_t j; /* observation index */
  size_t h; /* stratum index */
  size_t k; /* panel index */

  /* allocate memory for a matrix of strata by observation for each characteristic*/ 
  double *** V = cv_double3DMatrixCreate(dN,H,kN); 

  size_t xN = N * dN * kN;
  
  #pragma omp for private(d,k,j,h) 
  for(i = 0; i < xN; i++ ) { 

     d = i / (N*kN); /* commodity index */
     k = i % (kN);   /* panel index */
     j = i - d*N*kN; /* observation index */
     j = j /kN;

     h = I[j];       /* stratum */

     #pragma omp atomic
     V[d][h][k] += 
        (x[i] - mu[d][h][k]) * (x[i] - mu[d][h][k]) / ( (double) (Nh[h] - 1) ); 
  }

  return(V);
}


/* this is a function that does an online update of the mean and variance matricies */
/* this is a destructive function */
/* only the variance and mean will be updated, not I */
void cv_updateMatrix( 
    size_t * I,     /* this is the current stratification, not the modified one */ 
    size_t dN, 
    size_t N, 
    size_t kN, 
    size_t H, 
    double * x, 
    double *** mu, 
    double *** var, 
    size_t * Nh,
    size_t moveObs,       /* this is the observation to move */
    size_t moveObsStratum /* this is the stratum to move the observation to */
  ) {


  size_t d,k;
  double deltaPriorObs;
  double deltaMoveObs;

  /* get moveObs prior staratum */

  size_t priorObsStratum = I[moveObs];

  printf(" moveFrom: %d\n", (int) priorObsStratum);
  printf(" moveTo: %d\n",   (int) moveObsStratum);

  if( priorObsStratum == moveObsStratum ) {
    return;
  }

  for( d = 0; d < dN; d++){
    for( k = 0; k < kN; k++){

      /******* update the mean **********/
      /* 
        \bar{x}_{n+1} = \bar{x}_n + (x_i -\bar{x}_n)/(n+1) 
        
        \bar{x}_{n-1} = (n \bar{x}_{n} - x_i)/(n-1)  
      */
      deltaPriorObs = x[d*kN * kN*moveObs + k] - mu[d][priorObsStratum][k]; 
      mu[d][priorObsStratum][k] -= deltaPriorObs / ( (double) Nh[priorObsStratum] - 1 ); 


      /* update the mean of the new stratum */
      deltaMoveObs = x[d*kN * kN*moveObs + k] - mu[d][moveObsStratum][k]; 
      mu[d][moveObsStratum][k] += deltaMoveObs / (double) (Nh[moveObsStratum] + 1);
      
      
      /******* update the variance **********/
      /*     1   2   3   4   5
       * 1   0 d12 d13 d14 d15
       * 2 d21   0 d23 d24 d25
       * 3 d31 d32   0 d34 d35
       * 4 d41 d42 d43   0 d45
       * 5 d51 d52 d53 d54   0
       *
       * dij = (x_i - x_j)^2
       *
       * S^2_n = ( \sum_{i=1}^n \sum_{j=1}^n dij ) / (n *(n-1))
       *
       * S^2_{n+1} = (n-1)/n S^2_n + ( \sum{i=^n} 2 * di{n+1} )/(n*(n+1)) 
       */

      var[d][moveObsStratum][k] = 
        ((double) Nh[moveObsStratum] - 1) * var[d][moveObsStratum][k] +    /* M2 */ 
        deltaMoveObs *  (x[d*kN * kN*moveObs + k] - mu[d][moveObsStratum][k]); 

      var[d][moveObsStratum][k] /= (double) Nh[moveObsStratum]; 

      var[d][priorObsStratum][k] = 
        ((double) Nh[priorObsStratum] - 1) * var[d][priorObsStratum][k] -    /* M2 */ 
        deltaPriorObs *  (x[d*kN * kN*moveObs + k] - mu[d][priorObsStratum][k]); 

      var[d][priorObsStratum][k] /= (double) Nh[priorObsStratum] - 2; 


    }
  }
  
}

/**************** Functions for 2D Matricies **************************/
/* m - first index
 * n - second index
 */

/* function that deletes a matrix of type double */
double ** cv_doubleMatrixCreate( size_t m, size_t n) {
  size_t i,j;

  double ** x = malloc(sizeof(double *) * m);

  for( i = 0; i < m; i++) x[i] = calloc( n, sizeof(double) );   

  return(x);
}


/* function that deletes a matrix of type double */
void cv_doubleMatrixDelete( double ** x, size_t m) {
  size_t i;

  for( i = 0; i < m; i++) free(x[i]);

  free(x);
  return;
}


/* function that prints out a matrix of type double */
void cv_doubleMatrixPrint(double ** X , size_t m, size_t n ) {
  size_t i,j;

  for(i = 0; i < m; i++) {
    printf("%d:\t",(int) i);
    for(j = 0; j < n; j++) {
      printf("%0.4f\t",X[i][j]);
    }
    printf("\n");
  }

  return;
}


/**************** Functions for 3D Matricies **************************/
/* m - first index
 * n - second index
 * p - third index
 */


/* function to create double 3d matrix */
double *** cv_double3DMatrixCreate( size_t m, size_t n, size_t p) {

  size_t i,j,k;

  double *** x = malloc(sizeof(double **) * m);

  for( i = 0; i < m; i++) { 
    x[i] = malloc(sizeof(double *) * n );
    for( j = 0; j < n; j++) x[i][j] = calloc( p, sizeof(double) );   
  }

  return(x);
}


/* function to create double 3d matrix */
void cv_double3DMatrixDelete( double *** x, size_t m, size_t n) {

  size_t i,j;

  for( i = 0; i < m; i++) { 
    for( j = 0; j < n; j++) free(x[i][j]);
    free(x[i]);
  }

  free(x);
  return;
}


/* function that prints out the 3dMatrix */
void cv_double3DMatrixPrint( double ***x, size_t m, size_t n, size_t p) {
  size_t i,j;

  for( i = 0; i < m; i++) { 
    printf("\nV[%d]\n",(int) i),
    cv_doubleMatrixPrint(x[i], n, p ); 
  }

  return;
}


/***************** functions for totals ********************************/

/* create an array of totals */
double ** cv_createTotalMatrix( 
    size_t D, 
    size_t N, 
    size_t K, 
    size_t H,
    double ***mu, 
    size_t * Nh
    ) {

  size_t d,h,k;
  double KDouble = (double) K;

  /* created with calloc, so zero'ed out */
  double ** muArray = cv_doubleMatrixCreate(D,K); 

  for( d = 0; d < D; d++) {
    for( k = 0; k < K; k++ ) {
      for( h = 0; h < H; h++) {
        muArray[d][k] += Nh[h] * mu[d][h][k];
      }
    }
  }

  return(muArray);
}



/* objective function */
double cv_objectiveFunction( 
    size_t D, 
    size_t N, 
    size_t K, 
    size_t H, 
    double *** var, 
    size_t * Nh,
    double * nh,
    double ** Total,
    double * Target,
    double penalty
  ) {

  double result = 0;
  double KDouble = (double) K;
  double varSum[K];
  double varSumChar;

  size_t d,h,k;

  for( d = 0; d < D; d++) {
    varSumChar = 0;
    /* for each repeated observation of the survey population calculate the cv */
    for( k = 0; k < K; k++ ) {
      varSum[k] = 0;
      for( h = 0; h < H; h++) varSum[k] +=  (double) (Nh[h] * Nh[h]) /nh[h] * var[d][h][k];
   
      /* get the mean of the cv's for a characteristic */ 
//printf("cv[%d][%d]:\t%f\n", (int) d, (int) k,  sqrt(varSum[k])/Total[d][k]);
      varSumChar += sqrt(varSum[k])/( Total[d][k] * KDouble);
    }

//printf("cv[%d]:\t%f\n", (int) d, varSumChar);

    /* apply penalty */
    if( varSum[d] > Target[d] ) {
      result += (varSumChar - Target[d]) * penalty; 
    }
    result += varSumChar;
  }

  return( result );
}




/* test function */
int main() {

  double x[]  = {
    1.0, 2.0, 3.0, 
    4.0, 4.0, 5.0, 
    5.0, 5.0, 1.0, 
    2.0, 3.0, 4.0, 

    4.0, 5.0, 5.0, 
    5.0, 1.0, 2.0, 
    3.0, 4.0, 4.0, 
    5.0, 5.0, 5.0, 

    0.0, 2.0, 3.0, 
    4.0, 5.0, 5.0, 
    5.0, 5.0, 0.0, 
    2.0, 3.0, 4.0, 
    
    4.0, 6.0, 5.0, 
    5.0, 0.0, 2.0, 
    3.0, 4.0, 4.0, 
    5.0, 7.0, 5.0 };

  size_t I[]  = {  0,   0,   0,   0,   1,   1,   1,   1 };
  size_t Nh[] = {  4,   4 };
  double nh[] = {  2,   2 };
  double target[] = { 1, 1 };

  size_t dN = 2;
  size_t N  = 8;
  size_t kN = 3;
  size_t H  = 2;


  double penalty = 0;
  double ** total;

  printf("Mean\n");
  double *** mu = cv_createMeanMatrix(I, dN, N, kN, H, x, Nh); 
  cv_double3DMatrixPrint( mu, dN, H, kN); 

  printf("Variance\n");
  double *** var = cv_createVarMatrix(I,dN, N, kN, H, x, mu, Nh); 
  cv_double3DMatrixPrint( var, dN, H, kN); 

  total =  cv_createTotalMatrix( dN, N, kN, H, mu, Nh);
  printf("Total\n");
  cv_doubleMatrixPrint(total, dN, kN ); 

  printf("obj: %f\n", 
    cv_objectiveFunction( dN, N, kN, H, var, Nh, nh, total, target, penalty )
    ); 
 

  /* perform update */
  cv_updateMatrix( I, dN, N, kN, H, x, mu, var, Nh, 0, 1 ); 
  
  printf("Updated Mean\n");
  cv_double3DMatrixPrint( mu, dN, H, kN); 

  printf("Updated Variance\n");
  cv_double3DMatrixPrint( var, dN, H, kN); 
 

  /* clean up */ 
  printf(" Deleteing Mean\n");
  cv_double3DMatrixDelete( mu, dN, H); 
  printf(" Deleteing Var\n");
  cv_double3DMatrixDelete( var, dN, H); 


  free(total);

  printf("Finished\n");
  return( 0 ) ;
}



