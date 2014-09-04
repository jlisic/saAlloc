/* source for defined distances */

#include "dist.h"


/* function that creates the initial distance matrix */
/* example x:  1 2 3 4 5 6 7 8 9 10 11 12
 * for n = 3, k = 2, k = 2
 * v1,1 = (1 2) v2,2 = (3 4) v1,3(5 6) 
 * v2,1 = (7 8) v2,2 = (9 10) v2,3(11 12)
 * 
 * The first distance matrix is constructed from v1,* and 
 * the second is constructed from v2,*.
 */
   
  double * createDistMatrix(
    double * x, /* vector of values */ 
    size_t k,   /* length of vectors embedded in the values */ 
    size_t d,   /* number of distance matricies */
    size_t N,   /* the number of elements within a single distance matrix */ 
    double  (*dist)( double *, double*, size_t),
              /* the distance function used, it operates on vectors */ 
    size_t * size 
              /* the number of observations each item represents, this
               * is shared over all matricies */
) {

  double * D = NULL;
  double * iVector = NULL;
  double * jVector = NULL;
  size_t i, j, l, di, dOffset;
  /* run some sanity checks */
  /*
  assert(N > 1);
  assert(d > 0);
  */

  D = malloc(sizeof(double) * N*(N-1)/2 *d );
  //assert(D != NULL);

  /* here are some quick notes because what is going on down here might not be that obvious
   *
   * 1   X |                              this 
   * 2   X | X |                          way 
   * 3   X | X | X ...                    is j
   * .   . V . V .                        | 
   * .   .   .   .                        | 
   * .   .   .   .                        V
   *     0   1   2 ...
   *   ---------------->  this way is i
   *
   *  So by the above drawing we can see that we only need (N)(N-1)/2
   *
   *  The index for i < j would be N(N-1)/2 - (N-i)(N-i-1)/2 + (j-i-1) 
   *
   *  The way we handle the data for below with &(x[]) is to find the jth and ith vector
   *  from this big long string of vectors, recall that the size of each vector is just k
   *
   *  So if k=3 and we have vectors a,b,c,... we can see that there elements are easily
   *  covered by iVector and jVector
   *  
   *  0 1 2 3 4 5 6 7 8 . . .
   *  a a a b b b c c c . . .
   *
   *  For multiple distance matricies, we net d*N * (N-1)/2 be the starting index for 
   *  distance matrix d. 
   */


  l = 0; /* set initial increment to make things simple */

  for (di =0; di < d; di++) {

    dOffset = di*(N*k);

    for( i = 0; i < N-1; i++) { /* iterate over columns */
      iVector = &(x[dOffset + i*k]);
  
      for( j = i+1; j < N; j++) {  /* iterate over rows */
  
        jVector = &(x[dOffset + j*k]);

        if( size != NULL) D[l] = dist( iVector, jVector, k) * size[i] * size[j];
        else D[l] = dist( iVector, jVector, k);
       
        l++;
      }
  
    }
  }

  /* that was it, we are done */
  return(D);
}



/* function that destructively updates a dist matrix */
/* NOT UPDATED */
void  initUpdateDistMatrix(
    double * D, 
    double * E, 
    double * SS, 
    size_t d,
    size_t N, 
    size_t * size, 
    void (*objFunc) ( size_t , size_t , size_t , double * , double * , double * ,size_t, size_t , size_t * )
    /* this is an update function */
    ) 
{

  size_t i, j;
  /* run some sanity checks */
  //assert(N > 1);


  for( i = 0; i < N-1; i++) { /* iterate over columns */
    for( j = i+1; j < N; j++) {  /* iterate over rows */
      if((size[i] > 1) | (size[j] > 1)) objFunc( j,i,i,D,E,SS,d,N,size) ;
     }
  }
  /* that was it, we are done */
}



/* function that gets the Matrix (array) index of a specified index in the Matrix */
size_t getIndex( size_t i, size_t j, size_t d, size_t N) {
  size_t k;
  if ( i > j) {
    k = i;
    i = j;
    j = k;
  } 
  /* so now i < j */
  return( (d+1) * N*(N-1)/2 - (N-i)*(N-i-1)/2 + (j-i-1)); 
}



/* print out the dist matrix */
void printDistMatrix( double * D, size_t d, size_t N) {
  size_t i,j;

  #ifdef CLI
  printf("\n\tMatrix Dim:");
  printf(F_SIZET, d);
  printf("\n");

  if( D == NULL)
  {
    printf("\nNothing to print.\n");
    return;
  }

  for( i = 0; i < N-1; i++)
  {
    printf("\tcolumn ");
    printf(F_SIZET,i);
  }
  printf("\n");

  for( i = 1; i < N; i++) { 
    printf(F_SIZET,i);
    printf(":\t"); 

    for( j = 0; j < i; j++) {  
      printf("%f\t", getDist(i,j,D,d,N) );  
    }

    printf("\n"); 
  }
  #endif
 
  #ifndef CLI
  Rprintf("\n\tMatrix Dim:");
  Rprintf(F_SIZET, d);
  Rprintf("\n");

  if( D == NULL)
  {
    Rprintf("\nNothing to print.\n");
    return;
  }

  for( i = 0; i < N-1; i++)
  {
    Rprintf("\tcolumn ");
    Rprintf(F_SIZET,i);
  }
  Rprintf("\n");

  for( i = 1; i < N; i++) { 
    Rprintf(F_SIZET,i);
    Rprintf(":\t"); 

    for( j = 0; j < i; j++) {  
      Rprintf("%f\t", getDist(i,j,D,d,N) );  
    }

    Rprintf("\n"); 
  }
  #endif

  return;
}




/* get a particular value from the distance matrix */
double getDist( size_t i, size_t j, double * D, size_t d, size_t N) {
 
  size_t k;

  /* note that D is reflexive d(i,i) = 0 */ 
  if ( i == j) return (0.0);

  /* D is symmetric d(i,j) = d(j,i) */ 
  if ( i > j) {
    k = i;
    i = j;
    j = k;
  } 

  /* assume i < j */
  /*
  assert( i < N - 1 );
  assert( j < N     ); 
  */

  return(D[(d+1) * N*(N-1)/2 - (N-i)*(N-i-1)/2 + (j-i-1)]); 
}

/* get a distance without using the design matrix */
double getDistX( 
    size_t i, 
    size_t j, 
    double * x, 
    size_t k, 
    size_t d, 
    size_t N, 
    double (*dist)(double *, double *, size_t)  
  ) {

  /*
  printf(" N*k*d = %d, i*k = %d,  x[N*k*d + i*k] = %f\n",
      (int) N*k*d, (int) i*k, x[N*k*d + i*k] );
  printf(" N*k*d = %d, j*k = %d,  x[N*k*d + j*k] = %f\n",
      (int) N*k*d, (int) j*k, x[N*k*d + i*k] );
      */
      
 
  return( 
    dist( 
      &x[N*k*d + i*k],  
      &x[N*k*d + j*k],  
      k
    )

  );

}


/* this is our standare Euclidian distance function squared */
double squaredEuclidianDist( double * d1, double * d2, size_t k ) {

  size_t i;
  double sum = 0;

  for( i = 0; i < k; i ++)
	  sum += (d1[i] - d2[i])*(d1[i] - d2[i]);

  return( sum );
}

/* this is our standare Euclidian distance function squared */
double squaredEuclidianMeanDist( double * d1, double * d2, size_t k ) {

  size_t i;
  double sum = 0;

  for( i = 0; i < k; i ++)
	  sum += (d1[i] - d2[i])*(d1[i] - d2[i]);

  return( sum/k );
}

/* this is a dissimilarity measure based on equality */
double xorDist( double * d1, double * d2, size_t k ) {

  size_t i;
  double sum = 0;

  for( i = 0; i < k; i ++)
	  sum += (d1[i] != d2[i]);

  return( 2 * sum );
}


/* update the sum of squares index */
void EUpdate( size_t h, size_t i, size_t j, double * E, size_t d, size_t N, size_t * size) {

 double y; 

  y = getDist(h,i,E,d,N) + getDist(h,j,E,d,N);
  E[getIndex(j,h,d,N)]= y; 

}


/* This is an efficient update method for sum of squared distances
 * h is the new index
 * i,j have been merged but not updated
 * D is the prior distance matrix
 */
void JaeDUpdate( size_t h, size_t i, size_t j, double * D, double * E, double * SS, size_t d, size_t N, size_t * size) {

  size_t k;
  double y; 


  if ( i > j) {
    k = i;
    i = j;
    j = k;
  } 
  /* so now i < j */

  /* get the last set of sum of squares */

 y = 
   1/sqrt(size[h] + size[j]) * (getDist(j,h,E,d,N) + SS[j] + SS[h])  
   - 1/sqrt(size[j]) * SS[j] 
   - 1/sqrt(size[h]) *  SS[h]; 
/* This is the sqrt based method from earlier */
/*
 y = 
   (getDist(j,h,E,N) + SS[j] + SS[h])  
   - SS[j] 
   - SS[h]; 

 printf("h: %d\t i: %d\t j: %d\t nij: %d\t nh: %d\t, :",h,i,j, size[j], size[h]); 
 printf("%f, %f, %f, (%f)\n", 
     1/sqrt(size[j] + size[h])*( getDist(j,h,E,N) + SS[j] + SS[h] ) ,
     1/sqrt(size[j])* SS[j]  ,
     1/sqrt(size[h])* SS[h]  ,
     y);  
*/

 D[getIndex(j,h,d,N)] = y; 

}


/* update a single linkage method */
void SingleDUpdate( size_t h, size_t i, size_t j, double * D, double * E, double * SS, size_t d,size_t N, size_t * size ) {

  size_t k;
  double y1, y2; 

  if ( i > j) {
    k = i;
    i = j;
    j = k;
  } 
  /* so now i < j */

  /* get the last set of sum of squares */

  y1 = getDist(j,h,E,d,N);   
  y2 = getDist(i,h,E,d,N);   

  if (y1 > y2) 
    D[getIndex(j,h,d,N)] = y2; 
  else 
    D[getIndex(j,h,d,N)] = y1; 

}



