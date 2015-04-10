#include "mda.h"


/**************** Functions for 2D Matricies **************************/
/* m - first index
 *  * n - second index
 *   */

/* function that deletes a matrix of type size_t */
size_t ** cv_size_tMatrixCreate( size_t m, size_t n) {
  size_t i;

  size_t ** x = malloc(sizeof(size_t *) * m);

  for( i = 0; i < m; i++) x[i] = calloc( n, sizeof(size_t) );   

  return(x);
}


/* function that deletes a matrix of type size_t */
void cv_MatrixDelete( void ** x, size_t m) {
  size_t i;

  for( i = 0; i < m; i++) free(x[i]);

  free(x);
  x=NULL;

  return;
}


/* function that prints out a matrix of type size_t */
void cv_size_tMatrixPrint(size_t ** X , size_t m, size_t n ) {
  size_t i,j;

  for(i = 0; i < m; i++) {
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < n; j++) {
      Rprintf("%d\t",(int) X[i][j]);
    }
    Rprintf("\n");
  }

  return;
}

/* function that deletes a matrix of type int */
int ** cv_intMatrixCreate( size_t m, size_t n) {
  size_t i,j;

  int ** x = malloc(sizeof(int *) * m);

  for( i = 0; i < m; i++) x[i] = calloc( n, sizeof(int) );   

  return(x);
}




/* function that prints out a matrix of type int */
void cv_intMatrixPrint(int ** X , size_t m, size_t n ) {
  size_t i,j;

  for(i = 0; i < m; i++) {
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < n; j++) {
      Rprintf("%d\t",X[i][j]);
    }
    Rprintf("\n");
  }

  return;
}

/* function that deletes a matrix of type double */
double ** cv_doubleMatrixCreate( size_t m, size_t n) {
  size_t i;

  double ** x = malloc(sizeof(double *) * m);

  for( i = 0; i < m; i++) x[i] = calloc( n, sizeof(double) );   

  return(x);
}




/* function that prints out a matrix of type double */
void cv_doubleMatrixPrint(double ** X , size_t m, size_t n ) {
  size_t i,j;

  for(i = 0; i < m; i++) {
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < n; j++) {
      Rprintf("%4.8f\t",X[i][j]);
    }
    Rprintf("\n");
  }

  return;
}


/**************** Functions for 3D Matricies **************************/
/* m - first index
 * n - second index
 * p - third index
 */


/* function to create size_t 3d matrix */
size_t *** cv_size_t3DMatrixCreate( size_t m, size_t n, size_t p) {

  size_t i,j;

  size_t *** x = malloc(sizeof(size_t **) * m);

  for( i = 0; i < m; i++) { 
    x[i] = malloc(sizeof(size_t *) * n );
    for( j = 0; j < n; j++) x[i][j] = calloc( p, sizeof(size_t) );   
  }

  return(x);
}


/* function to create size_t 3d matrix */
void cv_3DMatrixDelete( void *** x, size_t m, size_t n) {

  size_t i,j;

  for( i = 0; i < m; i++) { 
    for( j = 0; j < n; j++) free(x[i][j]);
    free(x[i]);
  }

  free(x);
  x = NULL;
  return;
}


/* function that prints out the 3dMatrix */
void cv_size_t3DMatrixPrint( size_t ***x, size_t m, size_t n, size_t p) {
  size_t i;

  for( i = 0; i < m; i++)  cv_size_tMatrixPrint(x[i], n, p ); 

  return;
}


/* function to create int 3d matrix */
int *** cv_int3DMatrixCreate( size_t m, size_t n, size_t p) {

  size_t i,j;

  int *** x = calloc(m, sizeof(int **));

  for( i = 0; i < m; i++) { 
    x[i] = calloc(n , sizeof(int *) );
    for( j = 0; j < n; j++) x[i][j] = calloc( p, sizeof(int) );   
  }

  return(x);
}




/* function that prints out the 3dMatrix */
void cv_int3DMatrixPrint( int ***x, size_t m, size_t n, size_t p) {
  size_t i;

  for( i = 0; i < m; i++)  cv_intMatrixPrint(x[i], n, p ); 

  return;
}


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




/* function that prints out the 3dMatrix */
void cv_double3DMatrixPrint( double ***x, size_t m, size_t n, size_t p) {
  size_t i;

  for( i = 0; i < m; i++)  cv_doubleMatrixPrint(x[i], n, p ); 

  return;
}
        


/**************** Functions for 4D Matricies **************************/
/* m - first index
 * n - second index
 * p - third index
 * q - fourth index
 */


/* function to create size_t 4d matrix */
size_t **** cv_size_t4DMatrixCreate( size_t m, size_t n, size_t p, size_t q) {

  size_t i;

  size_t **** x = malloc(sizeof(size_t ***) * m);

  for( i = 0; i < m; i++) x[i] = cv_size_t3DMatrixCreate( n, p, q);

  return(x); 
}


/* function to create size_t 4d matrix */
void cv_4DMatrixDelete( void **** x, size_t m, size_t n, size_t p) {

  size_t i;

  for( i = 0; i < m; i++) cv_3DMatrixDelete( x[i],  n, p);

  free(x);
  x = NULL;
  return;
}


/* function that prints out the 4dMatrix */
void cv_size_t4DMatrixPrint( size_t ****x, size_t m, size_t n, size_t p, size_t q) {

  size_t i;

  for( i = 0; i < m; i++)  cv_size_t3DMatrixPrint(x[i], n, p, q); 

  return;
}



/* function to create int 4d matrix */
int **** cv_int4DMatrixCreate( size_t m, size_t n, size_t p, size_t q) {

  size_t i;

  int **** x = malloc(sizeof(int ***) * m);

  for( i = 0; i < m; i++) x[i] = cv_int3DMatrixCreate( n, p, q);

  return(x); 
}




/* function that prints out the 4dMatrix */
void cv_int4DMatrixPrint( int ****x, size_t m, size_t n, size_t p, size_t q) {

  size_t i;

  for( i = 0; i < m; i++)  cv_int3DMatrixPrint(x[i], n, p, q); 

  return;
}



/* function to create double 4d matrix */
double **** cv_double4DMatrixCreate( size_t m, size_t n, size_t p, size_t q) {

  size_t i;

  double **** x = malloc(sizeof(double ***) * m);

  for( i = 0; i < m; i++) x[i] = cv_double3DMatrixCreate( n, p, q);

  return(x); 
}




/* function that prints out the 4dMatrix */
void cv_double4DMatrixPrint( double ****x, size_t m, size_t n, size_t p, size_t q) {

  size_t i;

  for( i = 0; i < m; i++)  cv_double3DMatrixPrint(x[i], n, p, q); 

  return;
}




/************************************ GENERIC FUNCTION ************************/
/* it is assumed that all pointers are of size size_t */
/* input */
/* dim - number of dimensions, max 4 
 * type - type to create
 *
 */

/* CREATE ***************************************************************************************/
void * createMDA( size_t type, size_t dim, ...) {

  va_list ap;
  size_t j;

  size_t dims[4];

  va_start(ap, dim);


  if( dim > 4) {
    Rprintf("Error: dim exceeds bounds\n");
    return(NULL);
  }
  if( dim < 1) {
    Rprintf("Error: dim exceeds bounds\n");
    return(NULL);
  }


  for( j = 0; j < dim; j++) {
    dims[j] = va_arg( ap, size_t) ;
    //printf("check: %d\n", (int) dims[j] );
  } 


  switch( type ) {

    /* size_t */
    case 0:
      switch( dim ) {

        case 1:
          return( (void *)  calloc( dims[0], sizeof(size_t) ) ); 
        
        case 2:
          return( (void *)  cv_size_tMatrixCreate( dims[0], dims[1]) ); 

        case 3:
          return( (void *)  cv_size_t3DMatrixCreate( dims[0], dims[1], dims[2]) ); 
        
        case 4:
          return( (void *)  cv_size_t4DMatrixCreate( dims[0], dims[1], dims[2], dims[3]) ); 

        default: return(NULL);

      } 
      /* integer */
    case 1:
      switch( dim ) {

        case 1:
          return( (void *)  calloc( dims[0], sizeof(int) ) ); 
        
        case 2:
          return( (void *)  cv_intMatrixCreate( dims[0], dims[1]) ); 

        case 3:
          return( (void *)  cv_int3DMatrixCreate( dims[0], dims[1], dims[2]) ); 

        case 4:
          return( (void *)  cv_int4DMatrixCreate( dims[0], dims[1], dims[2], dims[3]) ); 

        default: return(NULL);

      } 
      /* double */
    case 2:
      switch( dim ) {

        case 1:
          return( (void *)  calloc( dims[0], sizeof(double) ) ); 
        
        case 2:
          return( (void *)  cv_doubleMatrixCreate( dims[0], dims[1]) ); 

        case 3:
          return( (void *)  cv_double3DMatrixCreate( dims[0], dims[1], dims[2]) ); 
        
        case 4:
          return( (void *)  cv_double4DMatrixCreate( dims[0], dims[1], dims[2], dims[3] ) ); 

        default: return(NULL);

      } 

   
      /* unspecified */     
    default: return(NULL);

  } 


  /* on error condition */
  return NULL;
}




/* DELETE ***************************************************************************************/
void deleteMDA( void * x, size_t dim, ...) {

  va_list ap;
  size_t j;

  size_t dims[4];

  va_start(ap, dim);


  if( dim > 3) {
    Rprintf("Error: dim exceeds bounds\n");
    return;
  }


  for( j = 0; j < dim; j++) {
    dims[j] = va_arg( ap, size_t) ;
    //printf("check: %d\n", (int) dims[j] );
  } 


  switch( dim ) {

    case 1:
      free(x);
      break; 
     
    case 2:
      cv_MatrixDelete( (void **) x, dims[0]); 
      break;

    case 3:
      cv_3DMatrixDelete( (void ***) x, dims[0], dims[1]); 
      break;
    
    case 4:
      cv_4DMatrixDelete( (void ****) x, dims[0], dims[1], dims[2]); 
      break;

    default: return;

  } 

  return;
}







/* PRINT ***************************************************************************************/
void printMDA( void * x, size_t type, size_t dim, ...) {

  va_list ap;
  size_t i,j;

  size_t dims[4];

  va_start(ap, dim);


  if( dim > 4) {
    Rprintf("Error: dim exceeds bounds");
    return;
  }
  if( dim < 1) {
    Rprintf("Error: dim exceeds bounds");
    return;
  }

  for( j = 0; j < dim; j++) {
    dims[j] = va_arg( ap, size_t) ;
    //printf("check: %d\n", (int) dims[j] );
  } 


  switch( type ) {

    /* size_t */
    case 0:
      switch( dim ) {

        case 1:
          for( i = 0; i < dims[0]; i++) Rprintf("%d ", (int) ((size_t *) x)[i]  ); 
          Rprintf("\n");
          Rprintf("\n");
          break;
        
        case 2:
          cv_size_tMatrixPrint( (size_t **) x, dims[0], dims[1]); 
          break;

        case 3:
          cv_size_t3DMatrixPrint( (size_t ***) x, dims[0], dims[1], dims[2]); 
          break;
        
        case 4:
          cv_size_t4DMatrixPrint( (size_t ****) x, dims[0], dims[1], dims[2], dims[3] ); 
          break;

        default: return;

      } 
      break;
      /* integer */
    case 1:
      switch( dim ) {

        case 1:
          for( i = 0; i < dims[0]; i++) Rprintf("%d ", ((int *) x)[i]  ); 
          Rprintf("\n");
          Rprintf("\n");
          break;
        
        case 2:
          cv_intMatrixPrint( (int **) x, dims[0], dims[1]); 
          break;

        case 3:
          cv_int3DMatrixPrint( (int ***) x, dims[0], dims[1], dims[2]); 
          break;
        
        case 4:
          cv_int4DMatrixPrint( (int ****) x, dims[0], dims[1], dims[2], dims[3] ); 
          break;

        default: return;
      } 
      break;
      /* double */
    case 2:
      switch( dim ) {

        case 1:
          for( i = 0; i < dims[0]; i++) Rprintf("%4.8f ", (float)  ((double *) x)[i]  ); 
          Rprintf("\n");
          Rprintf("\n");
          break;
        
        case 2:
          cv_doubleMatrixPrint( (double **) x, dims[0], dims[1]); 
          break;

        case 3:
          cv_double3DMatrixPrint( (double ***) x, dims[0], dims[1], dims[2]); 
          break;
        
        case 4:
          cv_double4DMatrixPrint( (double ****) x, dims[0], dims[1], dims[2], dims[3] ); 
          break;

        default: return;

      } 
      break;

   
      /* unspecified */     
    default: return;

  } 


  /* on error condition */
  return;
}














/* test function */
/*
int main( ) {

 int *** x;

 x = (int ***) createMDA( MDA_INT , 3, 2, 3, 4);

 printMDA( (void *) x , MDA_INT,  3, 2, 3, 4);
 
 deleteMDA( (void *) x , 3, 2, 3);


 return 0;
}
*/

