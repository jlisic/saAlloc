/* This is a simple Multi-Dimensional Array Library written to support saAlloc */

/* In the absence of C generic functions, there are three types that can be
 * specified by the MDA library:
 * 
 * index  ENUM       class 
 * 0      MDA_SIZE_T size_t 
 * 1      MDA_INT    int 
 * 2      MDA_DOUBLE double 
 *
 * All of these MDA's are constructed by sequential calloc's for each level of 
 * the pointers.  Therefore a three dimensional array can be accessed through
 * the standard A[a][b][c].
 *
 * There are three operations supported in this library
 * - construction
 * - deconstruction
 * - printing
 *
 * Varargs are used to allow for a single command to work for all types.
 *
 *
 * void * createMDA( size_t type, dimensions, ...)
 *
 * where type is a type such as double listed above (0,1,2)
 * dimensions is the number of dimensions. 
 * ... is one less than the max index, of size at least 1. 
 *
 * Note that this function returns void * so a cast is required to get the desired
 * type.
 *
 *
 * void deleteMDA( void * x, dimensions, ...)
 *
 * x is the object being deleted (must be cast to void *)
 * dimensions is the number of dimensions. 
 * ... is one less than the max index, of size at least 1, do not include the last 
 * max index. 
 *
 * 
 * void printMDA( void * x, dimensions, ...)
 *
 * x is the object being printed (must be cast to void *)
 * dimensions is the number of dimensions. 
 * ... is one less than the max index, of size at least 1.
 *
 *
 * An example is provided below 


 // create an MDA in 3 dimensions 2x3x4
 int main( ) {

 int *** x;

 // create MDA 
 x = (int ***) createMDA( MDA_INT , 3, 2, 3, 4);

 // print MDA 
 printMDA( (void *) x , MDA_INT, 3, 2, 3, 4);
 
 // delete MDA  
 deleteMDA( (void *) x , 3, 2, 3);


 return 0;
 }


 */ 

#ifndef HEADER_MDA
#define HEADER_MDA

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


#define MDA_SIZE_T 0 
#define MDA_INT 1 
#define MDA_DOUBLE 2



/* Matrix helper functions */
size_t ** cv_size_tMatrixCreate( size_t m, size_t n); 
void cv_size_tMatrixPrint(size_t ** X , size_t m, size_t n ); 

int ** cv_intMatrixCreate( size_t m, size_t n); 
void cv_intMatrixPrint(int ** X , size_t m, size_t n ); 

double ** cv_doubleMatrixCreate( size_t m, size_t n); 
void cv_doubleMatrixPrint(double ** X , size_t m, size_t n ); 

void cv_MatrixDelete( void ** x, size_t m); 


/* 3D Matrix helper functions */
size_t *** cv_size_t3DMatrixCreate( size_t m, size_t n, size_t p); 
void cv_size_t3DMatrixPrint( size_t ***x, size_t m, size_t n, size_t p); 

int *** cv_int3DMatrixCreate( size_t m, size_t n, size_t p); 
void cv_int3DMatrixPrint( int ***x, size_t m, size_t n, size_t p);

double *** cv_double3DMatrixCreate( size_t m, size_t n, size_t p); 
void cv_double3DMatrixPrint( double ***x, size_t m, size_t n, size_t p); 

void cv_3DMatrixDelete( void *** x, size_t m, size_t n); 


/* 4D Matrix helper functions */
size_t **** cv_size_t4DMatrixCreate( size_t m, size_t n, size_t p, size_t q); 
void cv_size_t4DMatrixPrint( size_t ****x, size_t m, size_t n, size_t p, size_t q); 

int **** cv_int4DMatrixCreate( size_t m, size_t n, size_t p, size_t q); 
void cv_int4DMatrixPrint( int ****x, size_t m, size_t n, size_t p, size_t q); 

double **** cv_double4DMatrixCreate( size_t m, size_t n, size_t p, size_t q); 
void cv_double4DMatrixPrint( double ****x, size_t m, size_t n, size_t p, size_t q); 

void cv_4DMatrixDelete( void **** x, size_t m, size_t n, size_t p); 



/* create MDA function */
void * createMDA( size_t type, size_t dim, ...); 

/* delete MDA function */
void deleteMDA( void * x, size_t dim, ...); 

/* print MDA function */
void printMDA( void * x, size_t type, size_t dim, ...); 

#endif
