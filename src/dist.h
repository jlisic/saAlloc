/* header for defined distances */
#ifndef HEADER_DIST

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define HEADER_DIST

#ifdef CLI
#define MATHLIB_STANDALONE
#include "Rmath.h"
#else
#include "R.h"
#endif

/* this needs to be updated for other platforms if not using gcc */
/* and also for running on 32 bit archs */
#ifndef F_SIZET
#define F_SIZET "%zu"
#endif 

/*** Distance matrix functions ***/

/* function that creates the initial distance matrix */
double * createDistMatrix(double * x, size_t k, size_t d, size_t N, double (*dist)( double *, double*, size_t), size_t * ); 

/* updates the distance matrix for duplicate work */
void initUpdateDistMatrix(
  double * ,
  double * ,
  double * ,
  size_t ,
  size_t ,
  size_t * ,
  void (*objFunc) ( size_t , size_t , size_t , double * , double * , double * , size_t, size_t , size_t * )
  ); 

/* function that gets the Matrix (array) index of a specified index in the Matrix */
size_t getIndex( size_t i, size_t j, size_t d, size_t N);

/* function that prints out the dist matrix for debugging */
void printDistMatrix( double * D, size_t d, size_t N); 

/* function that gets a particular value from the distance matrix */
double getDist( size_t i, size_t j, double * D, size_t d, size_t N); 

/* function that gets a particular value without the distance matrix */
double getDistX( size_t i, size_t j, double * x, size_t k, size_t d, size_t N, double (*dist)(double * , double *, size_t)  ); 


/*** specific vector difference functions ***/

/* calculates the basic euclidian distance between two k-vectors */
double squaredEuclidianDist( double * d1, double * d2, size_t k ); 

/* calculates the basic euclidian distance between two k-vectors, and divides by k */
double squaredEuclidianMeanDist( double * d1, double * d2, size_t k ); 

/* calculates the xor between two k-vectors */
double xorDist( double * d1, double * d2, size_t k ); 

/*** group update functions ***/
void EUpdate( size_t h, size_t i, size_t j, double * E, size_t d, size_t N , size_t * size); 
void JaeDUpdate( size_t h, size_t i, size_t j, double * D, double * E, double * SS, size_t d, size_t N, size_t * size ); 
void SingleDUpdate( size_t h, size_t i, size_t j, double * D, double * E, double * SS, size_t d, size_t N, size_t * size ); 


#endif
