#ifndef HEADER_CV
#define HEADER_CV

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"


/******* FUNCTIONS FOR 2D ARRAYS *********/

/* function to create matrix (2d array) */
/* for x bounded by x[m][n] */
double ** cv_doubleMatrixCreate(size_t m, size_t n );

/* function to delete matrix (2d array) */
/* for x bounded by x[m][n] */
void cv_doubleMatrixDelete(double ** x, size_t m);

/* function to print matrix (2d array) */
/* for x bounded by x[m][n] */
void cv_doubleMatrixPrint(double ** x , size_t m, size_t n );


/******* FUNCTIONS FOR 3D ARRAYS *********/

/* function to create MDA */
/* for x bounded by x[m][n][p] */
double *** cv_double3DMatrixCreate( size_t m, size_t n, size_t p); 

/* function to delete MDA */
/* for x bounded by x[m][n][p] */
void cv_double3DMatrixDelete( double ***, size_t m, size_t n); 

/* function to print MDA */
/* for x bounded by x[m][n][p] */
void cv_double3DMatrixPrint( double ***, size_t m, size_t n, size_t p); 


/******* FUNCTIONS FOR MEAN AND VARIANCE *********/

/* function that creates a mean MDA (matrix) [commodity][strata][obs] */
double *** cv_createMeanMatrix( size_t * I, size_t dN, size_t N, size_t kN, size_t H, double * x, size_t * Nh); 

/* function that creates a variance MDA (matrix) [commodity][strata][obs] */
double *** cv_createVarMatrix( size_t * I, size_t dN, size_t N, size_t kN, size_t H, double * x, double *** mu, size_t * Nh) ; 


/******* FUNCTIONS FOR UPDATING MEAN AND VARIANCE *********/

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
  ); 


/******* FUNCTION FOR CALCULATING TOTAL *********/
double ** cv_createTotalMatrix( 
    size_t D, 
    size_t N, 
    size_t K, 
    size_t H,
    double ***mu, 
    size_t * Nh
    ); 


/******* FUNCTION FOR CALCULATING OBJECTIVE FUNCTION *********/
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
  ); 


#endif
