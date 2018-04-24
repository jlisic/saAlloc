#ifndef HEADER_ALLOC
#define HEADER_ALLOC

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mda.h"
#include "cv.h"
#include "R.h"
#include "Rmath.h"


#ifndef SA_GETINDEX
  #define SA_GETINDEX(x) (size_t) ( runif(0.0,1.0) * (double) (x) ) 
#endif

/*
double fake_runif( double a, double b ); 
int fake_rsample( int b ); 
*/

/* get a strata for i to move to */
/* the selected strata can be equal to the current strata */
/* the function either returns H (no move possible */
/* or the selected stratum */
size_t alloc_getMoveStrata( 
    size_t i, 
    size_t * I, 
    double * probMatrix,
    size_t N,
    size_t H
    ); 


void alloc_sampleSizeChange (
    double * cvInit,         
    size_t N, 
    size_t K, 
    size_t H, 
    size_t R, 
    size_t J, 
    size_t *** Domain,
    double *** var, 
    size_t * Nh,
    double * nh,
    double ** Total,
    double *** locationAdj,
    double *** scaleAdj,
    double * Target,
    double * penalty,
    double p,
    size_t preserveSatisfied,
    size_t iter,
    double * a,
    size_t fpc
    ); 


void R_sampleAlloc (
  double * totalDouble,    /* the column major listing of means      J*R  */
  double * varDouble,   /* the column major listing of variances  J*H*R  */
  int * JInt,      /* number of variables */
  int * HInt,      /* number of strata */
  int * RInt,      /* number of observations of each PSU (e.g. years) */
  int * iterInt,   /* number of interations                  1      */
  int * domainInt, /*                                        H*J  */
  double * target,   /* target CV*/
  double * locationAdjDouble,  
  double * scaleAdjDouble,
  double * pDouble,
  double * penalty,  /* length J */
  double * nh,      /* sample size by stratum */
  int * NhInt,      /* pop size by stratum */
  double * a,  // what we return
  double * cooling,
  int * fpc
);               


#endif
