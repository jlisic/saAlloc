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
    double * a
    ); 



#endif
