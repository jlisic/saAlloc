#ifndef HEADER_CV
#define HEADER_CV

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
 #include "omp.h"
#endif

#include "mda.h"



/******* FUNCTIONS FOR MEAN AND VARIANCE *********/

/* function that creates a mean MDA [coi][strata][obs] */
double *** cv_createMeanMatrix( 
    size_t * I,  /* index */
    size_t N,    /* PSUs  */
    size_t K,    /* CoI   */
    size_t H,    /* Strata */
    size_t R,    /* obs */
    double * x,  /* data CoI x I x obs */ 
    size_t * Nh  /* strata size */
    ); 

/* function that creates a var MDA [coi][strata][obs] */
double *** cv_createVarMatrix( 
    size_t * I, 
    size_t N, 
    size_t K, 
    size_t H, 
    size_t R, 
    double * x, 
    double *** mu, 
    size_t * Nh 
    ); 


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
    size_t J,
    size_t ***Domain,
    double ***mu, 
    size_t * Nh
    ); 

/******* FUNCTION FOR APPLYING LOCATION SCALE ADJUSTMENTS  *********/
double *** cv_createDomainMDA( 
    size_t D, 
    size_t N, 
    size_t K, 
    size_t H,
    size_t J,
    size_t ***Domain,
    double ***x, 
    double *** locationAdj, /* if NULL then not used */
    double *** scaleAdj     /* if NULL then not used */
    ); 


/******* FUNCTION FOR CALCULATING THE CV *********/
double * cv_calcCV(
    double * cv,   /* if cv != NULL create it, otherwise handle in place */
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
    double *** locationAdj, /* IF NULL then not used */
    double *** scaleAdj,    /* IF NULL then not used */
    size_t fpc
  ); 

/******* FUNCTION FOR CALCULATING OBJECTIVE FUNCTION *********/
double cv_objectiveFunction( 
    double * cv,   // if cv != NULL create it, otherwise handle in place 
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
    double *** locationAdj, /* IF NULL then not used */
    double *** scaleAdj,    /* IF NULL then not used */
    double * Target,       /* IF NULL then not used */
    double * penalty,      /* IF NULL then not used */
    double p,
    size_t evaluateOnly, // option to not construct CV, under this condition CV cannot be null 
    size_t fpc
  ); 


/* handle the objective function with a constraint on preserving prior met constraints */
double cv_objectiveFunctionCompare( 
    double * cv,         //if cv != NULL create it, otherwise handle in place 
    double * cvPrior,   
    double * obj,
    double * objPrior,
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
    size_t evaluateOnly, // option to not construct CV, under this condition CV cannot be null 
    size_t preserveSatisfied,
    size_t fpc
  ); 


/* update matrix with the addition of handling 2 additional means, this is needed by min CV */
void cv_updateMatrix2( 
    size_t * I,     /* this is the current stratification, not the modified one */ 
    size_t N, 
    size_t K, 
    size_t H, 
    size_t R, 
    double * x, 
    double * x1,
    double * x2,
    double *** mu, 
    double *** mu1, 
    double *** mu2, 
    double *** var, 
    size_t * Nh,
    size_t moveObs,       /* this is the observation to move */
    size_t moveObsStratum /* this is the stratum to move the observation to */
  ); 


#endif
