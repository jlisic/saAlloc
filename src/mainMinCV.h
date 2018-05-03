/* Jonathan Lisic */
/* standard libraries */

#ifndef MAINMINCV
#define MAINMINCV 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "R.h"

/* program headers */
#include "dist.h"
//#include "file.h"
#include "sa.h"

#include "minCV.h"

/* this needs to be updated for other platforms if not using gcc */
/* and also for running on 32 bit archs */
#define F_SIZET "%zu"




/**************************** MAIN R ******************************************/
/* read in input and output from R                                            */
/*  -- old --
1      double * x,       the column major listing of items      k*n*d  
2      int * kInt,       the length of an individual item       1     
3      int * dInt,       number of commodities                  1      
4      int * NInt,       number of individual items             1      
5      int * iterInt,    number of interations                  1      
6      int * IInt,       destructive: init/final allocation     n      
7      double * Q,       cost value of each commodity           d      
8      double * adminDataDbl,   / Administrative data double             ANDbl  
9      int * adminDataInt,      / Administrative data integer            ANInt  
10      int * ANIntInt,  / number of double admin recordds        1      
11      int * ANDblInt,  / number of integer admin recors         1      
12      int * dup,       / UNTESTED/UNUSED dup vector             n      
13      double * costChangeDbl,      / (costChange,T) per item          iter * 3 
14      double * doubleSampleSize,   / output sample size                        
15      int * auxFunctionIterInt,    / how often to run the auxiliary function 1 
16      int * costChangeSizeInt      / size of a cost change element           1 
*/
/**************************** MAIN R ******************************************/

void R_minCV (
    double * x,      /* the column major listing of items              K*N*R          1*/
    int * NInt,      /* number of PSUs                                                2*/
    int * JInt,      /* number of domains                                             3*/
    int * KInt,      /* number of variables                                           4*/
    int * HInt,      /* number of strata                                              5*/
    int * RInt,      /* number of observations of each PSU                            6*/
    int * iterSampleSizeInt,   /* number of interations of the sampling algorithm  1  7*/
    int * IInt,      /* destructive: init/final allocation               n            8*/
    int * domainInt, /*                                              H*K*J            9*/
    double * prob,       /*                                            N             10*/ 
    double * probMatrix, /*                                            N*H           11*/ 
    double * target,           /*                                                    12*/
    double * locationAdjDouble, /*                                                   13*/
    double * scaleAdjDouble,    /*                                                   14*/
    double * pDouble,           /*                                                   15*/
    double * penalty,           /*                                                   16*/
    double * nh,      /* sample size by stratum                             H        17*/
    double * costChangeDbl,      /* data saved at each iteration    5 * H + K        18*/
    double * temp,               /* temperature for SA                               19*/
    int * auxFunctionIterInt,    /* how often to run the auxiliary function 1        20*/
    int * costChangeSizeInt,     /* size of a cost change element           1        21*/
    int * iterInt,               /*                                                  22*/
    int * preservedSatisfiedInt,
    int * fpc                    /* finite population correction factor              23*/
  );

#endif
