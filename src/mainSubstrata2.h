/* Jonathan Lisic */
/* standard libraries */

#ifndef MAINSUBSTRATA2
#define MAINSUBSTRATA2

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "R.h"

/* program headers */
#include "dist.h"
//#include "file.h"
#include "sa.h"

#include "substrata2.h"

/* this needs to be updated for other platforms if not using gcc */
/* and also for running on 32 bit archs */
#define F_SIZET "%zu"

/**************************** MAIN R ******************************************/
/* read in input and output from R                                            */
/**************************** MAIN R ******************************************/
void R_substrata2 (
      double * x,      /* the column major listing of items      k*n*d  */
      int * kInt,      /* the length of an individual item       1      */
      int * dInt,      /* number of commodities                  1      */
      int * NInt,      /* number of individual items             1      */
      int * iterInt,   /* number of interations                  1      */
      int * IInt,      /* destructive: init/final stratification n      */
      double * Q,      /* cost value of each commodity           d      */
      double * adminDataDbl,   /* Administrative data double             ANDbl  */
      int * adminDataInt,      /* Administrative data integer            ANInt  */
      int * ANIntInt,  /* number of double admin recordds        1      */
      int * ANDblInt,  /* number of integer admin recors         1      */
      int * dup,       /* UNTESTED/UNUSED dup vector             n      */
      double * costChangeDbl, /* (costChange,T) per item          iter * 2 */
      int * auxFunctionIterInt, /* how often to run the auxiliary function 1 */
      int * costChangeSizeInt      /* size of a cost change element           1 */
    );


#endif
