/* Copyright (c) 2015-2016  Jonathan Lisic 
 * Last edit: 16/08/24 - 09:20:44
 * License: GPL (>=2) 
 */  

#include <stdio.h> 
#include <stdlib.h>

#include "R.h"
#include "Rinternals.h"
#include "Rmath.h"
#include <R_ext/Rdynload.h>

#include "mainSubstrata2.h"
#include "mainMinCV.h"
#include "alloc.h"

/***********************************/
/* Register SO's                   */



static R_NativePrimitiveArgType R_minCV_t[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, 
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, 
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_substrata2_t[] = {
      REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, 
      INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_sampleAlloc_t[] = {
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, 
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, 
  INTSXP 
};

static const R_CMethodDef cMethods[] = {
     {"R_minCV", (DL_FUNC) &R_minCV, 24, R_minCV_t},
     {"R_substrata2", (DL_FUNC) &R_substrata2, 15, R_substrata2_t},
     {"R_sampleAlloc", (DL_FUNC) &R_sampleAlloc, 17, R_sampleAlloc_t},
        {NULL, NULL, 0, NULL}
};

void R_init_myLib(DllInfo *info)
{
     R_registerRoutines(info, cMethods, NULL, NULL, NULL);
     R_useDynamicSymbols(info, TRUE); 
}


