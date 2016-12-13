#include "alloc.h"


/*
//place holder
double fake_runif( double a, double b ) {
  return rand() / (RAND_MAX + 1.0) * (b - a) + a;
}

//place holder
int fake_rsample( int b ) {
  return ( (int) rand() ) % b ;
}
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
    ) {

  double totalProb = 0;
  double total;
  size_t Hj;

  /* get the total weight of all the strata */
  for( Hj=0; Hj < H; Hj++) { 
    totalProb += probMatrix[i*H + Hj];
  }
  /* if the total probability remaininig is 0 quit */
  if( totalProb == 0 ) return( H );   

  /* randomly select the strata proportional to its weight */
  //totalProb = runif(0,totalProb); 
  totalProb = runif(0.0,totalProb); 
  total = 0;

  /* find the strata that corresponds to the selected sampling weight */ 
  for( Hj=0; Hj < H; Hj++) {
    
    total += probMatrix[i*H + Hj];

    if( total > totalProb) break; /* unless all the probability is 0 then units are sampled on an 'open' interval */ 
  }
  
  return( Hj );
}



/* function to select an index */
/*
size_t minCV_getIndex( double * prob, double totalProbability ) {

  size_t index = 0;
  double total = 0;
  double searchProbability = runif(0,totalProbability); // select a value uniform across the sum of probabilities 
      
  total = prob[0]; 
  while( total < searchProbability ){ 
      index++;
      total += prob[index];
  }

  return( index );
}
*/


/* This is a function to update the sample size                              */
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
    ) { 


  double * cv = calloc( J, sizeof(double) );  // a place to store changes
  double * test_nh = calloc( H, sizeof(double) );  // a place to store changes
  size_t i,h,j;         // iteration
  double minSampleSize = 2.0; // minimum sample size
  double delta;
  size_t Hi, Hj;

  /* sanity check */
  double nhSumStart;
  double nhSumStop;

  for(h = 0, nhSumStart = 0; h < H; h++) nhSumStart += nh[h]; 
 

  // if there is nothing to do, do nothing 
  if( iter == 0) return;
 
  // copy values over
  for(h = 0; h < H; h++) test_nh[h] = nh[h]; 
  
  
  // iterate a fixed number of times 
  for( i = 0 ; i < iter; i++ ) {

    // 1.0 randomly select a strata to move from 
    Hi = SA_GETINDEX(H);


    // only proceed if stratum Hi can be made smaller 
    if( nh[Hi] < minSampleSize + 1 ) {
      if(a != NULL) a[i] = nan("");
      continue;
    }
    
      
    // 2.0 calculate objective function change for moving to each strata 
    Hj = SA_GETINDEX(H);
 
    // if the exchange would make the sample size too big, we don't do it 
    if( nh[Hj] + 1 > Nh[Hj] ) {
      if(a != NULL) a[i] = nan("");
      continue; 
    }


    // if they match don't evaluate
    if( Hj != Hi) {
        
      test_nh[Hi] -= 1.0; // decrement stratum Hi 
      test_nh[Hj] += 1.0; // increment stratum Hj
      
      //printf("%d:",(int) i);  
      //for(j =0; j <H; j++) printf("%d,",(int) test_nh[j]);

      delta = cv_objectiveFunctionCompare( 
        cv, 
        cvInit, 
        NULL,
        NULL,
        N, 
        K, 
        H, 
        R, 
        J, 
        Domain, 
        var, 
        Nh, 
        test_nh, 
        Total, 
        locationAdj, 
        scaleAdj, 
        Target, 
        penalty, 
        p, 
        0, 
        preserveSatisfied,
        fpc
      ); 
      //printf("\n");


      if( delta <= 0 ) { 
        // update the cv and strata assignment
        for(j = 0; j < J; j++) cvInit[j] = cv[j]; 
        nh[Hi] = test_nh[Hi];
        nh[Hj] = test_nh[Hj];
      } else {
        // change back
        test_nh[Hi] = nh[Hi];
        test_nh[Hj] = nh[Hj];
      }
      
      // cover the Hj == Hi case
    } else {
      delta = 0;
    }
    

   // sanity check 
    for(h = 0, nhSumStop = 0; h < H; h++) nhSumStop += nh[h]; 
    if( nhSumStop != nhSumStart ) Rprintf("alloc: nhSumStop = %f , nhSumStart = %f\n", nhSumStop, nhSumStart);

    // r
    // return min_delta
    if( a != NULL ) a[i] = delta;
  }
    

  // copy values over
  //for(h = 0; h < H; h++) nh[h] = test_nh[h]; 
  free(test_nh);
  free(cv);

  return;
}









/**************** R FUNC ****************************/


void R_sampleAlloc (
  double * x,      /* the column major listing of items      K*N*R  */
  int * NInt,      /* number of PSUs */
  int * JInt,      /* number of domains */
  int * KInt,      /* number of variables */
  int * HInt,      /* number of strata */
  int * RInt,      /* number of observations of each PSU */
  int * iterInt,   /* number of interations                  1      */
  int * IInt,      /* destructive: init/final allocation     n      */
  int * domainInt, /*                                        H*K*J  */
  double * probMatrix, /*                                    N*H    */ 
  double * target,
  double * locationAdjDouble,
  double * scaleAdjDouble,
  double * pDouble,
  double * penalty, 
  double * nh,      /* sample size by stratum */
  double * a,  // what we return
  double * cooling,
  int * fpc
) {              

  size_t N = *NInt;
  size_t J = *JInt;
  size_t K = *KInt;
  size_t H = *HInt;
  size_t R = *RInt;

  size_t i,h,j,k,r;
  size_t iter = *iterInt;
  double p = * pDouble;
  double ** total;
  double * cvInit;

  size_t * Nh = (size_t *) createMDA( MDA_SIZE_T, 1, H); 
  size_t * I = (size_t *) createMDA( MDA_SIZE_T, 1, N); 
  
  for(i = 0; i < N; i++) {
    I[i] = IInt[i];
    Nh[I[i]]++;
  }

  
  double *** locationAdj;
  double *** scaleAdj; 

  // Mean
  size_t *** domain = (size_t ***) createMDA( MDA_SIZE_T, 3, H, K, J); 
  for( h=0; h < H; h++)
    for( k=0; k < K; k++)
      for( j=0; j < J; j++) domain[h][k][j] = (size_t) domainInt[ h*K*J + k*J + j ];

  
  // location adjustement
  if( locationAdjDouble[0] >= 0 ) { 
    locationAdj = cv_createMeanMatrix(I, N, K, H, R, x, Nh); 
  } else {
    locationAdj = NULL;
  }
  
  // scale adjustmentf
  if( scaleAdjDouble[0] >= 0 ) { 
    scaleAdj = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
    for( k=0; k < K; k++) 
      for( h=0; h < H; h++) 
        for( r=0; r < R; r++) scaleAdj[k][h][r] = scaleAdjDouble[k*H*R + h*R + r];
  } else {
    scaleAdj = NULL;
  }

  // Mean
  double *** mu = cv_createMeanMatrix(I, N, K, H, R, x, Nh); 

  // Variance
  double *** var = cv_createVarMatrix(I, N, K, H, R,  x, mu, Nh); 

  // Total 
  total =  cv_createTotalMatrix( N, K, H, R, J, (size_t ***) domain, mu, Nh);

  
  // CV
  cvInit = cv_calcCV( NULL, N, K, H, R, J, domain, var, Nh, nh, total, locationAdj, scaleAdj, (size_t) *fpc ); 
 
  // get change in allocation 
  alloc_sampleSizeChange (
      cvInit, N, K, H, R, J, domain, var, Nh, nh, total, locationAdj, scaleAdj, target, penalty, p, 0, iter, a, (size_t) * fpc);
  

  
  printMDA( (void *) penalty, MDA_DOUBLE, 1, J); 

  // clean up 
  deleteMDA( (void * ) mu, 3, K, H); 
  deleteMDA( (void *) var, 3, K, H); 
  deleteMDA( (void *) total, 2, J); 
  deleteMDA( (void *) cvInit, 1); 
  deleteMDA( (void *) Nh, 1); 
  deleteMDA( (void *) I, 1); 

  if( locationAdjDouble[0] >= 0 ) {
    deleteMDA( (void * ) locationAdj, 3, K, H); 
  }
  if( scaleAdjDouble[0] >= 0 ) {
    deleteMDA( (void *) scaleAdj, 3, K,H); 
  }

  return;
}






