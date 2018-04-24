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
    size_t * Nh,   // what gets passed in is the candidate
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
    
  free(test_nh);
  free(cv);

  return;
}





/**************** R FUNC ****************************/


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
) {              

  size_t J = *JInt;
  size_t K = *JInt; /* not really used anymore */
  size_t H = *HInt;
  size_t R = *RInt;
  size_t N = 0;

  size_t i,h,j,k,r;
  size_t iter = *iterInt;
  double p = * pDouble;
  double * cvInit = NULL;

  size_t * Nh = (size_t *) createMDA( MDA_SIZE_T, 1, H); 
  for(h=0;h<H;h++) {
    Nh[h] = NhInt[h];
    N += Nh[h];
  }
  
  double *** locationAdj;
  double *** scaleAdj; 

  // Mean
  size_t *** domain = (size_t ***) createMDA( MDA_SIZE_T, 3, H, K, J); 
  for( h=0; h < H; h++)
    for( k=0; k < K; k++)
      for( j=0; j < J; j++) domain[h][k][j] = (size_t) domainInt[ h*K*J + k*J + j ];

  
  // location adjustment
  // k is domain (including commodities) 
  // h is strata
  // r is replicate
  // x_k,h,t
  //  x_0,0,t  x_0,0,t+1  x_0,0,t+2 
  if( locationAdjDouble[0] >= 0 ) { 
    locationAdj = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
    for( k=0; k < K; k++) 
      for( h=0; h < H; h++) 
        for( r=0; r < R; r++) locationAdj[k][h][r] = locationAdjDouble[k*H*R + h*R + r];
  } else {
    locationAdj = NULL;
  }
  
  // scale adjustment
  if( scaleAdjDouble[0] >= 0 ) { 
    scaleAdj = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
    for( k=0; k < K; k++) 
      for( h=0; h < H; h++) 
        for( r=0; r < R; r++) scaleAdj[k][h][r] = scaleAdjDouble[k*H*R + h*R + r];
  } else {
    scaleAdj = NULL;
  }

  // Variance (commodity, strata, rep)
  double ** total = (double **) createMDA( MDA_DOUBLE, 2, K, R); 
  for( k=0; k < K; k++) 
      for( r=0; r < R; r++) total[k][r] = totalDouble[k*R + r];

  // Variance (actually S^2)  (commodity, strata, rep)
  double *** var = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
  for( k=0; k < K; k++) 
    for( h=0; h < H; h++) 
      for( r=0; r < R; r++) var[k][h][r] = varDouble[k*H*R + h*R + r];

  /* for debug */
  /*
  printf("K = %d, J = %d, H= %d, R = %d\n", (int) K, (int) J, (int) H, (int) R);  
  printf("total\n"); 
  printMDA( (void *) total, MDA_DOUBLE, 2, J, R); 
  printf("total\n"); 
  printMDA( (void *) var, MDA_DOUBLE, 3, J, H, R); 
  */
  
  // CV 
  cvInit = cv_calcCV( NULL, N, K, H, R, J, domain, var, Nh, nh, total, locationAdj, scaleAdj, (size_t) *fpc ); 
  
  /* for debug */
  //printf("cvInit\n"); 
  //printMDA( (void *) cvInit, MDA_DOUBLE, 1, J); 

   
  // get change in allocation 
  alloc_sampleSizeChange (
      cvInit, N, K, H, R, J, domain, var, Nh, nh, total, locationAdj, scaleAdj, target, penalty, p, 0, iter, a, (size_t) * fpc);
  
  //printMDA( (void *) penalty, MDA_DOUBLE, 1, J); 
  

  // clean up 
  deleteMDA( (void *) var, 3, K, H); 
  deleteMDA( (void *) total, 2, J); 
  deleteMDA( (void *) cvInit, 1); 

  if( locationAdjDouble[0] >= 0 ) {
    deleteMDA( (void * ) locationAdj, 3, K, H); 
  }
  if( scaleAdjDouble[0] >= 0 ) {
    deleteMDA( (void *) scaleAdj, 3, K,H); 
  }

  return;
}





