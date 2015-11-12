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

  size_t Hi = I[i];
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
    double * a
    ) { 


  double * cv = calloc( J, sizeof(double) );  // a place to store changes
  double * test_nh = calloc( H, sizeof(double) );  // a place to store changes
  size_t i,h,j;         // iteration
  double minSampleSize = 2.0; // minimum sample size
  double min_delta, delta;
  size_t Hi, Hj, optHj;

  /* sanity check */
  double nhSumStart;
  double nhSumStop;

  for(h = 0, nhSumStart = 0; h < H; h++) nhSumStart += nh[h]; 
 

  // if there is nothing to do, do nothing 
  if( iter == 0) return;
 
  //for(h = 0; h < H; h++) nhSumStart += nh[h]; 

  // copy values over
  for(h = 0; h < H; h++) test_nh[h] = nh[h]; 
  
  
  // iterate a fixed number of times 
  for( i = 0 ; i < iter; i++ ) {
    
    // 1.0 randomly select a strata to move from 
    Hi = SA_GETINDEX(H);
    optHj = H;


    // only proceed if stratum Hi can be made smaller 
    if( test_nh[Hi] < minSampleSize + 1 ) {
      continue;
    }
    
    test_nh[Hi] -= 1.0; // perform first iteration
      
    // 2.0 calculate objective function change for moving to each strata 
    min_delta = 0.0;

    for( Hj = 0; Hj < H; Hj++ ) {
 
      // if the exchange would make the sample size too big, we don't do it 
      if( test_nh[Hj] + 1 > Nh[Hj] ) {

        continue;
      }

      // do not run over the selected state 
      if( Hj != Hi) {
        
        test_nh[Hj] += 1.0;
      
//        printf("test = %f\n", var[0][0][0]);
        /* calculate cv and write it to cv */ 
        // gcc is crazy and likes to delete this bit
        //void __attribute__((optimize("O0"))) foo(unsigned char data) {
        
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
          preserveSatisfied); 

        //}

        /*
        printf("%d delta = %f\n",(int) i, delta);
        printf("cv:  ");
        for(j = 0; j < J; j++) printf("%f,", cv[j]); 
        printf("\n");
        printf("cv - init:  ");
        for(j = 0; j < J; j++) printf("%f,", cvInit[j]); 
        printf("\n");
          
        //delta = 1;
        */ 
        
        test_nh[Hj] -= 1.0;

        //printf("%d, %f\n", (int) i, delta);

        // check if we are doing better 
        if( delta <= min_delta ) {
          optHj = Hj;
          min_delta = delta;
          for(j = 0; j < J; j++) cvInit[j] = cv[j]; 
        } 
      }
    }
    
    // 3.0 if there are reductions in objective function make a move to minimize CV 
  
    if( optHj < H ) {  // check if the final result minimized the objective function for all H 
      test_nh[optHj] += 1.0;
    } else {
      test_nh[Hi] += 1.0; // undo our change 
    }

   // sanity check 
    for(h = 0, nhSumStop = 0; h < H; h++) nhSumStop += test_nh[h]; 
    if( nhSumStop != nhSumStart ) printf("alloc: nhSumStop = %f , nhSumStart = %f\n", nhSumStop, nhSumStart);

    // return min_delta
    if( a != NULL ) a[i] = min_delta;
  }
    

  // copy values over
  for(h = 0; h < H; h++) nh[h] = test_nh[h]; 
  free(test_nh);
  free(cv);

  return;
}







#ifdef TESTALLOC

int main() {

  double probMatrix[]  = {
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,

    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,

    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2
  }; 

  double x[]  = {
    1.0, 2.0, 3.0, 
    4.0, 4.0, 5.0, 
    5.0, 5.0, 1.0, 
    2.0, 3.0, 4.0, 
    2.0, 3.0, 4.0, 

    4.0, 5.0, 5.0, 
    5.0, 1.0, 2.0, 
    3.0, 4.0, 4.0, 
    5.0, 5.0, 5.0, 
    5.0, 5.0, 5.0, 

    0.0, 2.0, 3.0, 
    4.0, 5.0, 5.0, 
    5.0, 5.0, 0.0, 
    2.0, 3.0, 4.0, 
    2.0, 3.0, 4.0, 
   

    4.0, 6.0, 5.0, 
    5.0, 0.0, 2.0, 
    3.0, 4.0, 4.0, 
    3.0, 4.0, 4.0, 
    5.0, 7.0, 5.0,

    0.0, 2.0, 3.0, 
    4.0, 5.0, 5.0, 
    5.0, 5.0, 0.0, 
    2.0, 3.0, 4.0, 
    2.0, 3.0, 4.0, 
    
    4.0, 6.0, 5.0, 
    5.0, 0.0, 2.0, 
    3.0, 4.0, 4.0, 
    3.0, 4.0, 4.0, 
    5.0, 7.0, 5.0 
  };

  size_t I[]  = {  0,   0,   0,   0,  0,  1,  1,   1,   1,   1, 2, 2, 2, 2, 2 };
  size_t Nh[] = {  5,   5, 5 };
  double nh[] = {  3,   3, 2 };
  double target[] = { 1.0, 1.1, 1.2, 1.4 };

  size_t K = 2;
  size_t N  = 15;
  size_t R = 3;
  size_t H  = 3;
  size_t J  = 4;
  size_t i;

  size_t domain00[] = { 1, 0, 1, 0 };
  size_t domain01[] = { 0, 1, 0, 0 };
  
  size_t domain10[] = { 1, 0, 0, 0 };
  size_t domain11[] = { 0, 1, 0, 0 };
 
  size_t domain20[] = { 1, 0, 1, 0 };
  size_t domain21[] = { 0, 1, 0, 1 };

  size_t * domain0[] = { domain00, domain01 };
  size_t * domain1[] = { domain10, domain11 };
  size_t * domain2[] = { domain20, domain21 };

  size_t ** domain[] = {domain0, domain1, domain2};


  double penalty[] = {10, 10, 10, 10};
  double ** total;
  double * cvInit;
  double p = 1;
   
  size_t iter = 5;


  // Mean
  Rprintf("Mean\n");
  double *** mu = cv_createMeanMatrix(I, N, K, H, R, x, Nh); 
  printMDA( (void *) mu, MDA_DOUBLE, 3, K, H, R); 

  // Variance
  Rprintf("Variance\n");
  double *** var = cv_createVarMatrix(I, N, K, H, R,  x, mu, Nh); 
  printMDA( (void *) var, MDA_DOUBLE, 3, K, H, R); 

  // Total 
  total =  cv_createTotalMatrix( N, K, H, R, J, (size_t ***) domain, mu, Nh);
  Rprintf("Total\n");
  printMDA( (void *) total, MDA_DOUBLE, 2, J,  R); 
  
  // CV
  Rprintf("CV\n");
  cvInit = cv_calcCV( NULL, N, K, H, R, J, domain, var, Nh, nh, total, NULL, NULL ); 
  printMDA( (void *) cvInit, MDA_DOUBLE, 1, J); 
  
  // Sample Size 
  Rprintf("nh - init\n");
  printMDA( (void *) nh, MDA_DOUBLE, 1, H); 
 
  // get change in allocation 
  alloc_sampleSizeChange (
      cvInit, N, K, H, R, J, domain, var, Nh, nh, total, NULL, NULL, target, penalty, p, 0, iter);
  
  // CV
  Rprintf("CV -final\n");
  printMDA( (void *) cvInit, MDA_DOUBLE, 1, J); 

  // Sample Size 
  Rprintf("nh - final\n");
  printMDA( (void *) nh, MDA_DOUBLE, 1, H); 


  // clean up 
  printf(" Deleteing Mean\n");
  deleteMDA( (void * ) mu, 3, K, H); 
  printf(" Deleteing Var\n");
  deleteMDA( (void *) var, 3, K, H); 
  printf(" Deleteing Total\n");
  deleteMDA( (void *) total, 2, K); 
  printf(" Deleteing CV\n");
  deleteMDA( (void *) cvInit, 1, J); 


  printf("Finished\n");
  return( 0 ) ;
}

#endif





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
  double * a
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

  
  /* create domain MDA */
  /* H x K x J */

  /*
  size_t domain00[] = { 1, 0, 1, 0 };
  size_t domain01[] = { 0, 1, 0, 0 };
  
  size_t domain10[] = { 1, 0, 0, 0 };
  size_t domain11[] = { 0, 1, 0, 0 };
 
  size_t domain20[] = { 1, 0, 1, 0 };
  size_t domain21[] = { 0, 1, 0, 1 };

  size_t * domain0[] = { domain00, domain01 };
  size_t * domain1[] = { domain10, domain11 };
  size_t * domain2[] = { domain20, domain21 };
  */
 
  double *** locationAdj;
  double *** scaleAdj; 

  // Mean
  Rprintf("Domain\n");
  size_t *** domain = (size_t ***) createMDA( MDA_SIZE_T, 3, H, K, J); 
  for( h=0; h < H; h++)
    for( k=0; k < K; k++)
      for( j=0; j < J; j++) domain[h][k][j] = (size_t) domainInt[ h*K*J + k*J + j ];
  printMDA( (void *) domain , MDA_SIZE_T, 3,  H, K, J); 

  
  // location adjustement
  if( locationAdjDouble[0] >= 0 ) { 
    Rprintf("Location Adjustment\n");
    locationAdj = cv_createMeanMatrix(I, N, K, H, R, x, Nh); 
    printMDA( (void *) locationAdj, MDA_DOUBLE, 3, K, H, R); 
  } else {
    locationAdj = NULL;
  }
  
  // scale adjustmentf
  if( scaleAdjDouble[0] >= 0 ) { 
    Rprintf("Scale Adjustment\n");
    scaleAdj = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
    for( k=0; k < K; k++) 
      for( h=0; h < H; h++) 
        for( r=0; r < R; r++) scaleAdj[k][h][r] = scaleAdjDouble[k*H*R + h*R + r];
    printMDA( (void *) scaleAdj, MDA_DOUBLE, 3, K, H, R ); 
  } else {
    scaleAdj = NULL;
  }

  // Mean
  Rprintf("Mean\n");
  double *** mu = cv_createMeanMatrix(I, N, K, H, R, x, Nh); 
  printMDA( (void *) mu, MDA_DOUBLE, 3, K, H, R); 

  // Variance
  Rprintf("Variance\n");
  double *** var = cv_createVarMatrix(I, N, K, H, R,  x, mu, Nh); 
  printMDA( (void *) var, MDA_DOUBLE, 3, K, H, R); 

  // Total 
  total =  cv_createTotalMatrix( N, K, H, R, J, (size_t ***) domain, mu, Nh);
  Rprintf("Total\n");
  printMDA( (void *) total, MDA_DOUBLE, 2, J,  R); 

  
  // CV
  Rprintf("CV\n");
  cvInit = cv_calcCV( NULL, N, K, H, R, J, domain, var, Nh, nh, total, locationAdj, scaleAdj ); 
  printMDA( (void *) cvInit, MDA_DOUBLE, 1, J); 
  
  // Sample Size 
  Rprintf("nh - init\n");
  printMDA( (void *) nh, MDA_DOUBLE, 1, H); 
 
  // get change in allocation 
  alloc_sampleSizeChange (
      cvInit, N, K, H, R, J, domain, var, Nh, nh, total, locationAdj, scaleAdj, target, penalty, p, 0, iter, a);
  
  // CV
  Rprintf("CV -final\n");
  printMDA( (void *) cvInit, MDA_DOUBLE, 1, J); 

  // Sample Size 
  Rprintf("nh - final\n");
  printMDA( (void *) nh, MDA_DOUBLE, 1, H); 


  // clean up 
  Rprintf(" Deleteing Mean\n");
  deleteMDA( (void * ) mu, 3, K, H); 
  Rprintf(" Deleteing Var\n");
  deleteMDA( (void *) var, 3, K, H); 
  Rprintf(" Deleteing Total\n");
  deleteMDA( (void *) total, 2, J); 
  Rprintf(" Deleteing CV\n");
  deleteMDA( (void *) cvInit, 1); 
  Rprintf(" Deleteing Nh\n");
  deleteMDA( (void *) Nh, 1); 
  Rprintf(" Deleteing I\n");
  deleteMDA( (void *) I, 1); 

  if( locationAdjDouble[0] >= 0 ) {
    Rprintf(" Deleteing Location Adj\n");
    deleteMDA( (void * ) locationAdj, 3, K, H); 
  }
  if( scaleAdjDouble[0] >= 0 ) {
    Rprintf(" Deleteing Scale Adj\n");
    deleteMDA( (void *) scaleAdj, 3, K,H); 
  }

  Rprintf("Finished\n");
  return;
}






