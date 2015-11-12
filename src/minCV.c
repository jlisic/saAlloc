#include "minCV.h"

/************************** MINCV PROBLEM FUNCTIONS ***********************/
void minCV_init (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t K,              /* number of variables            */
            size_t N               /* number of elements within a state       */
            ) { 

  size_t i,j;
  minCV_adminStructPtr a;

  
  /* variable used to sum up the total over each stratum */ 

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 

    Q[0] <- 
      cv_objectiveFunction( 
      a->cv,
      N,
      a->K,
      a->H,
      a->R,
      a->J,
      a->domain,
      a->var,
      a->Nh,
      a->nh,
      a->Total,
      a->locationAdj_mu,
      a->scaleAdj_mu,
      a->T,
      a->penalty,
      a->p,
      a->preserveSatisfied  
    );
  

  return;
}
  
/* program flow */
/*
 * within SA
 *
 * costChange
 *
 *
 */









/* helper function  */
double maxExclude( double * x, size_t j, size_t H ) {

  size_t i;
  double maxx = x[0];

  for( i = 0; i < H; i++) {
    if( i == j) continue;
    if( maxx < x[i] ) maxx = x[i];  
  }

  return( maxx );
}


/* packSubstrata is a function that build all required data sets
 * and converts the output to the adminSubstrataStructPtr to
 * a void pointer to pass to sa. 
 */
void * minCV_pack( 
  double * x,               // A link to the Data Set 
  size_t * I,                 // assignments
  size_t N,                 // Number of strata 
  size_t K,                 // Number of strata 
  size_t H,                 // Number of strata 
  size_t R,                 // The number of observations of a PSU                    
  size_t J,                 // The number of targeted characteristics                    
  double * nh,              // sample size
  double * T,               // target CV 
  double * locationAdj,
  double * scaleAdj,
  double temp,              // temp denominator for cooling function
  double * prob,            // vector of sampling weights 
  double * probMatrix,      // matrix of sampling weights , one weight for each stratum 
  double * penalty,         // penalty coefficient for objective function 
  double p,                 // exponent for penalty  
  size_t iterSampleSize,    // number of iterations to optimize the sample size
  size_t preserveSatisfied, //    1 - do not go above a prior met constraint, 
                            //    0 - allow moving above a prior met constraint 
  int * domain              // domain HxKxJ data structure that will be turned into an MDA
  ) {

  size_t i,h,k,j,r;
   
  // allocate struct 
  minCV_adminStructPtr packedStruct = malloc( sizeof( minCV_adminStruct ) ); 
  
  //  target CV
  packedStruct->T = T;
  
  // H is the number of labels 
  packedStruct->H = H;
  
  // R is the number of observations of each PSU 
  packedStruct->R = R;
  
  // J is the number of domains 
  packedStruct->J = J;
  
  // K is the number of variables 
  packedStruct->K = K;
  
  // p is the exponent for the penalty function 
  packedStruct->p = p;

  // penalty is the coefficient of the penalty function
  packedStruct->penalty = penalty;
  
  // temperature  
  packedStruct->temp = temp;

  // population size by strata 
  packedStruct->Nh = calloc( H, sizeof(size_t));
  packedStruct->candidate_Nh = calloc( H, sizeof(size_t));
  for(i=0; i < N; i++) {
  //  printf("strata size [%d][%d]: %d from %d \n", (int) H, (int) i, (int) packedStruct->Nh[I[i]], (int) I[i]); 
    packedStruct->Nh[ I[i] ]++; 
  }
  for(i=0; i < H; i++) {
   // printf("strata size [%d]: %d\n", (int) i, (int) packedStruct->Nh[i]); 
    packedStruct->candidate_Nh[i] = packedStruct->Nh[i];  
  }

  // created to take care of the convert issue 
  packedStruct->prob       = prob;       // prob is the selection prob per obs 
  packedStruct->probMatrix = probMatrix; // probMatrix is the selection prob per obsxstratum 

  //sums fo prob
  packedStruct->totalStrataProb = calloc( H, sizeof(double));
  packedStruct->totalProb = 0;
  for(i=0; i < N; i++) packedStruct->totalStrataProb[ I[i] ] += prob[i];
  for(i=0; i < H; i++) packedStruct->totalProb += packedStruct->totalStrataProb[i]; 


  // sample size
  packedStruct->nh = nh;
  packedStruct->candidate_nh = calloc(H, sizeof(double));
  for(i=0; i<H; i++) {
    packedStruct->nh[i] = nh[i];
    packedStruct->candidate_nh[i] = nh[i];
  }
    
  // x - attach the data 
  packedStruct->x = x;
  
  // sample size iteration 
  packedStruct->iterSampleSize = iterSampleSize;

  // preserve satisfied constraints 
  // preserveSatisfied;     1 - do not go above a prior met constraint, 
  //                        0 - allow moving above a prior met constraint 
  packedStruct->preserveSatisfied = preserveSatisfied;
  
  // Domain Conversion 
  packedStruct->domain = (size_t ***) createMDA( MDA_SIZE_T, 3, H, K, J); 
  for( h=0; h < H; h++) 
    for( k=0; k < K; k++) 
      for( j=0; j < J; j++) 
        packedStruct->domain[h][k][j] = (size_t) domain[ h*K*J + k*J + j ];
  
  // location adjustement
  if( locationAdj[0] >= 0 ) { 
    packedStruct->locationAdj_mu = cv_createMeanMatrix(I, N, K, H, R, locationAdj, packedStruct->Nh); 
    packedStruct->locationAdj = locationAdj;

    packedStruct->candidate_locationAdj_mu = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
    for( k=0; k < K; k++) 
      for( h=0; h < H; h++) 
        for( r=0; r < R; r++) 
          packedStruct->candidate_locationAdj_mu[k][h][r] = packedStruct->locationAdj_mu[k][h][r];
  } else {
    packedStruct->locationAdj_mu = NULL;
    packedStruct->candidate_locationAdj_mu = NULL;
    packedStruct->locationAdj = NULL;
  }
  
  // scale adjustment
  if( scaleAdj[0] >= 0 ) { 
    packedStruct->scaleAdj_mu = cv_createMeanMatrix(I, N, K, H, R, scaleAdj, packedStruct->Nh); 
    packedStruct->scaleAdj = scaleAdj; 

    packedStruct->candidate_scaleAdj_mu = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
    for( k=0; k < K; k++) 
      for( h=0; h < H; h++) 
        for( r=0; r < R; r++) 
          packedStruct->candidate_scaleAdj_mu[k][h][r] = packedStruct->scaleAdj_mu[k][h][r];

  } else {
    packedStruct->scaleAdj_mu = NULL;
    packedStruct->candidate_scaleAdj_mu = NULL;
    packedStruct->scaleAdj = NULL;
  }

  // Mean
  packedStruct->mu = cv_createMeanMatrix(I, N, K, H, R, x, packedStruct->Nh); 
  packedStruct->candidate_mu = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
  for( k=0; k < K; k++) 
    for( h=0; h < H; h++) 
      for( r=0; r < R; r++) 
        packedStruct->candidate_mu[k][h][r] = packedStruct->mu[k][h][r];

  // Variance
  packedStruct->var = cv_createVarMatrix(I, N, K, H, R, x, 
    packedStruct->mu, packedStruct->Nh); 
  packedStruct->candidate_var = (double ***) createMDA( MDA_DOUBLE, 3, K, H, R); 
  for( k=0; k < K; k++) 
    for( h=0; h < H; h++) 
      for( r=0; r < R; r++) 
        packedStruct->candidate_var[k][h][r] = packedStruct->var[k][h][r];

  // Total 
  packedStruct->Total =  cv_createTotalMatrix( N, K, H, R, J, packedStruct->domain, packedStruct->mu, packedStruct->Nh);

  // CV
  packedStruct->cv = cv_calcCV( NULL, N, K, H, R, J, 
    packedStruct->domain, 
    packedStruct->var, 
    packedStruct->Nh, 
    packedStruct->nh, 
    packedStruct->Total, 
    packedStruct->locationAdj_mu, 
    packedStruct->scaleAdj_mu 
  ); 
  packedStruct->candidate_cv = calloc(J, sizeof( double ) );
  for(j = 0; j < J; j++) packedStruct->candidate_cv[j] = packedStruct->cv[j];

  return( (void *) packedStruct );
}




/* clean up for the substrata administrative data */
void minCV_delete( minCV_adminStructPtr A, size_t K, size_t N ) {

  // allocate struct 
  minCV_adminStructPtr packedStruct =  (minCV_adminStructPtr) A;
  
  // H is the number of labels 
  size_t H = packedStruct->H;
  
  // R is the number of observations of each PSU 
  size_t R = packedStruct->R;
  
  // J is the number of domains 
  size_t J = packedStruct->J;
 
 // K handled 
//  printf("delete\n"); 
 
//  printf("Nh\n"); 
  deleteMDA( (void *) packedStruct->Nh, 1) ;
//  printf("can_Nh\n"); 
  deleteMDA( (void *) packedStruct->candidate_Nh, 1) ;
//  printf("can_nh\n"); 
  deleteMDA( (void *) packedStruct->candidate_nh, 1) ;
  
//  printf("domain\n"); 
  deleteMDA( (void *) packedStruct->domain, 3, H, K) ;
 
  if( packedStruct->locationAdj_mu != NULL ) { 
//    printf("locationAdj\n"); 
    deleteMDA( (void *) packedStruct->locationAdj_mu, 3, K, H) ;
    deleteMDA( (void *) packedStruct->candidate_locationAdj_mu, 3, K, H) ;
  }
  if( packedStruct->scaleAdj_mu != NULL ) { 
//    printf("scaleAdj\n"); 
    deleteMDA( (void *) packedStruct->scaleAdj_mu, 3, K, H) ;
    deleteMDA( (void *) packedStruct->candidate_scaleAdj_mu, 3, K, H) ;
  }

//  printf("mu\n"); 
  deleteMDA( (void *) packedStruct->mu, 3, K, H) ;
//  printf("candidate mu\n"); 
  deleteMDA( (void *) packedStruct->candidate_mu, 3, K, H) ;

//  printf("var\n"); 
  deleteMDA( (void *) packedStruct->var, 3, K, H) ;
//  printf("candidate var\n"); 
  deleteMDA( (void *) packedStruct->candidate_var, 3, K, H) ;

//  printf("Total\n"); 
  deleteMDA( (void *) packedStruct->Total, 2, J) ;

//  printf("cv\n"); 
  deleteMDA( (void *) packedStruct->cv, 1) ;
  
//  printf("candidate cv\n"); 
  deleteMDA( (void *) packedStruct->candidate_cv, 1) ;

//  printf("total strata probability\n"); 
  deleteMDA( (void *) packedStruct->totalStrataProb, 1) ;

  // finally free the struct
  free( packedStruct);
  packedStruct = NULL;

}



/* function to check status */
void minCV_diag( 
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            void * A,   /* administrative data */
            size_t K,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
    ) {

  size_t h = 0;
 
  // allocate struct 
  minCV_adminStructPtr packedStruct =  (minCV_adminStructPtr) A;


  if( N < 30 ) {
  Rprintf("I:\n");
  for(h = 0; h < N; h++) Rprintf("[%d] %d\n", (int) h, (int) I[h] ); 
  }


  
  // H is the number of labels 
  size_t H = packedStruct->H;
  
  // R is the number of observations of each PSU 
  size_t R = packedStruct->R;
  
  // J is the number of domains 
  size_t J = packedStruct->J;
 
  // K handled 
  Rprintf("H = %d, K = %d, J = %d, R = %d, N = %d\n", (int) H, (int) K, (int) J, (int) R, (int) N);
  Rprintf("iterSampleSize = %d\n", packedStruct->iterSampleSize); 
  Rprintf("preserveSatisfied = %d\n", packedStruct->preserveSatisfied); 
  Rprintf("p = %f, temp = %f\n", (float) packedStruct->p, (float) packedStruct->temp);

  // created to take care of the convert issue 
  //packedStruct->prob       = prob;       // prob is the selection prob per obs 
  //packedStruct->probMatrix = probMatrix; // probMatrix is the selection prob per obsxstratum 
  
  Rprintf("Q\n");
  printMDA( (void *) Q, MDA_DOUBLE, 1, 1); 
  
  Rprintf("T\n"); //targetCV
  printMDA( (void *) packedStruct->T, MDA_DOUBLE, 1, J); 
  
  Rprintf("Total Strata Probability\n"); //targetCV
  printMDA( (void *) packedStruct->totalStrataProb, MDA_DOUBLE, 1, H); 
  Rprintf("Total Probability = %4.2f\n", packedStruct->totalProb); //targetCV
  
  Rprintf("cv\n");
  printMDA( (void *) packedStruct->cv, MDA_DOUBLE, 1, J); 
  Rprintf("candidate cv\n");
  printMDA( (void *) packedStruct->candidate_cv, MDA_DOUBLE, 1, J); 

  Rprintf("Nh\n");
  printMDA( (void *) packedStruct->Nh, MDA_SIZE_T, 1, H); 

  Rprintf("candidate Nh\n");
  printMDA( (void *) packedStruct->candidate_Nh, MDA_SIZE_T, 1, H); 
  
  Rprintf("nh\n");
  printMDA( (void *) packedStruct->nh, MDA_DOUBLE, 1, H); 
  Rprintf("candidate nh\n");
  printMDA( (void *) packedStruct->candidate_nh, MDA_DOUBLE, 1, H); 

  Rprintf("domain\n");
  printMDA( (void *) packedStruct->domain, MDA_SIZE_T, 3, H, K, J); 
 
  Rprintf("penalty\n");
  printMDA( (void *) packedStruct->penalty, MDA_DOUBLE, 1, J); 


  // location adjustement
  if( packedStruct->locationAdj != NULL ) { 
    Rprintf("Location Adjustment\n");
    printMDA( (void *) packedStruct->locationAdj_mu, MDA_DOUBLE, 3, K, H, R); 
  } else {
    Rprintf("No Location Adjustment\n");
  }
  
  // scale adjustment
  if( packedStruct->scaleAdj != NULL ) { 
    Rprintf("Scale Adjustment\n");
    printMDA( (void *) packedStruct->scaleAdj_mu, MDA_DOUBLE, 3, K, H, R ); 
  } else {
    Rprintf("No Scale Adjustment\n");
  }

  // Mean
  Rprintf("Mean\n");
  printMDA( (void *) packedStruct->mu, MDA_DOUBLE, 3, K, H, R); 
  Rprintf("Candidate Mean\n");
  printMDA( (void *) packedStruct->candidate_mu, MDA_DOUBLE, 3, K, H, R); 

  // Variance
  Rprintf("Variance\n");
  printMDA( (void *) packedStruct->var, MDA_DOUBLE, 3, K, H, R); 
  Rprintf("Candidate Variance\n");
  printMDA( (void *) packedStruct->candidate_var, MDA_DOUBLE, 3, K, H, R); 

  // Total 
  Rprintf("Total\n");
  printMDA( (void *) packedStruct->Total, MDA_DOUBLE, 2, J,  R); 

  return;
}
  


/* get a new random state */
size_t minCV_randomState (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * IFin,         /* new state                               */
            double * QFin,         /* new cost                                */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t K,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
    ) { 

  minCV_adminStructPtr a;

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  double * nh = a->nh;
  size_t * Nh = a->Nh;
  size_t i =N+1;
  size_t H = a->H;
  size_t trys = 0;  /* number of times to try */

  size_t h ;

//  for(h=0;h<N;h++) printf("%4.2f ", a->prob[h] );
//  printf("\n");

  while( i >= N ) {

    /* generate possible move */
    i = minCV_getIndex( 
        a->prob, 
        a->totalProb, 
        a->totalStrataProb
      );
 
    /* check if there are units that I can move */ 
    if( Nh[ I[i] ] > 2 ) {
      
      a->Hj = minCV_getMoveStrata( i, I, a->prob, a->probMatrix, N, H);
  
      if( a->Hj == N+1 ) i = N+1;
  
  
    } else {
      i = N+1;
    }

    /* fail on lack of matches */
    if( trys >= 10000 ) {
      Rprintf(" 10,000 failures on finding a possible move, giving up\n"); 
      break;
    }

    trys++;
  }

  // record Hi
  a->Hi = I[i];

//  printf("Moving %d [%d] to [%d]\n", (int) i, (int) a->Hi, (int) a->Hj);

  /* return index */
  return(i);
}



/* function to select an index */
size_t minCV_getIndex( double * prob, double totalProbability, double * totalStrataProbability ) {

  size_t index = 0;
//  size_t strataIndex = 0;
  double total = 0;

//  size_t h = 0;

  /* select a value uniform across the sum of probabilities */
  double searchProbability = runif(0,totalProbability); 
  
//  printf("totalTest = %4.2f\n", totalStrataProbability[0]); 

//  printf("searchProbability = %4.2f, totalProbability = %4.2f\n",
//    searchProbability, totalProbability);
//  for( h = 0; h < 3; h++)  printf("%d: %4.2f\n", (int) h, totalStrataProbability[h]); 

   
  /* find the selected strata */ 
//  total = totalStrataProbability[0]; 
//  while( total < searchProbability ) {
//    strataIndex++;
//    total += totalStrataProbability[strataIndex];
//  }
    
  /* adjust search probability */ 
  //searchProbability = searchProbability - total + totalStrataProbability[strataIndex];

  total = prob[0]; 
  while( total < searchProbability ){ 
      index++;
      total += prob[index];
  }

  return( index );
}



/* get a strata for i to move to */
size_t minCV_getMoveStrata( 
    size_t i, 
    size_t * I, 
    double * prob, 
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
    //if( Hj != Hi ) totalProb += probMatrix[i*H + Hj];
    totalProb += probMatrix[i*H + Hj];
  }
  /* if the total probability remaininig is 0 quit */
  if( totalProb == 0 ) return( N + 1 );   

  /* randomly select the strata proportional to its weight */
  totalProb = runif(0,totalProb); 
  total = 0;

  /* find the strata that corresponds to the selected sampling weight */ 
  for( Hj=0; Hj < H; Hj++) {
    
    /* if the strata is in a different from i's strata then add total */ 
    //if( Hj != Hi ) 
    total += probMatrix[i*H + Hj];

    if( total > totalProb) break; /* unless all the probability is 0 then units are sampled on an 'open' interval */ 
  }
  
  return( Hj ) ;
}



/* This is a function to check the change in substrata cost                   */
double minCV_costChange (
            size_t i,              /* new Index                               */
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * IFin,         /* new state                               */
            double * QFin,         /* new cost                                */
            double * D,            /* distance matrix     (NOT USED)          */
            void * A,              /* administrative data                     */
            size_t K,              /* number of variables                     */
            size_t N               /* number of PSUs                          */
    ) { 

  /* cast A back to something useable */
  minCV_adminStructPtr a = (minCV_adminStructPtr) A; 
 
  size_t J = a->J;  // number of domains
  size_t H = a->H;  // number of strata
  size_t R = a->R;  // number of repeat observations for each PSU

  double delta = INFINITY; 
  double obj, candidate_obj;



  // apply the variance change 
  // note that we only update the cv if there is a change in strata 
  if( a->Hi != a->Hj ) 
  cv_updateMatrix2( I, N, K, H, R, 
      a->x, 
      a->locationAdj,
      a->scaleAdj,
      a->candidate_mu, 
      a->candidate_locationAdj_mu,
      a->candidate_scaleAdj_mu,
      a->candidate_var, 
      a->Nh, 
      i, 
      a->Hj ); 

  // update Nh
  a->candidate_Nh[I[i]]--;
  a->candidate_Nh[a->Hj]++;

  // calculate the CV 
  cv_calcCV( 
      a->candidate_cv, 
      N, K, H, R, J, 
      a->domain, 
      a->candidate_var, 
      a->candidate_Nh, 
      a->candidate_nh, 
      a->Total, 
      a->candidate_locationAdj_mu, 
      a->candidate_scaleAdj_mu
    ); 

  // get change in allocation 
  alloc_sampleSizeChange(
      a->candidate_cv, 
      N, K, H, R, J, 
      a->domain, 
      a->candidate_var, 
      a->candidate_Nh, 
      a->candidate_nh, 
      a->Total, 
      a->candidate_locationAdj_mu, 
      a->candidate_scaleAdj_mu, 
      a->T, 
      a->penalty, 
      a->p, 
      a->preserveSatisfied, 
      a->iterSampleSize, 
      NULL
    );
  

  //if( nhSum != nhSumCan ) 
  //printf(" nhSum = %f , nhSumCan = %f\n", nhSum, nhSumCan);


  delta =  cv_objectiveFunctionCompare( 
    a->candidate_cv,         //if cv = NULL create it, otherwise handle in place 
    a->cv,   
    &obj,
    &candidate_obj,
    N, K, H, R, J, 
    NULL,     // only used if evaluateOnly != 1
    NULL,        // only used if evaluateOnly != 1 
    NULL,
    NULL,
    NULL,       // only used if evaluateOnly != 1
    NULL, // only used if evaluateOnly != 1
    NULL, // only used if evaluateOnly != 1
    a->T,
    a->penalty,
    a->p,
    1, // option to not construct CV, under this condition CV cannot be null 
    a->preserveSatisfied
  );
    
//  printf("obj: %f objPrior %f\n", obj, candidate_obj); 
  
  return(delta);
}







/* update from change in i function */
/* within this step the following variables are updated
 * I - assignment index
 * C - contribution 3D arrray
 * V - variance
 * Nh - number of exchangable units in each stratum
 * NhSize - number of total units in each stratum
 * NhAcres - number of acres in each stratum
 * Total - Total acres in each stratum
 */

void minCV_update (
            size_t accept, /* 1 if accepted 0 if not, 2 handles additional functions */
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            size_t * IFin, /* new state, destructive */
            double * QFin, /* new cost, destructive */
            double * D, /* distance matrix */
            void * A,   /* administrative data */
            size_t K,  /* number of distance matricies */
            size_t N,   /* number of elements within a state */
            double * costChange /* cost Change */
            ) { 

  size_t k, j, r,h;  
  minCV_adminStructPtr a = (minCV_adminStructPtr) A; 
  size_t R = a->R;
  size_t H = a->H;
  size_t J = a->J;
  size_t Hi = a->Hi;
  size_t Hj = a->Hj; 
 

  double nhSumStop, nhSumStart;
  for(h = 0, nhSumStart = 0; h < H; h++) nhSumStart += a->nh[h]; 


  
  /* restore prior state candidates */
  if( accept == 0) {

    /* change back mu, var and adjustments */
    for( k = 0; k < K; k++) {
      for( r = 0; r < R; r++) {
        a->candidate_var[k][Hj][r] = a->var[k][Hj][r];
        a->candidate_mu[k][Hj][r] = a->mu[k][Hj][r];
        
        a->candidate_var[k][Hi][r] = a->var[k][Hi][r];
        a->candidate_mu[k][Hi][r] = a->mu[k][Hi][r];
      } 
      if( a->scaleAdj != NULL) { 
        for( r = 0; r < R; r++) {
          a->candidate_scaleAdj_mu[k][Hj][r] = a->scaleAdj_mu[k][Hj][r];
          a->candidate_scaleAdj_mu[k][Hi][r] = a->scaleAdj_mu[k][Hi][r];
        }
      }
      if( a->locationAdj != NULL) { 
        for( r = 0; r < R; r++) {
          a->candidate_locationAdj_mu[k][Hj][r] = a->locationAdj_mu[k][Hj][r];
          a->candidate_locationAdj_mu[k][Hi][r] = a->locationAdj_mu[k][Hi][r];
        }
      }
    } 

    return;
  }


  /* record objective function and sample sizes for diagnostics */
  if ( accept == 2) { 

    /* record what strata the selected PSU is moving to */
    costChange[4] = Hi;
    costChange[5] = Hj;

    /* add on sample size */
    for( h = 0; h < H; h++) costChange[6+h] = a->candidate_nh[h]; 
    /* add on current obj function */
    for( j = 0; j < J; j++) {
      costChange[6 + H +j] = a->candidate_cv[j]; 
    }
      
    /* change back nh & Nh */
    for( h = 0; h < H; h++) {
      a->candidate_nh[h] = a->nh[h]; 
      a->candidate_Nh[h] = a->Nh[h]; 
    }

//    for( j = 0; j < H + J + 6; j++) printf("%4.2f ", costChange[j]);
//    printf("\n");

    return;
  }
  
  /* optimize sample size */
  if ( accept == 3) { 
    //printf("accept == 3\n");
    //minCV_sampleSizeChange ( A, Q, R, dN, N, a->sampleIter);
    return;
  }
    
  //printf("accept == 1, i=%d, Hj=%d, I[i]%d\n", (int) i, (int) Hj, (int) I[i]);

  // check if there is no change 
  if( a->Hi == a->Hj) { 
    for(h=0; h < H; h++) if( a->nh[h] != a->candidate_nh[h] ) break;
    if( h == H) return; 
  }
  
  /* update the strata assignment */
  I[i] = Hj;

  /* adjust strata probabilities */
  // remove old prob
  a->totalStrataProb[Hi] -= a->prob[i];
  a->totalProb -= a->prob[i];
  
  // get new prob based on new stratum
  a->prob[i] = maxExclude( &( a->probMatrix[H * i] ), Hj, H );

  /* add on new prob */
  a->totalStrataProb[Hj] += a->prob[i];
  a->totalProb += a->prob[i];

  /* change back nh */
  for( h=0; h < H; h++) {
    a->nh[h] = a->candidate_nh[h]; 
    a->Nh[h] = a->candidate_Nh[h]; 
  }



  /* change back mu, var and adjustments */
  for( k = 0; k < K; k++) {
    /* change back cv */
    a->cv[k] = a->candidate_cv[k]; 

    for( r = 0; r < R; r++) {
      a->var[k][Hj][r] = a->candidate_var[k][Hj][r];
      a->mu[k][Hj][r] = a->candidate_mu[k][Hj][r];
      
      a->var[k][Hi][r] = a->candidate_var[k][Hi][r];
      a->mu[k][Hi][r] = a->candidate_mu[k][Hi][r];
    } 
    if( a->scaleAdj != NULL) { 
      for( r = 0; r < R; r++) {
        a->scaleAdj_mu[k][Hj][r] = a->candidate_scaleAdj_mu[k][Hj][r];
        a->scaleAdj_mu[k][Hi][r] = a->candidate_scaleAdj_mu[k][Hi][r];
      }
    }
    if( a->locationAdj != NULL) { 
      for( r = 0; r < R; r++) {
        a->locationAdj_mu[k][Hj][r] = a->candidate_locationAdj_mu[k][Hj][r];
        a->locationAdj_mu[k][Hi][r] = a->candidate_locationAdj_mu[k][Hi][r];
      }
    }
  } 
    
  for(h = 0, nhSumStop = 0; h < H; h++) nhSumStop += a->nh[h]; 
  if( nhSumStop != nhSumStart ) printf("update: nhSumStop = %f , nhSumStart = %f\n", nhSumStop, nhSumStart);
 
  return;
}



/* cooling schedule                                                           */
double minCV_cool ( 
            size_t iter,   /* iteration */
            double diff, /* cost change */ 
            size_t i,   /* new Index */ 
            size_t * I, /* current state */ 
            double * Q, /* current cost */
            double * D, /* distance matrix */
            void * A, /* administrative data */ 
            size_t dN,   /* number of distance matricies */
            size_t N    /* number of elements within a state */
            )  { 
  double temp;
  minCV_adminStructPtr a;
  
  a = (minCV_adminStructPtr) A;
  temp = a->temp;

  /*return( exp( (log(1+iter)) * diff /  (-1* temp) ) ); */ 
  return( exp( (1+iter) * diff /  (-1* temp) ) ) ; 
}  


