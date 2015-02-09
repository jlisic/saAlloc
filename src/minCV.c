#include "minCV.h"

/************************** MINCV PROBLEM FUNCTIONS ***********************/

/* This is a function to check initialize the structures and variables for
 * the sa proceedure                                                          */
void minCV_init (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
            ) { 

  size_t i,j;
  minCV_adminStructPtr a;

  
  /* variable used to sum up the total over each stratum */ 

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  double ** V = a->V;
  double ** RV = a->RV;
  double * T = a->T;
  double * Total = a->Total;
  size_t H = a->H;
  size_t * NhSize = a->NhSize;
  double * sampleSize = a->sampleSize;
  double * RSampleSize = a->RSampleSize;


  /* Q has been allocated but there is nothing in it yet */
  for(i=0; i < dN; i++) {
    Q[i] = 0;
    for(j=0; j < H; j++) {
      Q[i] += V[i][j] * NhSize[j]* NhSize[j]/( sampleSize[j] );

      /* copy over V to RV */ 
      RV[i][j] = V[i][j];
    }
    Q[i] = sqrt(Q[i]) / Total[i] - T[i];
 
  }

  /* copy over the initial sample sizes to a set of 'Reserved' Sample sizes */ 
  for( j=0; j < H; j++) {
    RSampleSize[j] = sampleSize[j]; 
  }

  return;
}


  

/* get a new random state */
size_t minCV_randomState (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * J,            /* new state                               */
            double * R,            /* new cost                                */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
    ) { 

  minCV_adminStructPtr a;


  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  double * prob = a->prob;
  double * probMatrix = a->probMatrix;
  double * sampleSize = a->sampleSize;
  size_t * NhSize = a->NhSize;
  double totalProbability = a->totalProbability;
  size_t i =N+1;
  size_t H = a->H;
  size_t trys = 0;  /* number of times to try */
  size_t d;

  /* first check if there is anything to do */
  for( d = 0; d < dN; d++) {
    if( Q[d] > 0 ) break;
  }

  /* if d = dN then all constraints are satisfied */
  if( d == dN) {
    Rprintf("All constraints satisfied, returning early\n");
    return(N+1);
  }


  while( i >= N ) {

    /* generate possible move */
    i = minCV_getIndex( prob, totalProbability );
 
    /* check if there are units that I can move */ 
    if( sampleSize[ I[i] ] < NhSize[ I[i] ] ) {
      
      a->Hj = minCV_getMoveStrata( i, I, prob, probMatrix, N, H);
  
      if( a->Hj == N+1 ) i = N+1; /* i = N+1 is used to signal to sa that no move can be made */
  
  
    } else {
      i = N+1;
    }

    /* fail on lack of matches */
    if( trys >= 1000 ) {
      Rprintf(" 1000 failures on finding a possible move, giving up\n"); 
      break;
    }

    trys++;
  }
  

  /* return index */
  return(i);
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

  /* get the total weight of all the strata except Hi */
  for( Hj=0; Hj < H; Hj++) 
    if( Hj != Hi ) totalProb += probMatrix[i*H + Hj];

  /* if the total probability remaininig is 0 quit */
  if( totalProb == 0 ) return( N + 1 );   

  /* randomly select the strata proportional to its weight */
  totalProb = runif(0,totalProb); 
  total = 0;

  /* find the strata that corresponds to the selected sampling weight */ 
  for( Hj=0; Hj < H; Hj++) {
    
    /* if the strata is in a different from i's strata then add total */ 
    if( Hj != Hi ) total += probMatrix[i*H + Hj];

    if( total >= totalProb) break;  
  }
  
  return( Hj ) ;
}



/* function to select an index */
size_t minCV_getIndex( double * prob, double totalProbability ) {

  size_t index = 0;
  double total = 0;
  double searchProbability = runif(0,totalProbability); /* select a value uniform across the sum of probabilities */
      
  total = prob[0]; 
  while( total < searchProbability ){ 
      index++;
      total += prob[index];
  }

  return( index );
}


/* This is a function to update the sample size                              */
void minCV_sampleSizeChange (
            void * A,              /* administrative data                     */
            double * R,
            size_t dN,             /* number of distance matricies            */
            size_t N,              /* number of elements within a state       */
            size_t iter
    ) { 

  minCV_adminStructPtr a;

  /* cast A back to something useable */
  a = (minCV_adminStructPtr) A; 
  
  size_t H = a->H;
  double * T = a->T;
  size_t * NhSize = a->NhSize;
  double ** RV = a->RV;
  double * Total = a->Total;
  double * sampleSize = a->sampleSize; /* sample Size */
  double * RSampleSize = a->RSampleSize; /* sample Size for R */
  size_t * RNhSize = a->RNhSize;
  size_t h;
  size_t i;
  size_t d;
  size_t Hi; 
  size_t Hj;
  size_t optHj;
  double maxRLocal;

  /* if there is nothing to do, do nothing */
  if( iter == 0) return; 

  double * sampleVar = malloc( sizeof(double) * H * dN);
  double * RLocal =    malloc( sizeof(double) * dN);
  double * RGlobal =   malloc( sizeof(double) * dN);

  double minSampleSize = 2;
    
  /* the goal here is fairly simple 
     * 0.0  we want to get the initial objective function, this is the typical max difference for CV's that violate the CV constraint
     * 1.0  now we randomly select a strata
     * 2.0  now we want to see which strata we can a 'draw' to.
   */
    
  /********* 0.0 pre calculate *************/

  /* 0.1 pre-calculate the sample variance */ 
  for( d = 0; d < dN; d++) 
    for( h = 0; h < H; h++) 
      sampleVar[H * d + h] = RNhSize[h] * RNhSize[h] * RV[d][h];   

  /* 0.2 copy R to RGlobal */
  for( d = 0; d < dN; d++) RGlobal[d] = R[d];

  /* iterate a fixed number of times */
  for( i = 0 ; i < iter; i++ ) {
      
    /* 1.0 randomly select a strata */
    Hi = SA_GETINDEX(H);
    optHj = H;
 

    /* only proceed if stratum Hi can be made smaller */ 
    if( RSampleSize[Hi] < minSampleSize + 1 ) {
      continue; 
    }

    /* 2.0 calculate objective function change for moving to each strata */ 
    for( Hj = 0; Hj < H; Hj++ ) {
      maxRLocal = 0;
 
      /* if the exchange would make the sample size too big, we don't do it */ 
      if( RSampleSize[Hj] + 1 > NhSize[Hj] ) continue;

      /* do not run over the selected state */
      if( Hj != Hi) {
        for( d = 0; d < dN; d++) { 
          RLocal[d] = 0;
      
          for( h = 0; h < H; h++) {
            if( h == Hj ) { 
              RLocal[d] += sampleVar[H*d+h] / (RSampleSize[h] + 1); 
            } else if( h == Hi ) {
              RLocal[d] += sampleVar[H*d+h] / (RSampleSize[h] - 1);
            } else { 
              RLocal[d] += sampleVar[H*d+h] / RSampleSize[h];
            }
          }

          RLocal[d] = sqrt(RLocal[d])/Total[d] - T[d];
  
          /* check if we are doing better than RGlobal[d] for every d that does not satisfy the CV constraint */ 
          if( RLocal[d] > 0 ) { 
           
            if( RLocal[d] > RGlobal[d] ) {
              maxRLocal = INFINITY; 
              break;
            } 
          }
        }
        /* finished iterating over all commodities */


        /* check if we are doing better */ 
        if( maxRLocal < INFINITY ) {
          optHj = Hj;
          for( d = 0; d < dN; d++) RGlobal[d] = RLocal[d];
        } 

      }
      
    }
    
    /* 3.0 if there are reductions in objective function make a move to minimize CV */
  
    if( optHj < H ) {  /* check if the final result minimized the objective function for all H */
      for( d = 0; d < dN; d++) R[d] = RGlobal[d]; /* assign RGlobal[d] to R[d] */
      RSampleSize[optHj] += 1.0;
      RSampleSize[Hi] -= 1.0;
    }

  }

  free(RLocal);
  free(RGlobal);
  free(sampleVar);
  return;
}



/* This is a function to check the change in substrata cost                   */
double minCV_costChange (
            size_t i,              /* new Index                               */
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * J,            /* new state                               */
            double * R,            /* new cost                                */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
    ) { 

  /* cast A back to something useable */
  minCV_adminStructPtr a = (minCV_adminStructPtr) A; 
  
  double *** C = a->C;
  size_t H = a->H;
  double * T = a->T;
  size_t * size = a->size;
  size_t * NhSize = a->NhSize;
  size_t * RNhSize = a->RNhSize;
  double ** V = a->V;
  double ** RV = a->RV;
  double * Total = a->Total;
  size_t Hj = a->Hj;               /* get our chosen strata */
  double * sampleSize = a->sampleSize;
  double * RSampleSize = a->RSampleSize;
  double delta = 0;
  size_t h;
  double fixedVar;
  size_t  Hi, d; 
  size_t preserveSatisfied = a->preserveSatisfied;

  /* figure out change in cost */
  Hi = I[i]; 
    
  /* before we continue on we need to ensure that RNhSize and RSample Size reflect the current State */ 
  for( h = 0; h < H; h++) {
    RNhSize[h] = NhSize[h]; 
    RSampleSize[h] = sampleSize[h]; 
  }
    
  /* proposed new population size for Hj */
  RNhSize[Hj] = NhSize[Hj] + size[i]; 

  /* proposed new population size for Hi */
  RNhSize[Hi] = NhSize[Hi] - size[i]; 
       

  /* update the CV differences */ 
  for(d = 0; d < dN; d++) {
    /* get the variance total for the otehr strata */
    fixedVar = 0;

    for( h = 0; h < H; h++) {
     if( (h != Hi) & (h != Hj) ) {
       RV[d][h] = V[d][h]; /* copy over any prior changes */
       fixedVar += NhSize[h] * NhSize[h] * RV[d][h] / sampleSize[h];   
     } 
    }

    /* update RV */
    RV[d][Hj] =  (V[d][Hj] * NhSize[Hj] * (NhSize[Hj] - 1) + C[d][i][Hj] ) / (double) ( (RNhSize[Hj] - 1) * RNhSize[Hj] );
    RV[d][Hi] =  (V[d][Hi] * NhSize[Hi] * (NhSize[Hi] - 1) - C[d][i][Hi] ) / (double) ( (RNhSize[Hi] - 1) * RNhSize[Hi] );
 

    /* get the distance between */
    R[d] =
      sqrt( fixedVar + 
       (RNhSize[Hj] * RNhSize[Hj]) / sampleSize[Hj] * RV[d][Hj]  + 
       (RNhSize[Hi] * RNhSize[Hi]) / sampleSize[Hi] * RV[d][Hi] 
      ) / Total[d] - T[d];
  } 
  
  /* update the allocation */ 
  minCV_sampleSizeChange ( A, R, dN, N, a->sampleIter);

  /* get the max change in CV */
  for( d = 0; d < dN; d++) {
    if( R[d] > 0 ) {

      /* if preserveSatisfied == 1 then any change that would violate a met 
       * constraint will be avoided
       */
      if( ( Q[d] <= 0 ) & (preserveSatisfied == 1) ) return( INFINITY );

      if( R[d] - Q[d] > delta ) {
        delta =R[d] - Q[d];
      }
    } 
  }

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
            size_t * J, /* new state, destructive */
            double * R, /* new cost, destructive */
            double * D, /* distance matrix */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N,   /* number of elements within a state */
            double * costChange /* cost Change */
            ) { 
  

  /* nothing to do, return */
  if( accept == 0) return;
 

  /* open the data structure and convert it to something useful */
  minCV_adminStructPtr a = (minCV_adminStructPtr) A; 
  size_t   d;
  double * sampleSize    = a->sampleSize;      /* sample size */
  size_t   H             = a->H;               /* number of strata */
  double * T             = a->T;               /* cv targets */

  /* record objective function and sample sizes for diagnostics */
  if ( accept == 2) { 
    /* record what strata the selected PSU is moving to */
    costChange[4] = a->Hj;

    /* add on sample size */
    for( d = 0; d < H; d++) costChange[5+d] = sampleSize[d]; 
    
    /* add on current obj function */
    for( d = 0; d < dN; d++) costChange[5 + H +d] = Q[d] + T[d]; 

    return;
  }
  
  /* optimize sample size */
  if ( accept == 3) { 
    //minCV_sampleSizeChange ( A, Q, R, dN, N, a->sampleIter);
    return;
  }
  

  size_t l, Hi; 
  double dil ; 

  double *** C = a->C;
  double **  V = a->V;
  double **  RV = a->RV;
  //size_t *   Nh = a->Nh;
  double *   NhAcres = a->NhAcres;
  size_t *   Nh = a->Nh;
  size_t *   NhSize = a->NhSize;
  size_t *   RNhSize = a->RNhSize;
  size_t *   size = a->size;
  double *   acres = a->acres;
  size_t     Hj = a->Hj;
  size_t **  L = a->L;
  size_t     k = a->k;
  double *   x = a->x;
  double * prob = a->prob;
  double * probMatrix = a->probMatrix;
  double totalProbability = a->totalProbability;
  double * RSampleSize = a->RSampleSize;

  /* figure out change in cost */
  Hi = I[i];
  
  /* update the strata assignment */
  I[i] = Hj;

  /* update prob */
  totalProbability = totalProbability - prob[i];
  prob[i] = 1 - probMatrix[ H * i + I[i] ];

  /* update total prob */
  a->totalProbability = totalProbability + prob[i];
  
  /* update L */
  /* need to rewrite this TODO, check terminal condition and dependencies on L*/
  l=0;
  while( L[Hi][l] != i ) l++;
  L[Hi][l] = N;

  l = 0;
  while( L[Hj][l] < N ) l++;
  L[Hj][l] = i;
 
  /* if it is an end point, push back the end point */ 
  if( L[Hj][l] == N+1 ) L[Hj][l+1] = N+1;  

  /* update sample size */
  for( d = 0; d < H; d++) sampleSize[d] = RSampleSize[d];

  /* update Q */ 
  for( d=0; d < dN; d++) Q[d] = R[d];

  if(Hi == Hj) return; 

  /* update R, C  and V */
  for( d=0; d < dN; d++) {
      
    /* update the variance, this must be done before C gets updated */
    V[d][Hj] = RV[d][Hj]; 
    V[d][Hi] = RV[d][Hi]; 

   for( l = 0; l < N; l++) {
     dil = getDistX(i,l,x,k,d,N,squaredEuclidianMeanDist); 
    
     /* Hj is i's new index */
     C[d][l][Hj] =  C[d][l][Hj] + dil; 
     C[d][l][Hi] =  C[d][l][Hi] - dil;  
    }
  }

  Nh[Hi]--;
  Nh[Hj]++;

  NhSize[Hi] = RNhSize[Hi];
  NhSize[Hj] = RNhSize[Hj];

  NhAcres[Hi] -= acres[i];
  NhAcres[Hj] += acres[i];
  
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


/* packSubstrata is a function that build all required data sets
 * and converts the output to the adminSubstrataStructPtr to
 * a void pointer to pass to sa. 
 */
void * minCV_packSubstrata( 
    size_t * I,     /* initial assignment                            */
    double * D,     /* distance matrix                               */
    double * x,     /* the data set */
    int * aInt,     /* integer admin data from file                  */
    double * aDbl,  /* double admin data from file                   */
    size_t dN,      /* number of distance matricies                  */
    size_t N,       /* number of elements within a state             */
    size_t k,       /* the size of an element of x                   */
    size_t NInt,    /* the number of items in the admin integer data */
    size_t NDbl     /* the number of items in the admin double data */
) {

  size_t i;
   
  /* allocate struct */
  minCV_adminStructPtr packedStruct = malloc( sizeof( minCV_adminStruct ) ); 
  
  /* H is the number of labels */
  size_t H =  minCV_labelCount( I, N );

  /* created to take care of the convert issue */ 
  size_t * size       = malloc(sizeof(size_t) * N );
  double * acres      = malloc(sizeof(double) * N );
  double * prob       = malloc(sizeof(double) * N );
  double * probMatrix = malloc(sizeof(double) * N * H );

  for( i=0; i < N; i++) {
    size[i]  =  (size_t) aInt[i];
    acres[i] =           aDbl[i];
  } 

  /* get the number of iterations to optimize sample size */
  packedStruct->sampleIter = (size_t) aInt[N];
  
  /* get the handling of changes that violate a met constraint */ 
  packedStruct->preserveSatisfied = (size_t) aInt[N + 1];

  /* Nh is the size of each label (vector) */
  size_t * Nh = minCV_labelTotalPSUs( I, N, H);
  
  /* NhSize is the total number of segments (vector) */
  size_t * NhSize = minCV_labelTotalSegments( I, N, H, size);
  size_t * RNhSize = malloc( sizeof(size_t) * H);

  /* NhSize is the total number of segments (vector) */
  double * NhAcres = minCV_labelTotalAcres( I, N, H, acres);
  
  /* NhMax is the largest label from Nh*/
  //NhMax = 2 * minCV_arrayMaxSize_t( Nh, H );
  size_t NhMax = N;

  /* creates a matrix of indexes with rows substrata, and colunns indexes */ 
  size_t ** L = minCV_labelCreateMaster( I, N, H, NhMax); 

  /* creates a matrix of differences between item i, and all the points in a 
   * stratum h.  Rows are items, and columns are strata.
   */

  /* get target variance values */
  double * T = malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    T[i] = aDbl[N + i]; 
 
  /* within variation */ 
  double * W = malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    W[i] = aDbl[N + dN + i]; 
   
  /* create variance matrix for each [commodity][strata] */
  double *** mu = minCV_createMeanMatrix( I, dN, N, k, H, x, L, Nh, NhSize);
  double *** V = minCV_createVarMatrix( I, dN, N, k, H, x, L, mu, Nh, NhSize);
 
  /* create a place to store temporary variances */ 
  double *** RV = double3DMatrix(dN,H,k); 

  for( i = 0; i < dN; i++) 
    for( h = 0; h < H; h++) 
      for( l = 0; l < k; l++) 
        RV[i][h][l] = V[i][h][l]; 


  /* create variance matrix for each [commodity]*/
  double * Total = malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    Total[i] = aDbl[N + 2*dN + i]; 
 
  /* get the sample Size */ 
  double * sampleSize = (double *) malloc( sizeof( double ) * H );
  double * RSampleSize = (double *) malloc( sizeof( double ) * H );
  for( i =0; i < H; i++) 
    sampleSize[i] = aDbl[N + 3*dN + i]; 
 
  /* get max prob */ 
  for( i =0; i < N; i++) 
    prob[i] = aDbl[N + 3*dN + H + i]; 
  packedStruct->prob = prob;

  /* get prob matrix */
  for( i =0; i < N*H; i++) 
    probMatrix[i] = aDbl[2*N + 3*dN + H + i]; 
  packedStruct->probMatrix = probMatrix;
  
  /* make assignments to struct */
  packedStruct->H = H;
  packedStruct->Nh = Nh;
  packedStruct->NhSize = NhSize;
  packedStruct->RNhSize = RNhSize;
  packedStruct->NhAcres = NhAcres;
  packedStruct->NhMax = NhMax;
  packedStruct->acres = acres;
  packedStruct->size = size;
  packedStruct->L = L;
  packedStruct->C = C;
  packedStruct->T = T;
  packedStruct->W = W;
  packedStruct->V = V;
  packedStruct->RV = RV;
  packedStruct->Total = Total;
  packedStruct->x = x;
  packedStruct->k = k;
  packedStruct->totalProbability = aDbl[NDbl - 3]; 
  packedStruct->temp = aDbl[NDbl-2];
  packedStruct->acreDif = aDbl[NDbl-1 ];
  packedStruct->sampleSize = sampleSize;
  packedStruct->RSampleSize = RSampleSize;
  
  return( (void *) packedStruct );
}




/* clean up for the substrata administrative data */
void minCV_deleteSubstrata( minCV_adminStructPtr  a, size_t dN, size_t N )
{
  int i;
  

  /* first the vectors */
  free( a->Nh );
  a->Nh = NULL;
  
  free( a->prob );
  a->prob = NULL;
  
  free( a->probMatrix );
  a->probMatrix = NULL;
  
  free( a->NhSize );
  a->NhSize = NULL;

  free( a->NhAcres );
  a->NhAcres = NULL;
  
  free( a->size );
  a->size = NULL;
  
  free( a->acres );
  a->acres = NULL;
  
  free( a->T );
  a->T = NULL;
 
 free( a->Total); 
  a->Total = NULL; 
 
  free( a->sampleSize); 
  a->sampleSize = NULL; 
  
  free( a->RSampleSize); 
  a->RSampleSize = NULL; 

  /* now the MDA */
  freeMDA( (void*) a->L, a->H);
  a->L = NULL; 
  
  freeMDA( (void*) a->V, dN);
  a->V = NULL; 
  
  freeMDA( (void*) a->RV, dN);
  a->RV = NULL; 
  
  for(i = 0; i < dN; i++)
    freeMDA( (void*) (a->C)[i], N);
  free( a->C );
  a->C = NULL;

  /* now the struct */
  free( a );
  a =  NULL;
}


/* function to check status */
void minCV_diag( 
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
    ) {

  double *** V;
  double *** mu;
  size_t H, NhMax, d;
  size_t * Nh;
  double * T;
  double * Total;
  double * acres;
  size_t * size;
  size_t * NhSize;
  double * NhAcres;
  minCV_adminStructPtr a; 

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  H = a->H;
  Nh = a->Nh;
  NhMax = a->NhMax;
  T = a->T;
  V = a->V;
  mu = a->mu;
  Total = a->Total;


  acres = a->acres;
  size = a->size;
  NhSize = a->NhSize;
  NhAcres = a->NhAcres;
  double * sampleSize = a->sampleSize; 


  Rprintf("\n************************* i = %d **************************\n",(int) i);
  Rprintf("\nNhMax: %d\n", (int) NhMax);
 
  for( d =0; d < dN; d++) 
    Rprintf("\nV[%d]\n",(int) d),
    printMatrixFullDbl(V[d], N, H ); 
  
  for( d =0; d < dN; d++) 
    Rprintf("\nmu[%d]\n",(int) d),
    printMatrixFullDbl(mu[d], N, H ); 


  Rprintf("\nQ\n");
  Rprintf("sqrt(n) * CV_j - TCV_j\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, Q[d]); 
  
  
  /* probMatrix */ 
  /*
  Rprintf("\nprobMatrix\n");
  for( d =0; d < N; d++) {
    Rprintf("%d:  ",(int) d);
  
    for( f =0; f < H; f++) 
      Rprintf("%f, ", probMatrix[H*d + f] );
    
    Rprintf("\n");
  } 
  */ 
  /* prob */ 
  /*
  Rprintf("\nprob\n");
  for( d =0; d < N; d++) 
    Rprintf("%d:  %f\n",(int) d,  prob[d]); 
   */ 
  
  Rprintf("\nNh\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %d\n",(int) d, (int) Nh[d]); 

  Rprintf("\nNhSize\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %d\n",(int) d, (int) NhSize[d]); 
  
  Rprintf("\nNhAcres\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %f\n",(int) d, NhAcres[d]); 
  
  Rprintf("\nSample Size\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %f\n",(int) d, sampleSize[d]); 
  
  Rprintf("\nT\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, T[d]); 
  
  Rprintf("\nTotal\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, Total[d]); 

}
  

/* SA helper functions */

/* count the number of lables */
size_t minCV_labelCount( size_t * label, size_t N ) {

  int i,H;

  /* enumerated labels are required with the maximum being the number of labels */
  H = label[0];
  for(i = 1; i < N; i++)
    if ( label[i] > H ) H = label[i];

  return(H + 1);
}



/* this is a function that creates and retrurns an array of sizes for the labels */
size_t * minCV_labelTotalPSUs ( size_t * label, size_t N, size_t H ) {

  size_t i;
  size_t * labelCount;

  labelCount = malloc( sizeof(size_t) * H);

  /* initialize */
  for(i = 0; i < H; i++) {
    labelCount[i] = 0;
  }

  /* aggregate */
  for(i = 0; i < N; i++) {
    labelCount[label[i]] ++;
  }

  return(labelCount);

}


/* this is a function that creates and retrurns an array of sizes for the labels */
size_t * minCV_labelTotalSegments ( size_t * label, size_t N, size_t H, size_t * size ) {

  size_t i;
  size_t * labelCount;

  labelCount = malloc( sizeof(size_t) * H);

  /* initialize */
  for(i = 0; i < H; i++) {
    labelCount[i] = 0;
  }

  /* aggregate */
  for(i = 0; i < N; i++) {
    labelCount[label[i]] += size[i];
  }

  return(labelCount);

}


/* this is a function that creates and retrurns an array of sizes for the labels */
double * minCV_labelTotalAcres ( size_t * label, size_t N, size_t H, double * acres ) {

  size_t i;
  double * labelCount;

  labelCount = malloc( sizeof(double) * H);

  /* initialize */
  for(i = 0; i < H; i++) {
    labelCount[i] = 0;
  }

  /* aggregate */
  for(i = 0; i < N; i++) {
    labelCount[label[i]] += acres[i];
  }

  return(labelCount);

}


/* quick function to determine the maximum int in an array */
size_t minCV_arrayMaxSize_t( size_t * a, size_t n ) {
  size_t i, max;


  max = a[0];

  for( i = 1; i < n; i++)
    if( a[i] > max ) max = a[i];

  return( max);
}


/* quick function to determine the maximum size an acceptible substrata can be */
size_t minCV_arrayMaxSubstrata( size_t * a, size_t n, size_t *segments, size_t *acres, size_t * NhSize , double * NhAcres ) {
  size_t i, max;


  max = a[0];

  for( i = 1; i < n; i++)
    if( a[i] > max ) max = a[i];

  return( max);
}


/* quick function to determine the maximum valued index in an array */
size_t minCV_arrayMaxIndexDbl( double * a, size_t n ) {
  size_t i, max;


  max = 0;

  for( i = 1; i < n; i++)
    if( a[i] > a[max] ) max = i;

  return( max);
}


/* this matrix identifies each item with all other items that can have labels assigned */
size_t ** minCV_labelCreateMaster( size_t * label, size_t N, size_t H, size_t NhMax ) {

  size_t i,j,a;

  /* allocate memory for InTo matrix */
  size_t ** L =  malloc(sizeof(size_t * ) * H );

  for( i = 0; i < H; i ++) {
    L[i] = malloc(sizeof(size_t) * NhMax );
  }

  /* set default */
  for( i = 0; i < H; i++ ) {
    for( j = 0; j < NhMax; j++) {
      L[i][j] = N+1;
    }
  }

  /* assign each item to index */
  for(i = 0; i < N; i++) {
    j = 0;
    a = label[i];
    while(L[a][j] != N+1) j=j+1;

    L[a][j] = i;
  }

  return(L);
}



