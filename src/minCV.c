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
  double ** V;
  double * T;
  double * Total;
  double * W;
  size_t H; 
  size_t * NhSize;
  double * sampleSize;

  
  /* variable used to sum up the total over each stratum */ 

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  V = a->V;
  T = a->T;
  Total = a->Total;
  W = a->W;
  H = a->H;
  NhSize = a->NhSize;
  sampleSize = a->sampleSize;


  /* value to get max Q */
  Q[dN] = -INFINITY;

  /* Q has been allocated but there is nothing in it yet */
  for(i=0; i < dN; i++) {
    Q[i] = 0;
    for(j=0; j < H; j++) {
      Q[i] += V[i][j] * NhSize[j]* NhSize[j]/( sampleSize[j] );
    }
    Q[i] = sqrt(Q[i]) / Total[i] - T[i];
    
    if( Q[i] > Q[dN] ) Q[dN] = Q[i];
  }


}


  
/* function to get possible moves under acreage constraints */
size_t minCV_getPossibleMovesAcres(
  size_t index,         /* index to move from */ 
  size_t * I,           /* indexs of assignments */
  size_t * possibleMoves,
  size_t H,
  double * acres,      /*  */
  double acreDif,      /*  */
  double * NhAcres,    /* Nh acres */
  size_t ** L,         /* Stratum assignment H x N */
  size_t * Nh         /* strata PSU count */
) {
  
  size_t  Hi, j;
  size_t k = 0; /* number of moves */

  double minAcres = NhAcres[0]; /* minimum acres for any stratum, initialize with the first Acres */
  
  /* get minimum acres */
  for( j = 0, minAcres = INFINITY;  j < H; j++)
    if( NhAcres[j] < minAcres ) minAcres = NhAcres[j];

  /* generate a canidate */ 
  Hi = I[index]; /* get stratum of canidate */

  /* if the change is going to be below the smallest number of acres  */
  /* use the acres of the proposed change to check acreDif instead of minAcres */
  if( NhAcres[Hi] - acres[index] < minAcres )  {
    for( j = 0; j < H; j++) {
      if(Hi == j) continue;
      /* add to possible moves if the move does not violate acre difference */
      if((1 - (minAcres - acres[index])/( NhAcres[j] + acres[index] )) < acreDif ) {
        possibleMoves[k]=j;
        k++;
      }
    }
  } else {
    for( j = 0; j < H; j++) {
      if(Hi == j) continue;
      /* add to possible moves if the move does not violate acre difference */
      if((1 - minAcres/( NhAcres[j] + acres[index] )) < acreDif ) {
        possibleMoves[k]=j;
        k++;
      }
    }
  }
      
  return(k);
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

  //size_t * possibleMoves; /* array to hold possible moves */

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  double * prob = a->prob;
  size_t * neighbors = a->neighbors;
  size_t nNeighbors = a->nNeighbors;
  double totalProbability = a->totalProbability;
  size_t i; 
  size_t j = N;

  while ( j >= N ) {
    /* generate possible move */
    i = minCV_getIndex( prob, totalProbability );
  
    /* can I make a potential move neighbors */
    j = minCV_getMoveNeighbor( i, I, prob, neighbors, nNeighbors, N);

  }

  /* save j */
  a->j = j;

  /* return index */
  return(i);
}



/* function to select neighbor */
/* if a neighbor cannot found it returns an integer > nNeighbors */
size_t minCV_getMoveNeighbor( 
    size_t i, 
    size_t * I, 
    double * prob, 
    size_t * neighbors, 
    size_t nNeighbors,
    size_t N
    ) {
  size_t j;
  size_t k = 0;
  size_t Hi = I[i];
  size_t currentNeighbor;
  double totalProb = 0;
  double total;

  /* get the number of neighbors with different stratum membership */
  for( j=0; j < nNeighbors; j++) { 

    /* get the current neighbor of i */
    currentNeighbor = neighbors[i*nNeighbors + j];

    /* get the sum of the weight and number of neighbors not in the same strata as i */
    if( I[ currentNeighbor ] != Hi ) {
      k++;
      totalProb += prob[currentNeighbor];
    }

  }

  /* no neighbors in different strata, quit */
  if( k == 0 ) return( N );   

  /* randomly select the neighbor proportional to its sampling weight */
  totalProb = runif(0,totalProb); 
  total = 0;

  /* find the neighbor that corresponds to the selected sampling weight */ 
  for( j=0; j < nNeighbors; j++) {

    /* get the first neighbor */
    currentNeighbor = neighbors[i*nNeighbors + j];
    
    /* if the neighbor is in a differest strata then add total */ 
    if( I[ currentNeighbor ] != Hi ) {
      total += prob[currentNeighbor];
    }
    

    if( total >= totalProb) break;  
  }
  
  return( currentNeighbor ) ;
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

  minCV_adminStructPtr a;


  /* cast A back to something useable */
  a = (minCV_adminStructPtr) A; 
  
  double *** C = a->C;
  double H = a->H;
  double * T = a->T;
  size_t * size = a->size;
  size_t * NhSize = a->NhSize;
  double ** V = a->V;
  double * Total = a->Total;
  size_t k = a->k;                 /* size of an element within a commodity */
  double * x = a->x;               /* the data */
  size_t j = a->j;                 /* get our chosen strata */
  double * sampleSize = a->sampleSize;
  double delta = 0;
  double dij;
  size_t h;
  double  fixedVar;
  size_t  Hi, Hj, d; 

  double NhSizeHj, NhSizeHi, nhSizeHi, nhSizeHj;

  /* figure out change in cost */
  Hi = I[i]; 
  Hj = I[j];

#ifdef DEBUG
  printf("\n(i=%d, j = %d)\n", (int) i, (int) j);
  printf("\n(Hi=%d,Hj=%d)\n", (int) Hi, (int) Hj);
#endif
  
  /*this calculates psuedo nh's for each stratum*/
  for(d = 0; d < dN; d++) {

#ifdef DEBUG 
   printf(" ==== d = %d ==== \n", (int) d);
   printf("dij = %f\n", dij);
   printf("Q[%zu] =%f\n",d, Q[d]);  
   printf("V[%zu][%zu] =%f\n",d,Hi, V[d][Hi]);  
   printf("V[%zu][%zu] =%f\n",d,Hj, V[d][Hj]);  
   printf("T[%zu] =%f\n",d, T[d]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,i,Hi, C[d][i][Hi]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,i,Hj, C[d][i][Hj]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,j,Hi, C[d][j][Hi]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,j,Hj, C[d][j][Hj]);  
   printf("size[%zu] =%zu\n",i, size[d]);  
   printf("NhSize[%zu] =%zu\n",Hi, NhSize[Hi]);  
   printf("NhSize[%zu] =%zu\n",Hj, NhSize[Hj]);  
#endif

    
   dij = getDistX(i,j,x,k,d,N,squaredEuclidianMeanDist); 

    /* get the variance total for the otehr strata */
    fixedVar = 0;

    for( h = 0; h < H; h++) {
     if( (h != Hi)  & (h != Hj) ) {
       fixedVar += NhSize[h] * NhSize[h] * V[d][h] / sampleSize[h];   
     } 
    }

    /* proposed new population size for Hj */
    NhSizeHj = (double) NhSize[Hj] + size[i] - size[j]; 

    /* proposed new population size for Hi */
    NhSizeHi = (double) NhSize[Hi] + size[j] - size[i]; 

    /* proposed new sample size for Hj */
    nhSizeHj = (double) sampleSize[Hj] + size[i] - size[j]; 

    /* proposed new sample size for Hi */
    nhSizeHi = (double) sampleSize[Hi] + size[j] - size[i]; 

/*
    printf("Proposed Var[%d][%d] = %f\n", (int) d, (int) Hj,  
        (V[d][Hj] * NhSize[Hj] * (NhSize[Hj] - 1) + C[d][i][Hj] - C[d][j][Hj] - dij )  /
         ( (NhSizeHj - 1) * NhSizeHj ) 
        ); 
    
    printf("Proposed Var[%d][%d] = %f\n", (int) d, (int) Hi,  
        (V[d][Hi] * NhSize[Hi] * (NhSize[Hi] - 1) - C[d][i][Hi] + C[d][j][Hi] -dij ) /
        ( (NhSizeHi - 1) * NhSizeHj )
        );
*/

    /* get the distance between */
    R[d] =
      sqrt( fixedVar + 
        (V[d][Hj] * NhSize[Hj] * (NhSize[Hj] - 1) + C[d][i][Hj] - C[d][j][Hj] - dij ) / 
        ( (NhSizeHj - 1) * (nhSizeHj) ) * NhSizeHj 
        + 
        (V[d][Hi] * NhSize[Hi] * (NhSize[Hi] - 1) - C[d][i][Hi] + C[d][j][Hi] -dij ) / 
        ( (NhSizeHi - 1) * (nhSizeHi) ) * NhSizeHi 
      ) / Total[d] - T[d];

    
    //Rprintf("%d: R[d]+T[d]=%f, Q[d]+T[d]=%f, R[d]=%f, Q[d]=%f, R[d]-Q[d]=%12.8f \n", (int) d, R[d] + T[d], Q[d] + T[d], -1 * R[d] + Q[d] );

    /*
    if( R[d] - Q[d] == 0 ) printf("R[d] - Q[d] == 0\n");
    if( R[d] - Q[d] < 0 )  printf("R[d] - Q[d] < 0\n");
    if( R[d] - Q[d] > 0 )  printf("R[d] - Q[d] > 0\n");
    */

    /* get the max change in CV */
    if( R[d] > 0 ) { 
      if( R[d] - Q[d] > 0 ) {
        delta++;
      }
    }
   
//    printf("R-Q[%zu] =%f\n",d, R[d] - Q[d] );  

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
            size_t accept, /* 1 if accepted 0 if not */
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            size_t * J, /* new state, destructive */
            double * R, /* new cost, destructive */
            double * D, /* distance matrix */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
            ) { 
  
  size_t l, Hi, Hj, d; 
  minCV_adminStructPtr a;
  double dil, djl, dij; 

  if( accept == 0) return;

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  double *** C = a->C;
  double **  V = a->V;
  //size_t *   Nh = a->Nh;
  double *   NhAcres = a->NhAcres;
  size_t *   NhSize = a->NhSize;
  size_t *   size = a->size;
  double *   acres = a->acres;
  double *   sampleSize = a->sampleSize;
  size_t     j = a->j;
  size_t **  L = a->L;
  size_t     k = a->k;
  size_t     H = a->H;
  double *   x = a->x;
  double * prob = a->prob;
  double * probMatrix = a->probMatrix;
  double totalProbability = a->totalProbability;

  /* figure out change in cost */
  Hi = I[i];
  Hj = I[j]; 
  
  /* update the strata assignment */
  I[i] = Hj;
  I[j] = Hi;

  /* update prob */

  totalProbability = totalProbability - prob[i] - prob[j];

  prob[i] = 1 - probMatrix[ H * i + I[i] ];
  prob[j] = 1 - probMatrix[ H * j + I[j] ];

  /* update total prob */
  a->totalProbability = totalProbability + prob[i] + prob[j];
  
  /* update L */
  l=0;
  while( L[Hi][l] != i ) l++;
//  Rprintf("replace i = %d, L[Hi][l] = %d\n", (int) i, (int) L[Hi][l]); 
  L[Hi][l] = j;

  l = 0;
  while( L[Hj][l] != j ) l++;
//  Rprintf("replace j = %d, L[Hi][l] = %d\n", (int) j, (int) L[Hj][l]); 
  L[Hj][l] = i;


#ifdef DEBUG
  printf(" Hi = %d, Hj = %d, i = %d, j = %d ",(int) Hi,(int) Hj, (int) i, (int) j); 
  printf("\nUPDATING!!!\n");
#endif
  


  if(Hi == Hj) return; 

  /* update the contribution */
  /* to update the contribution we need to subtract dik from the column that i's strata was
   * previously associated with for every k, we also need to add djk to that stratum.  This
   * needs to be done for each commodity */


  /* since it is an update we actually take the time to calculate things correctly */

  /* renember to update NhSize and Size */

  Q[dN]=R[dN]; 

  for( d=0; d < dN; d++) {
      
   Q[d] = R[d];


   //dij = getDist(i,j,D,d,N); 
   dij = getDistX(i,j,x,k,d,N,squaredEuclidianMeanDist); 
   /* update the variance, this must be done before C gets updated */


   V[d][Hj] = 
        (V[d][Hj] * NhSize[Hj] * (NhSize[Hj] - 1) + C[d][i][Hj] - C[d][j][Hj] - dij) / ( (double)  (NhSize[Hj] + size[i] -size[j] - 1) *(NhSize[Hj] + size[i] -size[j] ) );

    V[d][Hi] = 
        (V[d][Hi] * NhSize[Hi] * (NhSize[Hi] - 1) - C[d][i][Hi] + C[d][j][Hi] - dij) / ( (double)  (NhSize[Hi] + size[j] -size[i] - 1) * (NhSize[Hi] + size[j] -size[i] ) ); 



   for( l = 0; l < N; l++) {
     //dil = getDist(i,l,D,d,N); 
     //djl = getDist(j,l,D,d,N);
     dil = getDistX(i,l,x,k,d,N,squaredEuclidianMeanDist); 
     djl = getDistX(j,l,x,k,d,N,squaredEuclidianMeanDist); 
    

     /* Hj is j's old index */
     C[d][l][Hj] =  C[d][l][Hj] + dil; 
     C[d][l][Hi] =  C[d][l][Hi] - dil;  
     
     C[d][l][Hi] =  C[d][l][Hi] + djl; 
     C[d][l][Hj] =  C[d][l][Hj] - djl;  
    }
  }

  NhSize[Hi] += (size[j] - size[i]) ;
  NhSize[Hj] += (size[i] - size[j]) ;

  NhAcres[Hi] += (acres[j] - acres[i]);
  NhAcres[Hj] += (acres[i] - acres[j]);
  
  sampleSize[Hi] += (size[j] - size[i]) ;
  sampleSize[Hj] += (size[i] - size[j]) ;
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

  /*printf("(i,t,d): %d\t %f\t %f\n",(int) iter, temp, diff);*/

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
   
  size_t H;
  size_t * Nh;
  size_t * NhSize;
  double * NhAcres;
  size_t NhMax;
  size_t ** L;    /* label matrix */
  double *** C;   /* cost tensor */
  double * T;     /* Target Variance */
  double * W;     /* Within Variance */
  double ** V;    /* Matrix of Variance, [Commodity][strata] */ 
  double * Total; /* Matrix of Totals, [Commodity]*/ 
  size_t nNeighbors; /* number of neighbors */

  double * sampleSize;

  /* allocate struct */
  minCV_adminStructPtr packedStruct =  
    (minCV_adminStructPtr) malloc( sizeof( minCV_adminStruct ) ); 
  
  /* H is the number of labels */
  H =  minCV_labelCount( I, N );
  

  /* created to take care of the convert issue */ 
  size_t * size          = (size_t * ) malloc(sizeof(size_t) * N );
  double * acres         = (double * ) malloc(sizeof(size_t) * N );
  size_t * neighbors     = (size_t * ) malloc(sizeof(size_t) * N * aInt[N] );
  double * prob          = (double * ) malloc(sizeof(double) * N );
  double * probMatrix    = (double * ) malloc(sizeof(double) * N * H );

  for( i=0; i < N; i++) {
    size[i]  =  (size_t) aInt[i];
    acres[i] =           aDbl[i];
  } 
 
  /* nNeighbor is the number of neighbors */
  nNeighbors = (size_t) aInt[N];
  packedStruct->nNeighbors = nNeighbors;

  /* copy over neighbors */ 
  for( i=0; i < nNeighbors * N; i++) 
    neighbors[i] = (size_t) aInt[N + 1 + i];
  packedStruct->neighbors = neighbors;

  /* Nh is the size of each label (vector) */
  Nh = minCV_labelTotalPSUs( I, N, H);
  
  /* NhSize is the total number of segments (vector) */
  NhSize = minCV_labelTotalSegments( I, N, H, size);

  /* NhSize is the total number of segments (vector) */
  NhAcres = minCV_labelTotalAcres( I, N, H, acres);
  
  /* NhMax is the largest label from Nh*/
  NhMax = 2 * minCV_arrayMaxSize_t( Nh, H );

  /* creates a matrix of indexes with rows substrata, and colunns indexes */ 
  L = minCV_labelCreateMaster( I, N, H, NhMax); 

  /* creates a matrix of differences between item i, and all the points in a 
   * stratum h.  Rows are items, and columns are strata.
   */

  C = minCV_createContribMatrix( I, dN, N,k, H, x, L, Nh);
  
  /* get target variance values */
  T = (double *) malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    T[i] = aDbl[N + i]; 
 
  /* within variation */ 
  W = (double *) malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    W[i] = aDbl[N + dN + i]; 
   
  /* create variance matrix for each [commodity][strata] */
  V = minCV_createVarMatrix( I, dN, N, k, H, x, L, Nh, NhSize);
  
  /* create variance matrix for each [commodity]*/
  Total = (double *) malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    Total[i] = aDbl[N + 2*dN + i]; 
  
  sampleSize = (double *) malloc( sizeof( double ) * H );
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
  packedStruct->NhAcres = NhAcres;
  packedStruct->NhMax = NhMax;
  packedStruct->acres = acres;
  packedStruct->size = size;
  packedStruct->L = L;
  packedStruct->C = C;
  packedStruct->T = T;
  packedStruct->W = W;
  packedStruct->V = V;
  packedStruct->Total = Total;
  packedStruct->x = x;
  packedStruct->k = k;
  packedStruct->totalProbability = aDbl[NDbl - 3]; 
  packedStruct->temp = aDbl[NDbl-2];
  packedStruct->acreDif = aDbl[NDbl-1 ];
  packedStruct->sampleSize = sampleSize;
  
  return( (void *) packedStruct );
}




/* clean up for the substrata administrative data */
void minCV_deleteSubstrata( minCV_adminStructPtr  a, size_t dN, size_t N )
{
  int i;
  

  /* first the vectors */
  free( a->Nh );
  a->Nh = NULL;
  
  free( a->neighbors );
  a->neighbors = NULL;
  
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
  

  /* now the MDA */
  freeMDA( (void*) a->L, a->H);
  a->L = NULL; 
  
  freeMDA( (void*) a->V, dN);
  a->V = NULL; 
  
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

  double *** C;
  size_t H, NhMax, d,f;
  size_t * Nh;
  double * T;
  double ** V;
  double * Total;
  double * acres;
  size_t * size;
  size_t * NhSize;
  double * NhAcres;
  minCV_adminStructPtr a; 

  /* cast A back to somethine useable */
  a = (minCV_adminStructPtr) A; 
  C = a->C;
  H = a->H;
  Nh = a->Nh;
  NhMax = a->NhMax;
  T = a->T;
  V = a->V;
  Total = a->Total;

  double * prob = a->prob;
  double * probMatrix = a->probMatrix;

  size_t nNeighbors = a->nNeighbors;
  size_t * neighbors = a->neighbors;

  acres = a->acres;
  size = a->size;
  NhSize = a->NhSize;
  NhAcres = a->NhAcres;
  double * sampleSize = a->sampleSize; 


  Rprintf("\n************************* i = %d **************************\n",(int) i);
  Rprintf("\nNhMax: %d\n", (int) NhMax);
 
 /* 
  for( d =0; d < dN; d++) 
    Rprintf("\nC[%d]\n",(int) d),
    printMatrixFullDbl(C[d], N, H ); 
*/

  Rprintf("\nQ\n");
  Rprintf("sqrt(n) * CV_j - TCV_j\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, Q[d]); 
  Rprintf("max{ sqrt(n) * CV_j - TCV_j }\n");
  Rprintf("%d:  %f\n",(int) dN, Q[dN]); 
  
  /* neighbors */ 
  /*
  Rprintf("\nneighbors\n");
  for( d =0; d < N; d++) {
    Rprintf("%d:  ",(int) d);
  
    for( f =0; f < nNeighbors; f++) 
      Rprintf("%d, ",(int) neighbors[nNeighbors*d + f] );

    Rprintf("\n");
  } 
  */
  
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
  
  Rprintf("\nV\n");
    printMatrixFullDbl(V, dN, H ); 

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
  size_t ** L;

  /* allocate memory for InTo matrix */
  L = (size_t **) malloc(sizeof(size_t * ) * H );

  for( i = 0; i < H; i ++) {
    L[i] = (size_t *) malloc(sizeof(size_t) * NhMax );
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


/* creates a matrix of differences between item i, and all the points in a 
 * stratum h.  Rows are items, and columns are strata.
 */ 
double *** minCV_createContribMatrix( size_t * label, size_t dN, size_t N, size_t k, size_t H, double * x, size_t ** L, size_t * Nh)
{

  size_t i,j,m,l;
  double *** C;

  /* allocate memory for C matrix */
  C = (double ***) malloc(sizeof(double ** ) * dN );

  for( i = 0; i < dN; i ++)
  {
    C[i] = (double **) malloc(sizeof(double * ) * N );

    for( j = 0; j < N; j ++) {
      C[i][j] = (double *) malloc(sizeof(double) * H );
    }
  }


  /* aggregate */
  for( i = 0; i < dN; i++)
    for( j = 0; j < N; j++)
      for( m = 0; m < H; m++)
      {
      C[i][j][m] = 0;
      for( l = 0; l < Nh[m]; l++)
        C[i][j][m] += getDistX(j,L[m][l],x,k,i,N,squaredEuclidianMeanDist); 
        //C[i][j][m] += getDist(j,L[m][l],D,i,N); 
      }

  return(C);
}




/* function that creates a variance MDA (matrix) [commodity][strata] */
double ** minCV_createVarMatrix( size_t * label, size_t dN, size_t N, size_t k, size_t H, double * x, size_t ** L, size_t * Nh, size_t * NhSize)
{

  size_t i,j,m,l;
  double ** V;

  /* allocate memory for C matrix */
  V = (double **) malloc(sizeof(double * ) * dN );


  for( i = 0; i < dN; i++) V[i] = (double *) malloc(sizeof(double) * H );

  /* aggregate */
  for( i = 0; i < dN; i++)
    for( j = 0; j < H; j++)
    {
    /* set initial sum to 0 */
    V[i][j] = 0;
    for( m = 0; m < Nh[j]; m++)
      for( l = m+1; l < Nh[j]; l++) 
        V[i][j] += getDistX(L[j][m],L[j][l],x,k,i,N,squaredEuclidianMeanDist); 
        //V[i][j] += getDist(L[j][m],L[j][l],D,i,N); 
      V[i][j] = V[i][j] / ((NhSize[j] -1) * NhSize[j]);
    }
    
  return(V);
}




