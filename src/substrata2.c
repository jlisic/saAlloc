#include "substrata2.h"

/************************** SUBSTRATA PROBLEM FUNCTIONS ***********************/

/* This is a function to check the change in substrata cost                   */
void substrata2_init (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
            ) { 

  size_t i,j;
  substrata2_adminStructPtr a;
  double ** V;
  double * T;
  double * W;
  size_t H; 
    
  /* cast A back to somethine useable */
  a = (substrata2_adminStructPtr) A; 
  V = a->V;
  T = a->T;
  W = a->W;
  H = a->H;

  /* Q has been allocated but there is nothing in it yet */
  for(i=0; i < dN; i++) {
    Q[i] = 0;
    for(j=0; j < H; j++) Q[i] += V[i][j];
    Q[i] = ( Q[i] + W[i] ) / T[i];
    
  }
}


/* get a new random state */
size_t substrata2_randomState (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * J,            /* new state                               */
            double * R,            /* new cost                                */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
    ) { 

  size_t  Hi, H, j, k, i,l; 
  substrata2_adminStructPtr a;
  double *  acres; 
  double * NhAcres;
  double minAcres;
  double acreDif;

  size_t * possibleMoves; /* array to hold possible moves */

  /* cast A back to somethine useable */
  a = (substrata2_adminStructPtr) A; 
  H = a->H;
  acres = a->acres; 
  acreDif = a->acreDif; 
  NhAcres = a->NhAcres; 


  /* create some space to hold possible moves */
  possibleMoves = malloc(sizeof(size_t) * H); 

  /* get minimum acres */
  minAcres = NhAcres[0];
  for( j = 0, minAcres = INFINITY;  j < H; j++)
    if( NhAcres[j] < minAcres ) { 
      minAcres = NhAcres[j];
    }

  /* keep generating candiate moves until you try N*H times, then give up */ 
  l = 0;
  k = 0;
  while(k < 1) {

    /* generate a canidate */ 
    i = SA_GETINDEX(N);
    Hi = I[i];
    k = 0; /* number of places to move to */ 

    /* if the change is going to be below the smallest number of acres  */
    /* use the acres of the proposed change to check acreDif instead of minAcres */
    if( NhAcres[Hi] - acres[i] < minAcres )  {
      for( j = 0; j < H; j++) {
        if(Hi == j) continue;
        /* add to possible moves if the move does not violate acre difference */
        if((1 - (minAcres - acres[i])/( NhAcres[j] + acres[i] )) < acreDif ) {
          possibleMoves[k]=j;
          k++;
        }
      }
    } else {
      for( j = 0; j < H; j++) {
        if(Hi == j) continue;
        /* add to possible moves if the move does not violate acre difference */
        if((1 - minAcres/( NhAcres[j] + acres[i] )) < acreDif ) {
          possibleMoves[k]=j;
          k++;
        }
      }
    }
      

    if( l > N * H ) {
      Rprintf("No Movement Likely Possible, Terminating\n");
      return(N+1);
    }
    l++;
  }


  /* add to our data structure our next move */ 
  a->j = possibleMoves[SA_GETINDEX(k)];

  free(possibleMoves);

  return(i);
}


/* This is a function to check the change in substrata cost                   */
double substrata2_costChange (
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

  size_t  Hi, Hj, d; 
  substrata2_adminStructPtr a;
  double *** C;
  double * T;
  size_t * NhSize;
  size_t * size;
  double delta;
  double ** V;


  /* cast A back to somethine useable */
  a = (substrata2_adminStructPtr) A; 
  C = a->C;
  T = a->T;
  size = a->size;
  NhSize = a->NhSize;
  V = a->V;

  Hj = a->j; /* get our chosen strata */

  /* figure out change in cost */
  Hi = I[i]; 


  /*this calculates psuedo nh's for each stratum*/
  for(d = 0; d < dN; d++) {

    /* get the distance between */
    R[d] = 
      Q[d] 
      - 
      (V[d][Hi] + V[d][Hj])/T[d] 
      +
      (
        (V[d][Hj] + C[d][i][Hj]/(NhSize[Hj]*(NhSize[Hj] - 1))) * 
          NhSize[Hj] * (NhSize[Hj] - 1) / ( (NhSize[Hj] + size[i]) * (NhSize[Hj] + size[i] - 1))
        + 
        (V[d][Hi] - C[d][i][Hi]/(NhSize[Hi]*(NhSize[Hi] - 1))) *
          NhSize[Hi] * (NhSize[Hi] - 1) / ( (NhSize[Hi] - size[i]) * (NhSize[Hi] - size[i] - 1))
      )/T[d];

   }
   
  /* get the change in sample size */
  delta = 
    R[substrata2_arrayMaxIndexDbl( R, dN )]
    -
    Q[substrata2_arrayMaxIndexDbl( Q, dN )] ;
  
    return(delta);
}


/* update from change in i function */
void substrata2_update (
            size_t accept, /* 1 if accepted 0 if not, 2 for aux function */
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
  
  if( accept == 0) return;
  
  /* cast A back to somethine useable */
  substrata2_adminStructPtr a = (substrata2_adminStructPtr) A; 
 
  if( accept == 2) {
    costChange[4] = a->j; /* add on what it is moving to */
    return;
  }
  
  double dik;
   
  /* map data from data structure */ 
  size_t k, d; 
  double *** C = a->C;
  double ** V = a->V;
  size_t * Nh = a->Nh;
  size_t * NhSize = a->NhSize;
  double * NhAcres = a->NhAcres;
  size_t * size = a->size;
  double * acres = a->acres;
  size_t j  = a->j;

  /* figure out change in cost */
  size_t Hi = I[i];
  size_t Hj = j; 
  
  /* update the strata assignment */
  I[i] = Hj;

  if(Hi == Hj) return; 

  /* update the contribution */
  /* to update the contribution we need to subtract dik from the column that i's strata was
   * previously associated with for every k, we also need to add djk to that stratum.  This
   * needs to be done for each commodity */


  /* since it is an update we actually take the time to calculate things correctly */

  /* renember to update NhSize and Size */

  for( d=0; d < dN; d++) {
  
    /* update the variance */
   V[d][Hj] =  (V[d][Hj] + C[d][i][Hj]/(NhSize[Hj]*(NhSize[Hj] - 1))) * 
       NhSize[Hj] * (NhSize[Hj] - 1) / ( (NhSize[Hj] + size[i]) * (NhSize[Hj] + size[i] - 1));

   V[d][Hi] =  (V[d][Hi] - C[d][i][Hi]/(NhSize[Hi]*(NhSize[Hi] - 1))) *
       NhSize[Hi] * (NhSize[Hi] - 1) / ( (NhSize[Hi] - size[i]) * (NhSize[Hi] - size[i] - 1));
  
   Q[d] = R[d];

   for( k = 0; k < N; k++) {
    
     dik = getDist(i,k,D,d,N); 

     /* Hj is j's old index */
     C[d][k][Hj] =  C[d][k][Hj] + dik; 
     C[d][k][Hi] =  C[d][k][Hi] - dik;  

     }
  }

  Nh[Hi]--;
  Nh[Hj]++;

  NhSize[Hi] -= size[i];
  NhSize[Hj] += size[i];

  NhAcres[Hi] -= acres[i];
  NhAcres[Hj] += acres[i];

  return;
}


/* cooling schedule                                                           */
double substrata2_cool ( 
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
  substrata2_adminStructPtr a = (substrata2_adminStructPtr) A;

  temp = a->temp;

  return( exp( (1+iter) * diff /  (-1* temp) ) ) ; 
}  


/* packSubstrata is a function that build all required data sets
 * and converts the output to the adminSubstrataStructPtr to
 * a void pointer to pass to sa. 
 */
void * substrata2_packSubstrata( 
    size_t * I,     /* initial assignment                            */
    double * D,     /* distance matrix                               */
    int * aInt,     /* integer admin data from file                  */
    double * aDbl,  /* double admin data from file                   */
    size_t dN,      /* number of distance matricies                  */
    size_t N,       /* number of elements within a state             */
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


  /* allocate struct */
  substrata2_adminStructPtr packedStruct = 
    malloc( sizeof( substrata2_adminStruct ) ); 

  /* created to take care of the convert issue */ 
  size_t * size   =  calloc(N, sizeof(size_t) );
  double * acres  =  calloc(N, sizeof(double ) );

  for( i=0; i < N; i++) 
  {
    size[i]  =  (size_t) aInt[i];
    acres[i] =           aDbl[i];
  } 

  /* H is the number of labels */
  H =  substrata2_labelCount( I, N );

  /* Nh is the size of each label (vector) */
  Nh = substrata2_labelTotalPSUs( I, N, H);
  
  /* NhSize is the total number of segments (vector) */
  NhSize = substrata2_labelTotalSegments( I, N, H, size);

  /* NhSize is the total number of segments (vector) */
  NhAcres = substrata2_labelTotalAcres( I, N, H, acres);
  
  /* NhMax is the largest label from Nh*/
  NhMax = 2 * substrata2_arrayMaxSize_t( Nh, H );

  /* creates a matrix of indexes with rows substrata, and colunns indexes */ 
  L = substrata2_labelCreateMaster( I, N, H, NhMax); 

  /* creates a matrix of differences between item i, and all the points in a 
   * stratum h.  Rows are items, and columns are strata.
   */ 
  C = substrata2_createContribMatrix( I, dN, N, H, D, L, Nh);
  
  /* get target variance values */
  T = malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    T[i] = aDbl[N + i]; 
 
  /* within variation */ 
  W = malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    W[i] = aDbl[N + dN + i]; 
   
  /* create variance matrix for each [commodity][strata] */
  V = substrata2_createVarMatrix( I, dN, N, H, D, L, Nh, NhSize);
  

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
  packedStruct->temp = aDbl[NDbl-2];
  packedStruct->acreDif = aDbl[NDbl-1 ];
  
  return( (void *) packedStruct );
}


/* clean up for the substrata administrative data */
void substrata2_deleteSubstrata( substrata2_adminStructPtr  a, size_t dN, size_t N )
{
  int i;

  /* first the vectors */
  free( a->Nh );
  a->Nh = NULL;

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
void substrata2_diag( 
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
    ) {

  size_t H, d;
  size_t * Nh;
  double * T;
  double ** V;
  size_t * NhSize;
  substrata2_adminStructPtr a; 

  /* cast A back to somethine useable */
  a = (substrata2_adminStructPtr) A; 
  H = a->H;
  Nh = a->Nh;
  T = a->T;
  V = a->V;

  NhSize = a->NhSize;
  
  Rprintf("\n************************* i = %d **************************\n",(int) i);

  
  Rprintf("\nQ\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, Q[d]); 
  
  Rprintf("\nNh\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %d\n",(int) d, (int) Nh[d]); 

  Rprintf("\nNhSize\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %d\n",(int) d, (int) NhSize[d]); 
  
  Rprintf("\nT\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, T[d]); 
  
  Rprintf("\nV\n");
    printMatrixFullDbl(V, dN, H ); 

}
  

/* SA helper functions */

/* count the number of lables */
size_t substrata2_labelCount( size_t * label, size_t N ) {

  int i,H;

  /* enumerated labels are required with the maximum being the number of labels */
  H = label[0];
  for(i = 1; i < N; i++)
    if ( label[i] > H ) H = label[i];

  return(H + 1);
}



/* this is a function that creates and retrurns an array of sizes for the labels */
size_t * substrata2_labelTotalPSUs ( size_t * label, size_t N, size_t H ) {

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
size_t * substrata2_labelTotalSegments ( size_t * label, size_t N, size_t H, size_t * size ) {

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
double * substrata2_labelTotalAcres ( size_t * label, size_t N, size_t H, double * acres ) {

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
size_t substrata2_arrayMaxSize_t( size_t * a, size_t n ) {
  size_t i, max;


  max = a[0];

  for( i = 1; i < n; i++)
    if( a[i] > max ) max = a[i];

  return( max);
}


/* quick function to determine the maximum size an acceptible substrata can be */
size_t substrata2_arrayMaxSubstrata( size_t * a, size_t n, size_t *segments, size_t *acres, size_t * NhSize , double * NhAcres ) {
  size_t i, max;


  max = a[0];

  for( i = 1; i < n; i++)
    if( a[i] > max ) max = a[i];

  return( max);
}


/* quick function to determine the maximum valued index in an array */
size_t substrata2_arrayMaxIndexDbl( double * a, size_t n ) {
  size_t i, max;


  max = 0;

  for( i = 1; i < n; i++)
    if( a[i] > a[max] ) max = i;

  return( max);
}


/* this matrix identifies each item with all other items that can have labels assigned */
size_t ** substrata2_labelCreateMaster( size_t * label, size_t N, size_t H, size_t NhMax ) {

  size_t i,j,a;
  size_t ** L;

  /* allocate memory for InTo matrix */
  L = malloc( sizeof(size_t * ) * H );

  for( i = 0; i < H; i ++) {
    L[i] = malloc( sizeof(size_t) * NhMax );
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
double *** substrata2_createContribMatrix( size_t * label, size_t dN, size_t N, size_t H, double * D, size_t ** L, size_t * Nh)
{

  size_t i,j,k,l;
  double *** C;

  /* allocate memory for C matrix */
  C = malloc(sizeof(double ** ) * dN );

  for( i = 0; i < dN; i ++)
  {
    C[i] = malloc(sizeof(double * ) * N );

    for( j = 0; j < N; j ++) {
      C[i][j] = malloc(sizeof(double) * H );
    }
  }

  /* aggregate */
  for( i = 0; i < dN; i++)
    for( j = 0; j < N; j++)
      for( k = 0; k < H; k++)
      {
      C[i][j][k] = 0;
      for( l = 0; l < Nh[k]; l++)
        C[i][j][k] += getDist(j,L[k][l],D,i,N); 
      }

  return(C);
}


/* function that creates a variance MDA (matrix) [commodity][strata] */
double ** substrata2_createVarMatrix( size_t * label, size_t dN, size_t N, size_t H, double * D, size_t ** L, size_t * Nh, size_t * NhSize)
{

  size_t i,j,k,l;
  double ** V;

  /* allocate memory for C matrix */
  V = malloc(sizeof(double * ) * dN );


  for( i = 0; i < dN; i++) V[i] = malloc(sizeof(double) * H );

  /* aggregate */
  for( i = 0; i < dN; i++)
    for( j = 0; j < H; j++) {
      /* set initial sum to 0 */
      V[i][j] = 0;
      for( k = 0; k < Nh[j]; k++)
        for( l = k+1; l < Nh[j]; l++) 
          V[i][j] += getDist(L[j][k],L[j][l],D,i,N); 

      V[i][j] = V[i][j] / ((NhSize[j] -1) * NhSize[j]);
    }
    
  return(V);
}


