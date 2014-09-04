#include "optCV.h"

/************************** SUBSTRATA PROBLEM FUNCTIONS ***********************/

/* This is a function to check initialize the structures and variables for
 * the sa proceedure                                                          */
void optCV_init (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
            ) { 

  size_t i,j;
  optCV_adminStructPtr a;
  double ** V;
  double * T;
  double ** Total;
  double * W;
  size_t H; 
  size_t * NhSize;

  
  /* variable used to sum up the total over each stratum */ 

  /* cast A back to somethine useable */
  a = (optCV_adminStructPtr) A; 
  V = a->V;
  T = a->T;
  Total = a->Total;
  W = a->W;
  H = a->H;
  NhSize = a->NhSize;

  /* Q has been allocated but there is nothing in it yet */
  for(i=0; i < dN; i++) {
    Q[i] = 0;
    Q[dN + i] = 0;
    for(j=0; j < H; j++) {
      Q[i] += V[i][j] * NhSize[j] ;
      Q[dN + i] += Total[i][j]; 
    }
    Q[dN+i] = Q[dN+i] * Q[dN+i];
  }
}




/* get a new random state */
size_t optCV_randomState (
            size_t * I,            /* current state                           */
            double * Q,            /* current cost                            */
            size_t * J,            /* new state                               */
            double * R,            /* new cost                                */
            double * D,            /* distance matrix                         */
            void * A,              /* administrative data                     */
            size_t dN,             /* number of distance matricies            */
            size_t N               /* number of elements within a state       */
    ) { 

  size_t  Hi, Hj, H, j, k, i,l; 
  optCV_adminStructPtr a;
  double *  acres; 
  double * NhAcres;
  double minAcres;
  double acreDif;
  size_t ** L;
  size_t * Nh;

  size_t * possibleMoves; /* array to hold possible moves */

  /* cast A back to somethine useable */
  a = (optCV_adminStructPtr) A; 
  H = a->H;
  acres = a->acres; 
  acreDif = a->acreDif; 
  NhAcres = a->NhAcres; 
  L = a->L;
  Nh = a->Nh;


  /* create some space to hold possible moves */
  possibleMoves = (size_t *) malloc(sizeof(size_t) * H); 

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
#ifdef DEBUG
    printf("k=%zu : canidate i=%zu, N= %zu, test=%f\n", k,i,N, runif(0.0,1.0)  );
#endif
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
      #ifndef CLI
      Rprintf("No Movement Likely Possible, Terminating\n");
      #endif 
      #ifdef CLI
      printf("No Movement Likely Possible, Terminating\n");
      #endif 
      
      return(N+1);
    }
    l++;
  }


  /* add to our data structure our next move */ 
  Hj = possibleMoves[SA_GETINDEX(k)];
 
  a->j = L[Hj][ SA_GETINDEX( Nh[Hj] ) ]; 

  free(possibleMoves);

  return(i);
}


/* This is a function to check the change in substrata cost                   */
double optCV_costChange (
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
  optCV_adminStructPtr a;
  double *** C;
  double * T;
  size_t * NhSize;
  size_t * size;
  size_t k;  /* size of an element within a commodity */
  size_t j;
  double delta = -1;
  double ** V;
  double ** Total;
  double *x; /* the data */
  double dij;


  /* cast A back to something useable */
  a = (optCV_adminStructPtr) A; 
  C = a->C;
  T = a->T;
  size = a->size;
  NhSize = a->NhSize;
  V = a->V;
  Total = a->Total;
  k = a->k;
  x = a->x;
  j = a->j; /* get our chosen strata */

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
   printf("Q[%zu] =%f\n",d, Q[d]);  
   printf("V[%zu][%zu] =%f\n",d,Hi, V[d][Hi]);  
   printf("V[%zu][%zu] =%f\n",d,Hj, V[d][Hj]);  
   printf("Total[%zu][%zu] =%f\n",d,Hi, Total[d][Hi]);  
   printf("Total[%zu][%zu] =%f\n",d,Hj, Total[d][Hj]);  
   printf("T[%zu] =%f\n",d, T[d]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,i,Hi, C[d][i][Hi]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,i,Hj, C[d][i][Hj]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,j,Hi, C[d][j][Hi]);  
   printf("C[%zu][%zu][%zu] =%f\n",d,j,Hj, C[d][j][Hj]);  
   printf("size[%zu] =%zu\n",i, size[d]);  
   printf("NhSize[%zu] =%zu\n",Hi, NhSize[Hi]);  
   printf("NhSize[%zu] =%zu\n",Hj, NhSize[Hj]);  
#endif

    /* update the Mean */
    /* well not really */
      R[d+dN] = Q[d+dN]; 

    dij = getDist(i,j,D,d,N); 

    /* get the distance between */
    R[d] =
      (V[d][Hj] * NhSize[Hj] * (NhSize[Hj] - 1) + C[d][i][Hj] - C[d][j][Hj] -dij ) / (double) (NhSize[Hj] + size[i] -size[j] - 1)
      + 
      (V[d][Hi] * NhSize[Hi] * (NhSize[Hi] - 1) - C[d][i][Hi] + C[d][j][Hi] -dij ) / (double) (NhSize[Hi] + size[j] -size[i] - 1)
      ; 

    /* get the max change in CV */
#ifdef DEBUG
    printf( 
        "R[%zu]/R[%zu + %zu] - Q[%zu]/Q[%zu + %zu] = %f/%f - %f/%f = %f\n", 
        d, dN,d, d,dN,d,
        R[d],R[dN + d], Q[d],Q[dN + d],
        R[d]/R[dN + d] - Q[d]/Q[dN + d]
        );
#endif
    if( R[d]/R[dN + d] - Q[d]/Q[dN + d] > delta )  
        delta = R[d]/R[dN + d] - Q[d]/Q[dN + d];

   }
 
#ifdef DEBUG 
  printf("\ndelta: %f\n",delta); 
#endif
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

void optCV_update (
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
  
  size_t j,l, Hi, Hj, d; 
  optCV_adminStructPtr a;
  double *** C;
  double ** V;
  size_t * Nh;
  size_t * NhSize;
  size_t * size;
  double * NhAcres;
  double * acres;
  double dil, djl, dij; 

  if( accept == 0) return;

  /* cast A back to somethine useable */
  a = (optCV_adminStructPtr) A; 
  C = a->C;
  V = a->V;
  Nh = a->Nh;
  NhAcres = a->NhAcres;
  NhSize = a->NhSize;
  size = a->size;
  acres = a->acres;
  j  = a->j;

  /* figure out change in cost */
  Hi = I[i];
  Hj = I[j]; 
  
  /* update the strata assignment */
  I[i] = Hj;
  I[j] = Hi;

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
 

  for( d=0; d < dN; d++) {
      
   Q[d] = R[d];


   dij = getDist(i,j,D,d,N); 
   /* update the variance, this must be done before C gets updated */
   V[d][Hj] = 
        (V[d][Hj] * NhSize[Hj] * (NhSize[Hj] - 1) + C[d][i][Hj] - C[d][j][Hj] - dij) / ( (double)  (NhSize[Hj] + size[i] -size[j] - 1) *(NhSize[Hj] + size[i] -size[j] ) );

    V[d][Hi] = 
        (V[d][Hi] * NhSize[Hi] * (NhSize[Hi] - 1) - C[d][i][Hi] + C[d][j][Hi] - dij) / ( (double)  (NhSize[Hi] + size[j] -size[i] - 1) * (NhSize[Hi] + size[j] -size[i] ) ); 



   for( l = 0; l < N; l++) {
      dil = getDist(i,l,D,d,N); 
      djl = getDist(j,l,D,d,N);
    

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
}


/* cooling schedule                                                           */
double optCV_cool ( 
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
  optCV_adminStructPtr a;
  
  a = (optCV_adminStructPtr) A;
  temp = a->temp;

  /*printf("(i,t,d): %d\t %f\t %f\n",(int) iter, temp, diff);*/

  /*return( exp( (log(1+iter)) * diff /  (-1* temp) ) ); */ 
  return( exp( (1+iter) * diff /  (-1* temp) ) ) ; 
}  


/* packSubstrata is a function that build all required data sets
 * and converts the output to the adminSubstrataStructPtr to
 * a void pointer to pass to sa. 
 */
void * optCV_packSubstrata( 
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
  double ** Total;/* Matrix of Totals, [Commodity][strata] */ 

  size_t * size;  /* number of assumed segments in a psu */
  double * acres; /* acres for each psu */


  optCV_adminStructPtr packedStruct;

  /* created to take care of the convert issue */ 
  size          = (size_t * ) malloc(sizeof(size_t) * N );
  acres          = (double * ) malloc(sizeof(size_t) * N );

  for( i=0; i < N; i++) 
  {
    size[i]  =  (size_t) aInt[i];
    acres[i] =           aDbl[i];
  } 

  /* H is the number of labels */
  H =  optCV_labelCount( I, N );

  /* Nh is the size of each label (vector) */
  Nh = optCV_labelTotalPSUs( I, N, H);
  
  /* NhSize is the total number of segments (vector) */
  NhSize = optCV_labelTotalSegments( I, N, H, size);

#ifdef DEBUG
  printf("\n Obs: %zu", N);
  printf("\n Vars: %zu", dN);
  printf("\n H: %zu",H); 

 printf("\nNhSize: ");
  for( i = 0; i < H; i++)
    printf(" %d ", (int) NhSize[i]);
  printf("\n");
#endif

  /* NhSize is the total number of segments (vector) */
  NhAcres = optCV_labelTotalAcres( I, N, H, acres);
  
  /* NhMax is the largest label from Nh*/
  NhMax = 2 * optCV_arrayMaxSize_t( Nh, H );

  /* creates a matrix of indexes with rows substrata, and colunns indexes */ 
  L = optCV_labelCreateMaster( I, N, H, NhMax); 

  /* creates a matrix of differences between item i, and all the points in a 
   * stratum h.  Rows are items, and columns are strata.
   */ 
  C = optCV_createContribMatrix( I, dN, N, H, D, L, Nh);
  
  /* get target variance values */
  T = (double *) malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    T[i] = aDbl[N + i]; 
 
  /* within variation */ 
  W = (double *) malloc( sizeof( double ) * dN );
  for( i =0; i < dN; i++) 
    W[i] = aDbl[N + dN + i]; 
   
  /* create variance matrix for each [commodity][strata] */
  V = optCV_createVarMatrix( I, dN, N, H, D, L, Nh, NhSize);
  
  /* create variance matrix for each [commodity][strata] */
  Total = optCV_createTotalMatrix( I, k, dN, N, H, x);
  
  /* allocate struct */
  packedStruct = 
    (optCV_adminStructPtr) malloc( sizeof( optCV_adminStruct ) ); 

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
  packedStruct->temp = aDbl[NDbl-2];
  packedStruct->acreDif = aDbl[NDbl-1 ];
  
  printf("\nPacking\n Temp = %f Acre dif = %f", packedStruct->temp, packedStruct->acreDif);

  return( (void *) packedStruct );
}


/* clean up for the substrata administrative data */
void optCV_deleteSubstrata( optCV_adminStructPtr  a, size_t dN, size_t N )
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
  
  freeMDA( (void*) a->Total, dN);
  a->Total = NULL; 

  for(i = 0; i < dN; i++)
    freeMDA( (void*) (a->C)[i], N);
  free( a->C );
  a->C = NULL;

  /* now the struct */
  free( a );
  a =  NULL;
}


/* function to check status */
void optCV_diag( 
            size_t i,   /* new Index */
            size_t * I, /* current state */
            double * Q, /* current cost */
            void * A,   /* administrative data */
            size_t dN,  /* number of distance matricies */
            size_t N    /* number of elements within a state */
    ) {

  double *** C;
  size_t H, NhMax, d;
  size_t * Nh;
  double * T;
  double ** V;
  double ** Total;
  double * acres;
  size_t * size;
  size_t * NhSize;
  double * NhAcres;
  optCV_adminStructPtr a; 

  /* cast A back to somethine useable */
  a = (optCV_adminStructPtr) A; 
  C = a->C;
  H = a->H;
  Nh = a->Nh;
  NhMax = a->NhMax;
  T = a->T;
  V = a->V;
  Total = a->Total;

  acres = a->acres;
  size = a->size;
  NhSize = a->NhSize;
  NhAcres = a->NhAcres;


  #ifdef CLI
  printf("\n************************* i = %d **************************\n",(int) i);
  printf("\nNhMax: %d\n", (int) NhMax);

  /*
  for( d =0; d < dN; d++) 
    printf("\nC[%d]\n",(int) d),
    printMatrixFullDbl(C[d], N, H ); 
  
  printf("\nI,size,acres\n");
  for( d =0; d < N; d++) 
    printf("%d:  %d  %d  %f\n",(int) d, (int) I[d], (int) size[d], acres[d]); 
 */
  
  printf("\nQ\n");
  printf("Var\n");
  for( d =0; d < dN; d++) 
    printf("%d:  %f\n",(int) d, Q[d]); 
  printf("total\n");
  for( d =dN; d < 2*dN; d++) 
    printf("%d:  %f\n",(int) d, Q[d]); 
  
  printf("\nNh\n");
  for( d =0; d < H; d++) 
    printf("%d:  %d\n",(int) d, (int) Nh[d]); 

  printf("\nNhSize\n");
  for( d =0; d < H; d++) 
    printf("%d:  %d\n",(int) d, (int) NhSize[d]); 
  
  printf("\nNhAcres\n");
  for( d =0; d < H; d++) 
    printf("%d:  %f\n",(int) d, NhAcres[d]); 
  
  printf("\nT\n");
  for( d =0; d < dN; d++) 
    printf("%d:  %f\n",(int) d, T[d]); 
  
  printf("\nV\n");
    printMatrixFullDbl(V, dN, H ); 

  printf("\nTotal\n");
    printMatrixFullDbl(Total, dN, H ); 
  #endif
  #ifndef CLI
  Rprintf("\n************************* i = %d **************************\n",(int) i);
  Rprintf("\nNhMax: %d\n", (int) NhMax);


  Rprintf("\nQ\n");
  Rprintf("Var\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, Q[d]); 
  Rprintf("total\n");
  for( d =dN; d < 2*dN; d++) 
    Rprintf("%d:  %f\n",(int) d, Q[d]); 
  
  Rprintf("\nNh\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %d\n",(int) d, (int) Nh[d]); 

  Rprintf("\nNhSize\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %d\n",(int) d, (int) NhSize[d]); 
  
  Rprintf("\nNhAcres\n");
  for( d =0; d < H; d++) 
    Rprintf("%d:  %f\n",(int) d, NhAcres[d]); 
  
  Rprintf("\nT\n");
  for( d =0; d < dN; d++) 
    Rprintf("%d:  %f\n",(int) d, T[d]); 
  
  Rprintf("\nV\n");
    printMatrixFullDbl(V, dN, H ); 

  Rprintf("\nTotal\n");
    printMatrixFullDbl(Total, dN, H ); 

  #endif

}
  

/* SA helper functions */

/* count the number of lables */
size_t optCV_labelCount( size_t * label, size_t N ) {

  int i,H;

  /* enumerated labels are required with the maximum being the number of labels */
  H = label[0];
  for(i = 1; i < N; i++)
    if ( label[i] > H ) H = label[i];

  return(H + 1);
}



/* this is a function that creates and retrurns an array of sizes for the labels */
size_t * optCV_labelTotalPSUs ( size_t * label, size_t N, size_t H ) {

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
size_t * optCV_labelTotalSegments ( size_t * label, size_t N, size_t H, size_t * size ) {

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
double * optCV_labelTotalAcres ( size_t * label, size_t N, size_t H, double * acres ) {

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
size_t optCV_arrayMaxSize_t( size_t * a, size_t n ) {
  size_t i, max;


  max = a[0];

  for( i = 1; i < n; i++)
    if( a[i] > max ) max = a[i];

  return( max);
}


/* quick function to determine the maximum size an acceptible substrata can be */
size_t optCV_arrayMaxSubstrata( size_t * a, size_t n, size_t *segments, size_t *acres, size_t * NhSize , double * NhAcres ) {
  size_t i, max;


  max = a[0];

  for( i = 1; i < n; i++)
    if( a[i] > max ) max = a[i];

  return( max);
}


/* quick function to determine the maximum valued index in an array */
size_t optCV_arrayMaxIndexDbl( double * a, size_t n ) {
  size_t i, max;


  max = 0;

  for( i = 1; i < n; i++)
    if( a[i] > a[max] ) max = i;

  return( max);
}


/* this matrix identifies each item with all other items that can have labels assigned */
size_t ** optCV_labelCreateMaster( size_t * label, size_t N, size_t H, size_t NhMax ) {

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
double *** optCV_createContribMatrix( size_t * label, size_t dN, size_t N, size_t H, double * D, size_t ** L, size_t * Nh)
{

  size_t i,j,k,l;
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
      for( k = 0; k < H; k++)
      {
      C[i][j][k] = 0;
      for( l = 0; l < Nh[k]; l++)
        C[i][j][k] += getDist(j,L[k][l],D,i,N); 
      }

  return(C);
}




/* function that creates a variance MDA (matrix) [commodity][strata] */
double ** optCV_createVarMatrix( size_t * label, size_t dN, size_t N, size_t H, double * D, size_t ** L, size_t * Nh, size_t * NhSize)
{

  size_t i,j,k,l;
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
    for( k = 0; k < Nh[j]; k++)
      for( l = k+1; l < Nh[j]; l++) 
        V[i][j] += getDist(L[j][k],L[j][l],D,i,N); 
      V[i][j] = V[i][j] / ((NhSize[j] -1) * NhSize[j]);
    }
    
  return(V);
}


/* function that creates a total MDA (matrix) [commodity][strata] */
double ** optCV_createTotalMatrix( size_t * label, size_t k, size_t dN, size_t N, size_t H, double * x)
{

  size_t i,j,l,offset;
  double ** Total;

  /* allocate memory for C matrix */
  Total = (double **) malloc(sizeof(double * ) * dN );


  for( i = 0; i < dN; i++) Total[i] = (double *) malloc(sizeof(double) * H );
  
 
  /* set initial sum to 0 */
  for( i = 0; i < dN; i++)
    for( j = 0; j < H; j++) Total[i][j] = 0;


  // iterate over commodities
  for( i = 0; i < dN; i++)  {
    offset = N * k;

    // iterate over 
    for( j = 0; j < N; j++)  
      for( l = 0; l < k; l++)  Total[i][ label[j] ] += x[offset*i + j*k + l]/k;

  }

  return(Total);
}


