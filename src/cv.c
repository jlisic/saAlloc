/* this is a class of functions to efficiently calculate and re-calculate the cv */

#include "cv.h"

/* index,  max_Index - 1,  Description 
 * j    ,  J            ,  Domain 
 * k    ,  K            ,  Characteristic of Interest
 * q    ,  Q            ,  Strata Group
 * h    ,  H            ,  Strata within Strata Group
 * r    ,  R            ,  Observations of each PSU
 *
 * I                       Index of length N
 * N                       Number of PSUs 
 * Domain                  Domain Index   [Strata Group][Commodity]
 * Nh                      Number of PSUs in within strata [Strata]
 *
 * Each within group strata is unique over the entire problem,
 * Strata group is a partitioning of these strata used for domains 
 */




/* function that creates a mean MDA [coi][strata][obs] */
double *** cv_createMeanMatrix( 
    size_t * I,  /* index */
    size_t N,    /* PSUs  */
    size_t K,    /* CoI   */
    size_t H,    /* Strata */
    size_t R,    /* obs */
    double * x,  /* data CoI x I x obs */ 
    size_t * Nh  /* strata size */
  ) {

  size_t r,h,k,i,j;

  /* allocate memory for a matrix of strata by observation for each characteristic*/ 
  double *** mu = (double ***) createMDA(MDA_DOUBLE,3,K,H,R); 
  
  size_t xN = N * K * R;

  /* note that the first K entries in the domain are the commodities */

//  #pragma omp for private(k,r,j,h) 
  for(i = 0; i < xN; i++ ) { 

     k = i / (N*R); /* commodity index */
     r = i % (R);   /* panel index */
     j = i - k*N*R; /* observation index */
     j = j /R;

     h = I[j];       /* stratum */

//     #pragma omp atomic
     mu[k][h][r] += x[i] / (double) Nh[h];  
  }
  
  return(mu);
}



/* function that creates a var MDA [coi][strata][obs] */
double *** cv_createVarMatrix( 
    size_t * I, 
    size_t N, 
    size_t K, 
    size_t H, 
    size_t R, 
    double * x, 
    double *** mu, 
    size_t * Nh 
  ) {

  size_t i; /* element in data structure index */

  size_t k; /* characteristic index */
  size_t j; /* observation index */
  size_t h; /* stratum index */
  size_t r; /* panel index */

  /* allocate memory for a matrix of strata by observation for each characteristic*/ 
  double *** V = (double ***) createMDA(MDA_DOUBLE, 3, K,H,R); 

  size_t xN = N * R * K;
  
//  #pragma omp for private(k,r,j,h) 
  for(i = 0; i < xN; i++ ) { 

     k = i / (N*R); /* commodity index */
     r = i % (R);   /* panel index */
     j = i - k*N*R; /* observation index */
     j = j /R;


     h = I[j];       /* stratum */

//     #pragma omp atomic
     V[k][h][r] += 
        (x[i] - mu[k][h][r]) * (x[i] - mu[k][h][r]) / ( (double) (Nh[h] - 1) ); 
  }

  return(V);
}


/* this is a function that does an online update of the mean and variance matricies */
/* this is a destructive function */
/* only the variance and mean will be updated, not I */
void cv_updateMatrix( 
    size_t * I,     /* this is the current stratification, not the modified one */ 
    size_t N, 
    size_t K, 
    size_t H, 
    size_t R, 
    double * x, 
    double *** mu, 
    double *** var, 
    size_t * Nh,
    size_t moveObs,       /* this is the observation to move */
    size_t moveObsStratum /* this is the stratum to move the observation to */
  ) {


  size_t k,r;
  size_t moveObsOffset;
  double deltaPriorObs;
  double deltaMoveObs;

  /* get moveObs prior staratum */

  size_t priorObsStratum = I[moveObs];

#ifdef CV_DEBUG
printf(" moveFrom: %d\n", (int) priorObsStratum);
printf(" moveTo: %d\n",   (int) moveObsStratum);
#endif

  /* return if there is no change in stratum */
  if( priorObsStratum == moveObsStratum ) return;
  

  for( k = 0; k < K; k++){
    for( r = 0; r < R; r++){

      /******* update the mean **********/
      /* 
        \bar{x}_{n+1} = \bar{x}_n + (x_i -\bar{x}_n)/(n+1) 
        
        \bar{x}_{n-1} = (n \bar{x}_{n} - x_i)/(n-1)  
      */

      /* offset for moveObs */
      moveObsOffset = k*R*N + R*moveObs; 

      deltaPriorObs = x[moveObsOffset + r] - mu[k][priorObsStratum][r]; 
      mu[k][priorObsStratum][r] -= deltaPriorObs / ( (double) Nh[priorObsStratum] - 1 ); 


      /* update the mean of the new stratum */
      deltaMoveObs = x[moveObsOffset + r] - mu[k][moveObsStratum][r]; 
      mu[k][moveObsStratum][r] += deltaMoveObs / (double) (Nh[moveObsStratum] + 1);
      
      
      /******* update the variance **********/
      /*     1   2   3   4   5
       * 1   0 d12 d13 d14 d15
       * 2 d21   0 d23 d24 d25
       * 3 d31 d32   0 d34 d35
       * 4 d41 d42 d43   0 d45
       * 5 d51 d52 d53 d54   0
       *
       * dij = (x_i - x_j)^2
       *
       * S^2_n = ( \sum_{i=1}^n \sum_{j=1}^n dij ) / (n *(n-1))
       *
       * S^2_{n+1} = (n-1)/n S^2_n + ( \sum{i=^n} 2 * di{n+1} )/(n*(n+1)) 
       */

      var[k][moveObsStratum][r] = 
        ((double) Nh[moveObsStratum] - 1) * var[k][moveObsStratum][r] +    /* M2 */ 
        deltaMoveObs *  (x[moveObsOffset + r] - mu[k][moveObsStratum][r]); 

      var[k][moveObsStratum][r] /= (double) Nh[moveObsStratum]; 

      var[k][priorObsStratum][r] = 
        ((double) Nh[priorObsStratum] - 1) * var[k][priorObsStratum][r] -    /* M2 */ 
        deltaPriorObs *  (x[moveObsOffset + r] - mu[k][priorObsStratum][r]); 

      var[k][priorObsStratum][r] /= (double) Nh[priorObsStratum] - 2; 


    }
  }
  
}


/* this is a function that does an online update of the mean and variance matricies */
/* this is a destructive function */
/* only the variance and mean will be updated, not I */
void cv_updateMatrix2( 
    size_t * I,     /* this is the current stratification, not the modified one */ 
    size_t N, 
    size_t K, 
    size_t H, 
    size_t R, 
    double * x, 
    double * x1,
    double * x2,
    double *** mu, 
    double *** mu1, 
    double *** mu2, 
    double *** var, 
    size_t * Nh,
    size_t moveObs,       /* this is the observation to move */
    size_t moveObsStratum /* this is the stratum to move the observation to */
  ) {


  size_t k,r;
  size_t moveObsOffset;
  double deltaPriorObs;
  double deltaMoveObs;
  double deltaPriorObs1;
  double deltaMoveObs1;
  double deltaPriorObs2;
  double deltaMoveObs2;

  /* get moveObs prior staratum */

  size_t priorObsStratum = I[moveObs];

  /* return if there is no change in stratum */
  if( priorObsStratum == moveObsStratum ) return;

  for( k = 0; k < K; k++){
    for( r = 0; r < R; r++){

      /******* update the mean and aux variables **********/

      /* offset for moveObs */
      moveObsOffset = k*R*N + R*moveObs; 

      deltaPriorObs = x[moveObsOffset + r] - mu[k][priorObsStratum][r]; 
      mu[k][priorObsStratum][r] -= deltaPriorObs / ( (double) Nh[priorObsStratum] - 1 ); 


      /* update the mean of the new stratum */
      deltaMoveObs = x[moveObsOffset + r] - mu[k][moveObsStratum][r]; 
      mu[k][moveObsStratum][r] += deltaMoveObs / (double) (Nh[moveObsStratum] + 1);
      

      /* handle adjustments */     
      if( mu1 != NULL) {
        deltaPriorObs1 = x1[moveObsOffset + r] - mu1[k][priorObsStratum][r]; 
        mu1[k][priorObsStratum][r] -= deltaPriorObs1 / ( (double) Nh[priorObsStratum] - 1 ); 
        deltaMoveObs1 = x1[moveObsOffset + r] - mu1[k][moveObsStratum][r]; 
        mu1[k][moveObsStratum][r] += deltaMoveObs1 / (double) (Nh[moveObsStratum] + 1);
      }
      
      if( mu2 != NULL) {
        deltaPriorObs2 = x2[moveObsOffset + r] - mu2[k][priorObsStratum][r]; 
        mu2[k][priorObsStratum][r] -= deltaPriorObs2 / ( (double) Nh[priorObsStratum] - 1 ); 
        deltaMoveObs2 = x2[moveObsOffset + r] - mu2[k][moveObsStratum][r]; 
        mu2[k][moveObsStratum][r] += deltaMoveObs2 / (double) (Nh[moveObsStratum] + 1);
      }

      /******* update the variance **********/
      var[k][moveObsStratum][r] = 
        ((double) Nh[moveObsStratum] - 1) * var[k][moveObsStratum][r] +    /* M2 */ 
        deltaMoveObs *  (x[moveObsOffset + r] - mu[k][moveObsStratum][r]); 

      var[k][moveObsStratum][r] /= (double) Nh[moveObsStratum]; 

      var[k][priorObsStratum][r] = 
        ((double) Nh[priorObsStratum] - 1) * var[k][priorObsStratum][r] -    /* M2 */ 
        deltaPriorObs *  (x[moveObsOffset + r] - mu[k][priorObsStratum][r]); 

      var[k][priorObsStratum][r] /= (double) Nh[priorObsStratum] - 2; 


    }
  }
  
}



/***************** functions for totals ********************************/


/**********************************


Domain[h][k][j]

Domain[h] = 
    0  1  2  ... K-1 ... J-1 
0   1  0  0      0 
1   0  1  0      0 
2   0  0  1      0 
...
K-1 0  0  0 ...  1

***********************************/

/* create an array of totals */
double ** cv_createTotalMatrix( 
    size_t N, 
    size_t K, 
    size_t H,
    size_t R, 
    size_t J,
    size_t *** Domain,
    double ***mu, 
    size_t * Nh
    ) {

  size_t k,h,j,r;

  /* created with calloc, so zero'ed out */
  double ** muArray = (double **) createMDA(MDA_DOUBLE,2,J,R); 

  for( k = 0; k < K; k++) {
    for( h = 0; h < H; h++) {
      for( j = 0; j < J; j++){ 
        if( Domain[h][k][j] ) {
          for( r = 0; r < R; r++ ) {
            muArray[j][r] += Nh[h] * mu[k][h][r];
          }
        }
      }
    }
  }

  return(muArray);
}


/* apply domains to variable */
double *** cv_createDomainMDA( 
    size_t N, 
    size_t K, 
    size_t H,
    size_t R, 
    size_t J,
    size_t *** Domain,
    double ***x,
    double *** locationAdj,   /* adjustments to location of x [k][h] */
    double *** scaleAdj       /* adjustments to scale of x[k][h] */
    ) {

  size_t k,h,j,r;

  /* created with calloc, so zero'ed out */
  double *** xArray = (double ***) createMDA(MDA_DOUBLE,3,J,H,R); 

  /* messy but efficient adjustment */
  if( (locationAdj == NULL) & (scaleAdj == NULL) ) {
  
    for( k = 0; k < K; k++) {
      for( h = 0; h < H; h++) {
        for( j = 0; j < J; j++){ 
          if( Domain[h][k][j] ) {
            for( r = 0; r < R; r++ ) {
              xArray[j][h][r] = x[k][h][r];
            }
          }
        }
      }
    }
    return(xArray);
  }
  
  if( locationAdj == NULL ) {
  
    for( k = 0; k < K; k++) {
      for( h = 0; h < H; h++) {
        for( j = 0; j < J; j++){ 
          if( Domain[h][k][j] ) {
            for( r = 0; r < R; r++ ) {
              xArray[j][h][r] = x[k][h][r] * scaleAdj[k][h][r];
            }
          }
        }
      }
    }
    return(xArray);
  }
 
  if( scaleAdj == NULL ) {
  
    for( k = 0; k < K; k++) {
      for( h = 0; h < H; h++) {
        for( j = 0; j < J; j++){ 
          if( Domain[h][k][j] ) {
            for( r = 0; r < R; r++ ) {
              xArray[j][h][r] = x[k][h][r]  + locationAdj[k][h][r];
            }
          }
        }
      }
    }
    return(xArray);
  }
  

  for( k = 0; k < K; k++) {
    for( h = 0; h < H; h++) {
      for( j = 0; j < J; j++){ 
        if( Domain[h][k][j] ) {
          for( r = 0; r < R; r++ ) {
            xArray[j][h][r] = x[k][h][r] * scaleAdj[k][h][r] + locationAdj[k][h][r];
          }
        }
      }
    }
  }

  return(xArray);
}



/* cv calc */
double * cv_calcCV(
    double * cv,   /* if cv = NULL create it, otherwise handle in place */
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
    size_t fpc
  ) {

  double *** varDomain;
  size_t j,h,r;
  double RDouble = (double) R;
  double varSum[R];

  /* get domain variance */ 
  varDomain =  cv_createDomainMDA( N, K, H, R, J, Domain, var, locationAdj, scaleAdj);

  if( cv == NULL) cv = calloc( J, sizeof(double) );

  if( fpc == 1 ) {
    for( j = 0; j < J; j++) {
      cv[j] = 0;
      for( r = 0; r < R; r++ ) {
        varSum[r] = 0;
        for( h = 0; h < H; h++) {
          varSum[r] +=  (double) (Nh[h] * Nh[h]) * (1- (double) nh[h]/Nh[h]) /nh[h] * varDomain[j][h][r];
        }
        cv[j] += sqrt(varSum[r])/( Total[j][r] * RDouble);
      }
    }
  } else {
    for( j = 0; j < J; j++) {
      cv[j] = 0;
      for( r = 0; r < R; r++ ) {
        varSum[r] = 0;
        for( h = 0; h < H; h++) {
          varSum[r] +=  (double) (Nh[h] * Nh[h]) * (1- (double) nh[h]/Nh[h]) /nh[h] * varDomain[j][h][r];
        }
        cv[j] += sqrt(varSum[r])/( Total[j][r] * RDouble);
      }
    }
  }
  
  deleteMDA( (void *) varDomain, 3, J, H); 

  return cv;
}


/* objective function */
double cv_objectiveFunction( 
    double * cv,   /* if cv = NULL create it, otherwise handle in place */
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
    size_t evaluateOnly, // option to not construct CV, under this condition CV cannot be null 
    size_t fpc
  ) {


  double result = 0;
  size_t j;
  double delta = 0;

  /* get the CV */
  if( evaluateOnly != 1 ) {
    cv = cv_calcCV( cv, N, K, H, R, J, Domain, var, Nh, nh, Total, locationAdj, scaleAdj,fpc); 
  }
    

  for( j = 0; j < J; j++) {


    if( Target != NULL) delta = fabs(cv[j]) - fabs(Target[j]);
 
    /* apply the penalty function */ 
    if( (delta <= 0) | (penalty == NULL) | (Target == NULL) ) {
      result+= pow(cv[j], p); 
    } else {
      result+= pow(delta, p) * penalty[j] + pow(cv[j],p); 
    }
  }

  return( result );
}



/* objective function */
/* this is an altenative cv function - this is designed for the SA algorithm */
/* this is a function to compare a prior objective function to a new one  */
double cv_objectiveFunctionCompare( 
    double * cv,         //if cv = NULL create it, otherwise handle in place 
    double * cvPrior,   
    double * obj,
    double * objPrior,
    size_t N, 
    size_t K, 
    size_t H, 
    size_t R, 
    size_t J, 
    size_t *** Domain,     // only used if evaluateOnly != 1
    double *** var,        // only used if evaluateOnly != 1 
    size_t * Nh,
    double * nh,
    double ** Total,       // only used if evaluateOnly != 1
    double *** locationAdj, // only used if evaluateOnly != 1
    double *** scaleAdj,    // only used if evaluateOnly != 1
    double * Target,
    double * penalty,
    double p,
    size_t evaluateOnly, // option to not construct CV, under this condition CV cannot be null 
    size_t preserveSatisfied,
    size_t fpc
  ) {
  
  double result = 0;
  double resultPrior = 0;
  double result_penalty = 0;
  double resultPrior_penalty=0;
  size_t j;
  double delta = 0;
  double deltaPrior = 0;
  size_t failFlag = 0;  // value set to 1 if there is a violation of preserve satisfied


  /* get the CV */
  if( evaluateOnly != 1 ) {
    cv = cv_calcCV( cv, N, K, H, R, J, Domain, var, Nh, nh, Total, locationAdj, scaleAdj,fpc ); 
  }

  for( j = 0; j < J; j++) {

    if( Target != NULL) {
      delta = cv[j] - Target[j];
      deltaPrior = cvPrior[j] - Target[j];

      // handle the condition on preserve satisfied
      if( 
          (deltaPrior <= 0) &
          (delta > 0) & 
          (preserveSatisfied == 1) 
        ) failFlag = 1; 
    }


    // apply the penalty function 
    if( (deltaPrior <= 0) | (penalty == NULL) | (Target == NULL) ) {
      resultPrior+= pow(cvPrior[j], p); 
    } else {
      resultPrior+= pow(cvPrior[j], p); 
      resultPrior_penalty+= deltaPrior * penalty[j]; 
    }
 
    // apply the penalty function 
    if( (delta <= 0) | (penalty == NULL) | (Target == NULL) ) {
      result+= pow(cv[j], p); 
    } else {
      result+= pow(cv[j], p); 
      result_penalty+= delta * penalty[j]; 
    }
  }

  result=sqrt(result) + result_penalty;
  resultPrior = sqrt(resultPrior) + resultPrior_penalty;

  if( obj != NULL ) *obj = result;
  if( objPrior != NULL ) *objPrior = resultPrior;
        
  
  if( failFlag == 1) return(INFINITY);

//printf("(%4.12f) (%4.12f)", result, resultPrior);
printf("%4.12f ", result);
printf("\n");
  return( result - resultPrior );
}  
  
  
  
  
  
