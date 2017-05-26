#include "CNVassoc.h"

SEXP loglikinterlogisticcov(SEXP N, double *par, SEXP Y, SEXP WW, SEXP Q, SEXP C){
  int i, j, k, l, *n = INTEGER(N), p = 3, p2 = p*p, snp1, snp2, snp12, *q = INTEGER(Q), numpar;
  double *y = REAL(Y), *ww = REAL(WW), *c = REAL(C), intercept, yi, x1, x2, etacov, eta, pred, hhi, ggi, ggi2, loglik, dif1, dif2,
    **hessian, **aux, *score, *dh_theta, *dh_thetai;
    
  numpar = 4 + *q;
  score = (double *)malloc(numpar*sizeof(double));
  dh_theta = (double *)malloc(numpar*sizeof(double));
  dh_thetai = (double *)malloc(numpar*sizeof(double));  
  hessian = (double **)malloc(numpar*sizeof(double *));
  aux = (double **)malloc(numpar*sizeof(double *));
  // Initialize hessian and score to 0
  for(i = 0; i<numpar; i++){
    hessian[i] = (double *)malloc(numpar*sizeof(double));
    aux[i] = (double *)malloc(numpar*sizeof(double));
    for(j = 0; j<numpar; j++){
      hessian[i][j] = 0;
      aux[i][j] = 0;
    }
    score[i] = 0;
    dh_theta[i] = 0;
    dh_thetai[i] = 0;
  }

  // Compute 'score', 'loglik' and 'hessian' (IN THE SAME LOOP!)
  // Compute 'score', 'loglik' and 'hessian' (IN THE SAME LOOP!)
  loglik = 0;
  intercept = par[0];  
  for(i = 0; i<*n; i++){

    ggi=0;
    yi = y[i];
    
    for(k = 0; k<numpar; k++){
      for(l = 0; l<numpar; l++)
        aux[k][l] = 0;
      dh_theta[k] = 0;
    }
    
    etacov = 0;
    for (k = 0; k<*q; k++)
      etacov += par[4+k]*c[i* *q + k];

    for(j = 0; j<p2; j++){

      snp1 = j%p;
      snp2 = floor(j/p);
      snp12 = snp1*snp2;    

      eta = intercept + snp1*par[1] + snp2*par[2] + snp12*par[3] + etacov;
      pred = 1/(1+exp(-eta));
      dif1 = yi - pred;
      dif2 = pred - pred*pred;

      hhi = (ww[i*p2 + j]*exp(eta*yi))/(1+exp(eta));
      ggi += hhi;


      for(k = 0; k<numpar; k++){
        if (k == 0)
          x1 = 1;  // constant
        else if (k == 1)
          x1 = snp1; // snp1
        else if (k ==2)
          x1 = snp2; // snp2
        else if (k ==3)
          x1 = snp12; // snp1*snp2
        else 
          x1 = c[i* *q + k-4];  //covariate 
        // first derivatives            
        dh_thetai[k] = hhi*x1*dif1;
        dh_theta[k] += dh_thetai[k];
      }
      
      for(k = 0; k<numpar; k++){
        if (k == 0)
          x1 = 1;  // constant
        else if (k == 1)
          x1 = snp1; // snp1
        else if (k == 2)
          x1 = snp2; // snp2
        else if (k == 3)
          x1 = snp12; // snp1*snp2
        else 
          x1 = c[i* *q + k-4]; // covariate               
        for(l = k; l<numpar; l++){
          if (l == 0)
            x2 = 1;
          else if (l == 1)
            x2 = snp1; // snp1
          else if (l == 2)
            x2 = snp2; // snp2
          else if (l == 3)
            x2 = snp12; // snp1*snp2
          else 
            x2 = c[i* *q + l-4];  //covariates
          // second derivatives     
          aux[k][l] += dh_thetai[l]*x1*dif1-hhi*x1*x2*dif2;
        }
      }

    }

    ggi2 = ggi*ggi;

    loglik += log(ggi);
                   
    for(k = 0; k<numpar; k++){
      for (l = 0; l<numpar; l++)
        hessian[k][l] += (aux[k][l]*ggi - dh_theta[k]*dh_theta[l])/ggi2;
      score[k] += dh_theta[k]/ggi;
    }

  }

  // RETURN:
  SEXP RES, LOGLIKE, SCORE, HESSIAN;
  double *loglike, *sc, *hess;

  PROTECT(RES = allocVector(VECSXP, numpar));
  PROTECT(LOGLIKE = allocVector(REALSXP, 1));
  loglike = REAL(LOGLIKE);
  PROTECT(SCORE = allocVector(REALSXP, numpar));
  sc = REAL(SCORE);
  PROTECT(HESSIAN = allocMatrix(REALSXP, numpar, numpar));
  hess = REAL(HESSIAN);

  *loglike = loglik;
  for(i = 0; i<numpar; i++){
    sc[i] = score[i];
    for(j = 0; j<numpar; j++){
      if(j<i)
	hess[j + numpar*i] = hessian[j][i];
      else
	hess[j + numpar*i] = hessian[i][j];
    }
  }
  UNPROTECT(3);
  SET_VECTOR_ELT(RES, 0, LOGLIKE);
  SET_VECTOR_ELT(RES, 1, SCORE);
  SET_VECTOR_ELT(RES, 2, HESSIAN);

  for(i = 0; i<numpar; i++)
    free(hessian[i]);
  free(hessian), free(score);

  UNPROTECT(1);
  return(RES);
}

