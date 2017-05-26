#include "CNVassoc.h"

// function than returns the loglikelihood value, derivatives and hessian, without covariates.               
SEXP logliklogistic(SEXP N, double *par, SEXP P, SEXP Y, SEXP WW){
  int i, j, k, l, *n = INTEGER(N), *p = INTEGER(P), numpar;
  double *y = REAL(Y), *ww = REAL(WW), eta, hhi, ggi, loglik, pred, dif1, dif2, 
    ggi2, x1, x2, yi, intercept,
    **hessian, **aux, *score, *dh_theta, *dh_thetai;
    
  numpar = 2;
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
  }

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

    for(j = 0; j< *p; j++){

      eta = intercept + j*par[1];
      pred = 1/(1+exp(-eta));

      hhi = (ww[i* *p + j]*exp(eta*yi))/(1+exp(eta));
      dif1 = yi - pred;
      dif2 = pred - pred*pred;

      ggi += hhi;

      for(k = 0; k<numpar; k++){
        if (k == 0)
          x1 = 1;  // constant
        else
          x1 = j; // snp
        dh_thetai[k] = hhi*x1*dif1;   
        dh_theta[k] += dh_thetai[k];
      }
      for(k = 0; k<numpar; k++){
        if (k == 0)
          x1 = 1;  // constant
        else
          x1 = j; // snp
        for(l = k; l<numpar; l++){
          if (l == 0)
            x2 = 1;
          else
            x2 = j;
          // second derivatives     
          aux[k][l] += dh_thetai[l]*x1*dif1 - hhi*x1*x2*dif2;
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
