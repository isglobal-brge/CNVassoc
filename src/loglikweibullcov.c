#include "CNVassoc.h"

// function than returns the loglikelihood value, derivatives and hessian, with covariates.               
SEXP loglikweibullcov(SEXP N, double *par, SEXP P, SEXP Y, SEXP CENS, SEXP WW, SEXP Q, SEXP C){
  int i, j, k, l, *n = INTEGER(N), *p = INTEGER(P), *q = INTEGER(Q), numpar;
  double *y = REAL(Y), *cens = REAL(CENS), *ww = REAL(WW), *c = REAL(C), eta, hhi, ggi, loglik, lambda, scale, 
    ggi2, etacov, x1, x2, censi, yi, ypow, ypowlambda, logy, intercept,
    **hessian, **aux, *score, *dh_theta, *dh_thetai;
    
  numpar = 3 + *q;
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
  loglik = 0;
  scale = par[numpar-1];
  intercept = par[0];
  
  for(i = 0; i<*n; i++){
    ggi=0;
    yi = y[i];
    censi = cens[i];
    ypow = pow(yi,scale);
    logy = log(yi);    
    
    for(k = 0; k<numpar; k++){
      for(l = 0; l<numpar; l++)
        aux[k][l] = 0;
      dh_theta[k] = 0;
    }

    for(j = 0; j< *p; j++){

      etacov = 0;
      for (k = 0; k<*q; k++)
        etacov += par[2+k]*c[i* *q + k];

      eta = intercept + j*par[1] + etacov;
      lambda = exp(eta);
      ypowlambda = ypow*lambda;

      if (censi == 0){
        hhi = (ww[i* *p + j]*exp(-ypowlambda));
      }else{ 
        hhi = (ww[i* *p + j]*lambda*scale*pow(yi,(scale-1))*exp(-ypowlambda));
      }

      ggi += hhi;

      for(k = 0; k<numpar; k++){
        if (k == 0)
          x1 = 1;  // constant
        else if (k == 1)
          x1 = j; // snp
        else 
          x1 = c[i* *q + k-2]; 
        // first derivatives            
        if (censi == 0){
          if (k == numpar-1)
            dh_thetai[k] = -hhi*ypowlambda*logy;
          else  
            dh_thetai[k] = hhi*x1*(-ypow*lambda);
        }else{
          if (k == numpar-1)            
            dh_thetai[k] = -hhi*(ypowlambda*logy-logy-1/scale);
          else
            dh_thetai[k] = hhi*x1*(1-ypow*lambda);
        }
        dh_theta[k] += dh_thetai[k];
      }
      for(k = 0; k<numpar; k++){
        if (k == 0)
          x1 = 1;  // constant
        else if (k == 1)
          x1 = j; // snp
        else 
          x1 = c[i* *q + k-2];                
        for(l = k; l<numpar; l++){
          if (l == 0)
            x2 = 1;
          else if (l == 1)
            x2 = j;
          else 
            x2 = c[i* *q + l-2];      
          // second derivatives     
          if (censi == 0){
            if (k==numpar-1 && l==numpar-1)
              aux[k][l] += -lambda*logy*ypow*(dh_thetai[k]+hhi*logy);
            if (l==numpar-1 && k<numpar-1) 
              aux[k][l] += x1*(-ypowlambda*(dh_thetai[l]+hhi*logy));
            if (l<numpar-1 && k<numpar-1)
              aux[k][l] += x1*(-ypow*lambda*(x2*hhi+dh_thetai[l]));
          }else{
            if (k==numpar-1 && l==numpar-1)
              aux[k][l] += -dh_thetai[k]*(ypow*logy*lambda-logy-1/scale)-hhi*(ypow*pow(logy,2)*lambda+1/pow(scale,2));
            if (k<numpar-1 && l==numpar-1) 
              aux[k][l] += x1*(dh_thetai[l]*(1-ypowlambda)-hhi*ypowlambda*logy);
            if (k<numpar-1 && l<numpar-1)
              aux[k][l] += x1*(dh_thetai[l]*(1-ypowlambda)-hhi*ypowlambda*x2);
          }
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

