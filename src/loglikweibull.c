#include "CNVassoc.h"

// function than returns the loglikelihood value, derivatives and hessian, without covariates.  
SEXP loglikweibull(SEXP N, double *par, SEXP P, SEXP Y, SEXP CENS, SEXP WW){
  int i, j, *n = INTEGER(N), *p = INTEGER(P);
  double *y = REAL(Y), *cens = REAL(CENS), *ww = REAL(WW), eta, hhi, ggi, loglik, lambda, scale, 
    dh_alphai, dh_alpha, dh_beta1i, dh_beta1, dh_scalei, dh_scale, ggi2, snp1, yi, censi, ypow, logy,
    **hessian, *score, aux11, aux12, aux13, aux22, aux23, aux33;
    

  score = (double *)malloc(3*sizeof(double));
  hessian = (double **)malloc(3*sizeof(double *));
  // Initialize hessian and score to 0
  for(i = 0; i<3; i++){
    hessian[i] = (double *)malloc(3*sizeof(double));
    for(j = 0; j<3; j++)
      hessian[i][j] = 0;
    score[i] = 0;
  }

  // Compute 'score', 'loglik' and 'hessian' (IN THE SAME LOOP!)
  loglik = 0;
  scale = par[2];
  
  for(i = 0; i<*n; i++){
    ggi=dh_alpha=dh_beta1=dh_scale=aux11=aux12=aux13=aux22=aux23=aux33=0;
    yi = y[i];
    censi = cens[i];
    ypow = pow(yi,scale);
    logy = log(yi); 
    
    for(j = 0; j< *p; j++){
      snp1 = j;

      eta = par[0] + snp1*par[1];
      lambda = exp(eta);

      if (censi == 0){
        hhi = (ww[i* *p + j]*exp(-lambda*ypow));
      }else{ 
        hhi = (ww[i* *p + j]*lambda*scale*pow(yi,(scale-1))*exp(-lambda*ypow));
      }

      ggi += hhi;

      if (censi == 0){
        dh_alphai = hhi*(-ypow*lambda);
        dh_beta1i = hhi*snp1*(-ypow*lambda);
        dh_scalei = -hhi*lambda*ypow*logy;
      }else{
        dh_alphai = hhi*(1-ypow*lambda);
        dh_beta1i = hhi*snp1*(1-ypow*lambda);
        dh_scalei = -hhi*(lambda*ypow*logy-logy-1/scale);
      }
      dh_alpha += dh_alphai;
      dh_beta1 += dh_beta1i;
      dh_scale += dh_scalei;
      
      if (censi == 0){
        aux11 += 1*(-ypow*lambda*(1*hhi+dh_alphai));
        aux12 += 1*(-ypow*lambda*(snp1*hhi+dh_beta1i));
        aux13 += 1*(-lambda*ypow*(dh_scalei+hhi*logy));
        aux22 += snp1*(-ypow*lambda*(snp1*hhi+dh_beta1i));
        aux23 += snp1*(-lambda*ypow*(dh_scalei+hhi*logy));
        aux33 += -lambda*logy*ypow*(dh_scalei+hhi*logy);
      }else{
        aux11 += 1*(dh_alphai*(1-lambda*ypow)-hhi*lambda*ypow*1);
        aux12 += 1*(dh_beta1i*(1-lambda*ypow)-hhi*lambda*ypow*snp1);
        aux13 += 1*(dh_scalei*(1-lambda*ypow)-hhi*lambda*ypow*logy);
        aux22 += snp1*(dh_beta1i*(1-lambda*ypow)-hhi*lambda*ypow*snp1);
        aux23 += snp1*(dh_scalei*(1-lambda*ypow)-hhi*lambda*ypow*logy);
        aux33 += -dh_scalei*(ypow*logy*lambda-logy-1/scale)-hhi*(ypow*pow(logy,2)*lambda+1/pow(scale,2)); 
      }
    }
                   
    score[0] += dh_alpha/ggi;
    score[1] += dh_beta1/ggi;
    score[2] += dh_scale/ggi;
    
    ggi2 = ggi*ggi;

    hessian[0][0] += (aux11*ggi - dh_alpha*dh_alpha)/ggi2;
    hessian[0][1] += (aux12*ggi - dh_alpha*dh_beta1)/ggi2;
    hessian[0][2] += (aux13*ggi - dh_alpha*dh_scale)/ggi2;
    hessian[1][1] += (aux22*ggi - dh_beta1*dh_beta1)/ggi2;
    hessian[1][2] += (aux23*ggi - dh_beta1*dh_scale)/ggi2;
    hessian[2][2] += (aux33*ggi - dh_scale*dh_scale)/ggi2;

    loglik += log(ggi);
  }

  // RETURN:
  SEXP RES, LOGLIKE, SCORE, HESSIAN;
  double *loglike, *sc, *hess;

  PROTECT(RES = allocVector(VECSXP, 3));
  PROTECT(LOGLIKE = allocVector(REALSXP, 1));
  loglike = REAL(LOGLIKE);
  PROTECT(SCORE = allocVector(REALSXP, 3));
  sc = REAL(SCORE);
  PROTECT(HESSIAN = allocMatrix(REALSXP, 3, 3));
  hess = REAL(HESSIAN);

  *loglike = loglik;
  for(i = 0; i<3; i++){
    sc[i] = score[i];
    for(j = 0; j<3; j++){
      if(j<i)
	hess[j + 3*i] = hessian[j][i];
      else
	hess[j + 3*i] = hessian[i][j];
    }
  }
  UNPROTECT(3);
  SET_VECTOR_ELT(RES, 0, LOGLIKE);
  SET_VECTOR_ELT(RES, 1, SCORE);
  SET_VECTOR_ELT(RES, 2, HESSIAN);

  for(i = 0; i<3; i++)
    free(hessian[i]);
  free(hessian), free(score);

  UNPROTECT(1);
  return(RES);
}

