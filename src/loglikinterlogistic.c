#include "CNVassoc.h"

SEXP loglikinterlogistic(SEXP N, double *par, SEXP Y, SEXP WW){
  int i, j, *n = INTEGER(N), p = 3, p2 = p*p, snp1, snp2, snp12; //P CAN BE MODIFIED TO ACCEPT MORE PROBABILITIES 
  double *y = REAL(Y), *ww = REAL(WW), *eta, *pred, hhi, ggi, loglik, dif1, dif2,
    dh_alphai, dh_alpha, dh_beta1i, dh_beta1, dh_beta2i, dh_beta2, dh_beta12i, dh_beta12, ggi2,
    **hessian, *score, aux11, aux12, aux13, aux14, aux22, aux23, aux24, aux33, aux34, aux44;
    
  eta = (double *)malloc(p2*sizeof(double));
  pred = (double *)malloc(p2*sizeof(double));
  score = (double *)malloc(4*sizeof(double));
  hessian = (double **)malloc(4*sizeof(double *));
  // Initialize hessian and score to 0
  for(i = 0; i<4; i++){
    hessian[i] = (double *)malloc(4*sizeof(double));
    for(j = 0; j<4; j++)
      hessian[i][j] = 0;
    score[i] = 0;
  }

  // Compute 'eta' and 'pred'
  for(i = 0; i<p; i++)
    for(j = 0; j<p; j++){
      eta[i*p + j] = par[0] + j*par[1] + i*par[2] + i*j*par[3];
      pred[i*p + j] = 1/(1+exp(-eta[i*p + j]));
    }

  // Compute 'score', 'loglik' and 'hessian' (IN THE SAME LOOP!)
  loglik = 0;
  for(i = 0; i<*n; i++){
    ggi=dh_alpha=dh_beta1=dh_beta2=dh_beta12=aux11=aux12=aux13=aux14=aux22=aux23=aux24=aux33=aux34=aux44=0;
    for(j = 0; j<p2; j++){
      snp1 = j%p;
      snp2 = floor(j/p);
      snp12 = snp1*snp2;
      dif1 = y[i] - pred[j];
      dif2 = pred[j] - pred[j]*pred[j];

      hhi = (ww[i*p2 + j]*exp(eta[j]*y[i]))/(1+exp(eta[j]));
      ggi += hhi;

      dh_alphai = hhi*dif1;
      dh_alpha += dh_alphai;
      dh_beta1i = hhi*snp1*dif1;
      dh_beta1 += dh_beta1i;
      dh_beta2i = hhi*snp2*dif1;
      dh_beta2 += dh_beta2i;
      dh_beta12i = hhi*snp12*dif1;
      dh_beta12 += dh_beta12i;
      
      aux11 += dh_alphai*dif1 - hhi*dif2;
      aux12 += dh_beta1i*dif1 - hhi*snp1*dif2;
      aux13 += dh_beta2i*dif1 - hhi*snp2*dif2;
      aux14 += dh_beta12i*dif1 - hhi*snp12*dif2;
      aux22 += dh_beta1i*snp1*dif1 - hhi*snp1*snp1*dif2;
      aux23 += dh_beta2i*snp1*dif1 - hhi*snp1*snp2*dif2;
      aux24 += dh_beta12i*snp1*dif1 - hhi*snp1*snp12*dif2;
      aux33 += dh_beta2i*snp2*dif1 - hhi*snp2*snp2*dif2;
      aux34 += dh_beta12i*snp2*dif1 - hhi*snp2*snp12*dif2;
      aux44 += dh_beta12i*snp12*dif1 - hhi*snp12*snp12*dif2;
    }
    score[0] += dh_alpha/ggi;
    score[1] += dh_beta1/ggi;
    score[2] += dh_beta2/ggi;
    score[3] += dh_beta12/ggi;

    ggi2 = ggi*ggi;

    hessian[0][0] += (aux11*ggi - dh_alpha*dh_alpha)/ggi2;
    hessian[0][1] += (aux12*ggi - dh_alpha*dh_beta1)/ggi2;
    hessian[0][2] += (aux13*ggi - dh_alpha*dh_beta2)/ggi2;
    hessian[0][3] += (aux14*ggi - dh_alpha*dh_beta12)/ggi2;
    hessian[1][1] += (aux22*ggi - dh_beta1*dh_beta1)/ggi2;
    hessian[1][2] += (aux23*ggi - dh_beta1*dh_beta2)/ggi2;
    hessian[1][3] += (aux24*ggi - dh_beta1*dh_beta12)/ggi2,
    hessian[2][2] += (aux33*ggi - dh_beta2*dh_beta2)/ggi2;
    hessian[2][3] += (aux34*ggi - dh_beta2*dh_beta12)/ggi2;
    hessian[3][3] += (aux44*ggi - dh_beta12*dh_beta12)/ggi2;

    loglik += log(ggi);
  }

  // RETURN:
  SEXP RES, LOGLIKE, SCORE, HESSIAN;
  double *loglike, *sc, *hess;

  PROTECT(RES = allocVector(VECSXP, 3));
  PROTECT(LOGLIKE = allocVector(REALSXP, 1));
  loglike = REAL(LOGLIKE);
  PROTECT(SCORE = allocVector(REALSXP, 4));
  sc = REAL(SCORE);
  PROTECT(HESSIAN = allocMatrix(REALSXP, 4, 4));
  hess = REAL(HESSIAN);

  *loglike = loglik;
  for(i = 0; i<4; i++){
    sc[i] = score[i];
    for(j = 0; j<4; j++){
      if(j<i)
	hess[j + 4*i] = hessian[j][i];
      else
	hess[j + 4*i] = hessian[i][j];
    }
  }
  UNPROTECT(3);
  SET_VECTOR_ELT(RES, 0, LOGLIKE);
  SET_VECTOR_ELT(RES, 1, SCORE);
  SET_VECTOR_ELT(RES, 2, HESSIAN);

  /*  for(i = 0; i<4; i++){
    for(j = 0; j<4; j++)
      Rprintf("%g ", hessian[i][j]);
    Rprintf("\n");
    }*/
  
  for(i = 0; i<4; i++)
    free(hessian[i]);
  free(hessian), free(eta), free(pred), free(score);

  UNPROTECT(1);
  return(RES);
}


