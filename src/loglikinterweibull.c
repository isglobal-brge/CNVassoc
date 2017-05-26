#include "CNVassoc.h"

SEXP loglikinterweibull(SEXP N, double *par, SEXP Y, SEXP CENS, SEXP WW){

  int i, j, *n = INTEGER(N), p = 3, p2 = p*p, snp1, snp2, snp12; //P CAN BE MODIFIED TO ACCEPT MORE PROBABILITIES 
  double *y = REAL(Y), *cens = REAL(CENS), *ww = REAL(WW), eta, lambda, scale, yi, censi, ypow, logy, ypowlambda, hhi, ggi, loglik,
    dh_alphai, dh_alpha, dh_beta1i, dh_beta1, dh_beta2i, dh_beta2, dh_beta12i, dh_beta12, dh_scalei, dh_scale, ggi2,
    **hessian, *score, aux11, aux12, aux13, aux14, aux15, aux22, aux23, aux24, aux25, aux33, aux34, aux35, aux44, aux45, aux55;
    
  score = (double *)malloc(5*sizeof(double));
  hessian = (double **)malloc(5*sizeof(double *));
  // Initialize hessian and score to 0
  for(i = 0; i<5; i++){
    hessian[i] = (double *)malloc(5*sizeof(double));
    for(j = 0; j<5; j++)
      hessian[i][j] = 0;
    score[i] = 0;
  }

  // Compute 'score', 'loglik' and 'hessian' (IN THE SAME LOOP!)

  loglik = 0;

  for(i = 0; i<*n; i++){
  
    yi = y[i];
    censi = cens[i];
    ypow = pow(yi,scale);
    logy = log(yi);

    ggi=dh_alpha=dh_beta1=dh_beta2=dh_beta12=dh_scale=aux11=aux12=aux13=aux14=aux15=aux22=aux23=aux24=aux25=aux33=aux34=aux35=aux44=aux45=aux55=0;

    for(j = 0; j<p2; j++){
    
      snp1 = j%p;
      snp2 = floor(j/p);
      snp12 = snp1*snp2;
      eta = par[0] + snp1*par[1] + snp2*par[2] + snp12*par[3];
      scale = par[4];
      lambda = exp(eta);
      ypowlambda = ypow*lambda;

      if (censi == 0){
        hhi = ww[i*p2 + j]*exp(-ypowlambda);
      }else{ 
        hhi = ww[i*p2 + j]*lambda*scale*ypow*exp(-ypowlambda)/yi;
      }      
      ggi += hhi;
      
      if (censi == 0){
        dh_alphai = hhi*(-ypowlambda);
        dh_beta1i = hhi*snp1*(-ypowlambda);
        dh_beta2i = hhi*snp2*(-ypowlambda);
        dh_beta12i = hhi*snp12*(-ypowlambda);
        dh_scalei = -hhi*ypowlambda*logy;
      }else{
        dh_alphai = hhi*(1-ypowlambda);
        dh_beta1i = hhi*snp1*(1-ypowlambda);
        dh_beta2i = hhi*snp2*(1-ypowlambda);
        dh_beta12i = hhi*snp12*(1-ypowlambda);        
        dh_scalei = -hhi*(ypowlambda*logy-logy-1/scale);
      }
      dh_alpha += dh_alphai;
      dh_beta1 += dh_beta1i;
      dh_beta2 += dh_beta2i;
      dh_beta12 += dh_beta12i;
      dh_scale += dh_scalei;
      
      if (censi == 0){
        aux11 += (-ypowlambda*(hhi+dh_alphai));
        aux12 += (-ypowlambda*(snp1*hhi+dh_beta1i));
        aux13 += (-ypowlambda*(snp2*hhi+dh_beta2i));
        aux14 += (-ypowlambda*(snp12*hhi+dh_beta12i));        
        aux15 += (-ypowlambda*(dh_scalei+hhi*logy));
        aux22 += snp1*(-ypowlambda*(snp1*hhi+dh_beta1i));
        aux23 += snp1*(-ypowlambda*(snp2*hhi+dh_beta2i));
        aux24 += snp1*(-ypowlambda*(snp12*hhi+dh_beta12i));        
        aux25 += snp1*(-ypowlambda*(dh_scalei+hhi*logy));
        aux33 += snp2*(-ypowlambda*(snp2*hhi+dh_beta2i));
        aux34 += snp2*(-ypowlambda*(snp12*hhi+dh_beta12i));        
        aux35 += snp2*(-ypowlambda*(dh_scalei+hhi*logy));
        aux44 += snp12*(-ypowlambda*(snp12*hhi+dh_beta12i));        
        aux45 += snp12*(-ypowlambda*(dh_scalei+hhi*logy));
        aux55 += -lambda*logy*ypow*(dh_scalei+hhi*logy);
      }else{
        aux11 += (dh_alphai*(1-ypowlambda)-hhi*ypowlambda);
        aux12 += (dh_beta1i*(1-ypowlambda)-hhi*ypowlambda*snp1);
        aux13 += (dh_beta2i*(1-ypowlambda)-hhi*ypowlambda*snp2);
        aux14 += (dh_beta12i*(1-ypowlambda)-hhi*ypowlambda*snp12);        
        aux15 += (dh_scalei*(1-ypowlambda)-hhi*ypowlambda*logy);
        aux22 += snp1*(dh_beta1i*(1-ypowlambda)-hhi*ypowlambda*snp1);
        aux23 += snp1*(dh_beta2i*(1-ypowlambda)-hhi*ypowlambda*snp2);
        aux24 += snp1*(dh_beta12i*(1-ypowlambda)-hhi*ypowlambda*snp12);
        aux25 += snp1*(dh_scalei*(1-ypowlambda)-hhi*ypowlambda*logy);
        aux33 += snp2*(dh_beta2i*(1-ypowlambda)-hhi*ypowlambda*snp2);
        aux34 += snp2*(dh_beta12i*(1-ypowlambda)-hhi*ypowlambda*snp12);
        aux35 += snp2*(dh_scalei*(1-ypowlambda)-hhi*ypowlambda*logy);
        aux44 += snp12*(dh_beta12i*(1-ypowlambda)-hhi*ypowlambda*snp12);
        aux45 += snp12*(dh_scalei*(1-ypowlambda)-hhi*ypowlambda*logy);
        aux55 += -dh_scalei*(ypow*logy*lambda-logy-1/scale)-hhi*(ypow*pow(logy,2)*lambda+1/pow(scale,2)); 
      }

    }
    score[0] += dh_alpha/ggi;
    score[1] += dh_beta1/ggi;
    score[2] += dh_beta2/ggi;
    score[3] += dh_beta12/ggi;
    score[4] += dh_scale/ggi;    

    ggi2 = ggi*ggi;

    hessian[0][0] += (aux11*ggi - dh_alpha*dh_alpha)/ggi2;
    hessian[0][1] += (aux12*ggi - dh_alpha*dh_beta1)/ggi2;
    hessian[0][2] += (aux13*ggi - dh_alpha*dh_beta2)/ggi2;
    hessian[0][3] += (aux14*ggi - dh_alpha*dh_beta12)/ggi2;
    hessian[0][4] += (aux15*ggi - dh_alpha*dh_scale)/ggi2;    
    hessian[1][1] += (aux22*ggi - dh_beta1*dh_beta1)/ggi2;
    hessian[1][2] += (aux23*ggi - dh_beta1*dh_beta2)/ggi2;
    hessian[1][3] += (aux24*ggi - dh_beta1*dh_beta12)/ggi2,
    hessian[1][4] += (aux25*ggi - dh_beta1*dh_scale)/ggi2,
    hessian[2][2] += (aux33*ggi - dh_beta2*dh_beta2)/ggi2;
    hessian[2][3] += (aux34*ggi - dh_beta2*dh_beta12)/ggi2;
    hessian[2][4] += (aux35*ggi - dh_beta2*dh_scale)/ggi2;
    hessian[3][3] += (aux44*ggi - dh_beta12*dh_beta12)/ggi2;
    hessian[3][4] += (aux45*ggi - dh_beta12*dh_scale)/ggi2;    
    hessian[4][4] += (aux55*ggi - dh_scale*dh_scale)/ggi2;    

    loglik += log(ggi);
  }

  // RETURN:
  SEXP RES, LOGLIKE, SCORE, HESSIAN;
  double *loglike, *sc, *hess;

  PROTECT(RES = allocVector(VECSXP, 3));
  PROTECT(LOGLIKE = allocVector(REALSXP, 1));
  loglike = REAL(LOGLIKE);
  PROTECT(SCORE = allocVector(REALSXP, 5));
  sc = REAL(SCORE);
  PROTECT(HESSIAN = allocMatrix(REALSXP, 5, 5));
  hess = REAL(HESSIAN);

  *loglike = loglik;
  for(i = 0; i<5; i++){
    sc[i] = score[i];
    for(j = 0; j<5; j++){
      if(j<i)
	hess[j + 5*i] = hessian[j][i];
      else
	hess[j + 5*i] = hessian[i][j];
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
  
  for(i = 0; i<5; i++)
    free(hessian[i]);
  free(hessian), free(score);

  UNPROTECT(1);
  return(RES);
}

