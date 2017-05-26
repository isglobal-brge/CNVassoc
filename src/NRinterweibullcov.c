#include "CNVassoc.h"

SEXP NRinterweibullcov(SEXP N, SEXP PARAM, SEXP Y, SEXP CENS, SEXP WW, SEXP Q, SEXP C, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE){
  int *max_iter = INTEGER(MAX_ITER), *q = INTEGER(Q), iter = 1, i, j, p, *verbose = INTEGER(VERBOSE);
  double *tol = REAL(TOL), error, *parIni = REAL(PARAM), *l, *s, *h, *par, *par_old, aux, *hessian;
  SEXP DERIVATIVES, L, S, H, H2;
  
  
  p = 5 + *q; // number of parameters
  
  double *id;
  id = (double *)malloc((p*p)*sizeof(double));
  par = (double *)malloc(p*sizeof(double));
  par_old = (double *)malloc(p*sizeof(double));
  for(i=0;i<p;i++)
    par[i] = parIni[i];

  // VARIABLES FOR INVERSE COMPUTATION:
  int *ipiv;
  ipiv = (int *)malloc(p*sizeof(int));
  int nr = p;
  int *pnr;
  pnr = &nr;
  *pnr = nr;
  int inf = 0;
  int *pinf;
  pinf = &inf;
  *pinf = inf;

  error = *tol + 0.5;
  // ITERATIVE PROCESS:
  while((error>*tol) && (iter<*max_iter)){
  
    // STORE OLD PARAMETERS:
    for(i=0; i<p; i++)
      par_old[i] = par[i];

    // SCORE AND HESSIAN COMPUTATION:
    DERIVATIVES = loglikinterweibullcov(N, par, Y, CENS, WW, Q, C);

    L = VECTOR_ELT(DERIVATIVES, 0);
    S = VECTOR_ELT(DERIVATIVES, 1);
    H = VECTOR_ELT(DERIVATIVES, 2);

    // STORE HESSIAN TO RETURN:
    PROTECT(H2 = allocMatrix(REALSXP, p, p));
    Rf_copyMatrix(H2, H, FALSE);
    UNPROTECT(1);
    l = REAL(L);
    s = REAL(S);
    h = REAL(H);

    // IDENTITY MATRIX:
    for(i = 0; i<p; i++){
      for(j = 0; j<p; j++){
      	if(i==j)
      	  id[i*p + j] = 1;
        else
      	  id[i*p + j] = 0;
      }
    }

    // COMPUTE INVERSE OF HESSIAN:
    F77_NAME(dgesv)(&p, pnr, h, &p, ipiv, id, &p, pinf);
    
    // UPDATE PARAMETERS AND COMPUTE ERROR:
    error = -1;
    for(i = 0; i<p; i++){
      for(j = 0; j<p; j++)
      	par[i] -= id[p*i + j]*s[j];
      aux = par[i] - par_old[i];
      if(aux>=0){
      	if(aux>error)
          error = aux;
      } else{
      	if(-aux>error)
        error = -aux;
      }
    }
    
    // VERBOSE:
    if(*verbose == 1){
      Rprintf("... iteration ... %d\n", iter);
      Rprintf("likelihood = ");
      Rprintf("%g ", l);
      Rprintf("\nparam = ");
      for(i=0; i<p; i++)
       	Rprintf("%g ", par[i]);
      Rprintf("\nscore = ");
      for(i=0; i<p; i++)
	     Rprintf("%g ", s[i]);
      hessian = REAL(H2);
      Rprintf("\nHessian\n:");
      for(i=0; i<p; i++){
	       for(j=0; j<p; j++)
	         Rprintf("%g ", hessian[p*i + j]);
	       Rprintf("\n");
      }
      Rprintf("error = %g\n\n", error);
    }
    iter++;
  }

  // RETURN:
  SEXP RES, RESPARAM, ITEROUT;
  int *iterOut;
  double *resparam;

  PROTECT(ITEROUT = allocVector(INTSXP, 1));
  iterOut = INTEGER(ITEROUT);
  *iterOut = iter;
  UNPROTECT(1);

  PROTECT(RESPARAM = allocVector(REALSXP, p));
  resparam = REAL(RESPARAM);
  for(i = 0; i<p; i++)
    resparam[i] = par[i];
  UNPROTECT(1);
  
  PROTECT(RES = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(RES, 0, RESPARAM);
  SET_VECTOR_ELT(RES, 1, L);
  SET_VECTOR_ELT(RES, 2, S);
  SET_VECTOR_ELT(RES, 3, H2);
  SET_VECTOR_ELT(RES, 4, ITEROUT);     
  UNPROTECT(1);
  
  // FREE ALLOCATED MEMORY:
  free(id), free(par_old), free(ipiv);

  return RES;
}

