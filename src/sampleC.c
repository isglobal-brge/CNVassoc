#include "CNVassoc.h"


SEXP sampleC(SEXP N, SEXP X){
  int *n = INTEGER(N), i, *res;
  double *x = REAL(X), u;
  SEXP RES;
  
  PROTECT(RES = allocVector(INTSXP, *n));
  res = INTEGER(RES);

  GetRNGstate();
  for(i = 0; i<*n; i++){
    u = runif(0,1);
    if(u<x[3*i])
      res[i] = 0;
    else if(u<(x[3*i]+x[3*i + 1]))
      res[i] = 1;
    else
      res[i] = 2;
  }
  PutRNGstate();
  
  UNPROTECT(1);
  return(RES);
}
