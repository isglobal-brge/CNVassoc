#include "CNVassoc.h"

SEXP maxC(SEXP N, SEXP X){
  int *n = INTEGER(N), i, *res;
  double *x = REAL(X);
  SEXP RES;
  
  PROTECT(RES = allocVector(INTSXP, *n));
  res = INTEGER(RES);

  for(i = 0; i<*n; i++){
    if(x[3*i + 1]>x[3*i])
      if(x[3*i + 2]>x[3*i + 1])
	res[i] = 2;
      else
	res[i] = 1;
    else if(x[3*i + 2]>x[3*i])
      res[i] = 2;
    else
      res[i] = 0;
  }

  UNPROTECT(1);
  return(RES);
}
