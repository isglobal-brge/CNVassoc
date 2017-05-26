#include "CNVassoc.h"


SEXP outerC(SEXP N, SEXP X, SEXP Y){
  int *n = INTEGER(N), i, j, k;
  double *x = REAL(X), *y = REAL(Y), *res;
  SEXP RES;
  
  PROTECT(RES = allocMatrix(REALSXP, 9, *n));
  res = REAL(RES);

  for(i = 0; i<*n; i++)
    for(j = 0; j<3; j++)
      for(k = 0; k<3; k++)
	res[9*i + 3*j + k] = x[3*i + k]*y[3*i + j];

  UNPROTECT(1);
  return(RES);
}
