#include "CNVassoc.h"

// function to simulate right-censored data from empirical distribution.
SEXP simcensrespC(SEXP SOBS, SEXP SCENS, SEXP TIMESOBS, SEXP TIMESCENS, SEXP TIMEMAX, SEXP TIMEMAXOBS, SEXP B, SEXP NSOBS, SEXP NSCENS, SEXP MINSOBS, SEXP MINSCENS){

  int *b = INTEGER(B), *nsobs = INTEGER(NSOBS), *nscens = INTEGER(NSCENS), i, iter;
  double *sobs = REAL(SOBS), *scens = REAL(SCENS), *timesobs = REAL(TIMESOBS), *timescens = REAL(TIMESCENS), *timemax = REAL(TIMEMAX), *timemaxobs = REAL(TIMEMAXOBS), *minsobs = REAL(MINSOBS), *minscens = REAL(MINSCENS), r, tt = .0, cc = .0, *res;
  SEXP RES;
  
  PROTECT(RES = allocMatrix(REALSXP, 2, *b));
  res = REAL(RES);

  GetRNGstate();
  for (iter = 0; iter<*b; iter++){	

    // observed 
    r = runif(0,1);
    if (r < *minsobs)
      tt = 10000000000;
    else{
      for (i = 1; i<*nsobs; i++){
        if (r <= sobs[i-1] && r >= sobs[i]){
          if (r >= (sobs[i-1] + sobs[i])/2)
            if (timesobs[i-1] >= timesobs[i])
              tt = timesobs[i];
            else
              tt = runif(timesobs[i-1],timesobs[i]);
          else
            if (timesobs[i] >= timesobs[i+1])
              tt = timesobs[i];
            else
              tt = runif(timesobs[i],timesobs[i+1]);            
          break;
        }
      }
    }
    // censored     
    r = runif(0,1);
    if (r < *minscens)
      cc = runif(*timemaxobs,*timemax);
    else{
      for (i = 1; i<*nscens; i++){
        if (r <= scens[i-1] && r >= scens[i]){
          if (r >= (scens[i-1]+scens[i])/2){
            if (timescens[i-1] >= timescens[i]) 
              cc = timescens[i];
            else
              cc = runif(timescens[i-1],timescens[i]);
          }else{
            if (timescens[i] >= timescens[i+1])
              cc = timescens[i];
            else 
              cc = runif(timescens[i],timescens[i+1]);            
          }
          break;
        }
      }
    }
    if (tt<cc){
      res[0+2*iter] = tt;
      res[1+2*iter] = 1.0;
    } else {
      res[0+2*iter] = cc;
      res[1+2*iter] = 0.0;
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return(RES);
  
}
