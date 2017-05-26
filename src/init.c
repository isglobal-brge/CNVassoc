#include "CNVassoc.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

R_CallMethodDef CallMethods[]  = {
  {"maxC", (DL_FUNC) &maxC, 3}, //maxC}
  {"sampleC", (DL_FUNC) &sampleC, 2}, //ff_NewtonRaphson
  {"NRlogistic", (DL_FUNC) &NRlogistic, 6}, //ff_NewtonRaphson
  {"NRweibull", (DL_FUNC) &NRweibull, 7}, //ff_NewtonRaphson
  {"NRlogisticcov", (DL_FUNC) &NRlogisticcov, 8}, //ff_NewtonRaphson
  {"NRweibullcov", (DL_FUNC) &NRweibullcov, 9}, //ff_NewtonRaphson
  {NULL, NULL, 0}
};


void R_init_test(DllInfo *ddl){
  R_registerRoutines(ddl, NULL, CallMethods, NULL, NULL);
}
