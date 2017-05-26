#include <R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

SEXP maxC(SEXP N, SEXP X);
SEXP sampleC(SEXP N, SEXP X);
SEXP simcensrespC(SEXP SOBS, SEXP SCENS, SEXP TIMESOBS, SEXP TIMESCENS, SEXP TIMEMAX, SEXP TIMEMAXOBS, SEXP B, SEXP NSOBS, SEXP NSCENS, SEXP MINSOBS, SEXP MINSCENS);

SEXP logliklogistic(SEXP N, double *par, SEXP P, SEXP Y, SEXP WW);
SEXP logliklogisticcov(SEXP N, double *par, SEXP P, SEXP Y, SEXP WW, SEXP Q, SEXP C);
SEXP NRlogistic(SEXP N, SEXP PARAM, SEXP P, SEXP Y, SEXP WW, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);
SEXP NRlogisticcov(SEXP N, SEXP PARAM, SEXP P, SEXP Y, SEXP WW, SEXP Q, SEXP C, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);

SEXP loglikweibull(SEXP N, double *par, SEXP P, SEXP Y, SEXP CENS, SEXP WW);
SEXP loglikweibullcov(SEXP N, double *par, SEXP P, SEXP Y, SEXP CENS, SEXP WW, SEXP Q, SEXP C);
SEXP NRweibull(SEXP N, SEXP PARAM, SEXP P, SEXP Y, SEXP CENS, SEXP WW, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);
SEXP NRweibullcov(SEXP N, SEXP PARAM, SEXP P, SEXP Y, SEXP CENS, SEXP WW, SEXP Q, SEXP C, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);

SEXP loglikinterlogistic(SEXP N, double *par, SEXP Y, SEXP WW);
SEXP loglikinterlogisticcov(SEXP N, double *par, SEXP Y, SEXP WW, SEXP Q, SEXP C);
SEXP NRinterlogistic(SEXP N, SEXP PARAM, SEXP Y, SEXP WW, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);
SEXP NRinterlogisticcov(SEXP N, SEXP PARAM, SEXP Y, SEXP WW, SEXP Q, SEXP C, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);

SEXP loglikinterweibull(SEXP N, double *par, SEXP Y, SEXP CENS, SEXP WW);
SEXP loglikinterweibullcov(SEXP N, double *par, SEXP Y, SEXP CENS, SEXP WW, SEXP Q, SEXP C);
SEXP NRinterweibull(SEXP N, SEXP PARAM, SEXP Y, SEXP CENS, SEXP WW, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);
SEXP NRinterweibullcov(SEXP N, SEXP PARAM, SEXP Y, SEXP CENS, SEXP WW, SEXP Q, SEXP C, SEXP TOL, SEXP MAX_ITER, SEXP VERBOSE);
     