#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP NRinterlogistic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NRinterlogisticcov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NRinterweibull(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NRinterweibullcov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NRlogistic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NRlogisticcov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NRweibull(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NRweibullcov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP outerC(SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(dqrls)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"NRinterlogistic",    (DL_FUNC) &NRinterlogistic,     7},
    {"NRinterlogisticcov", (DL_FUNC) &NRinterlogisticcov,  9},
    {"NRinterweibull",     (DL_FUNC) &NRinterweibull,      8},
    {"NRinterweibullcov",  (DL_FUNC) &NRinterweibullcov,  10},
    {"NRlogistic",         (DL_FUNC) &NRlogistic,          8},
    {"NRlogisticcov",      (DL_FUNC) &NRlogisticcov,      10},
    {"NRweibull",          (DL_FUNC) &NRweibull,           9},
    {"NRweibullcov",       (DL_FUNC) &NRweibullcov,       11},
    {"outerC",             (DL_FUNC) &outerC,              3},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"dqrls", (DL_FUNC) &F77_NAME(dqrls), 13},
    {NULL, NULL, 0}
};

void R_init_CNVassoc(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
