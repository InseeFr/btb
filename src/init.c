#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP btb_calculeQuantiles(SEXP, SEXP, SEXP);
extern SEXP btb_rcppLissage(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP btb_rcppLissageMedianSort(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"btb_calculeQuantiles",      (DL_FUNC) &btb_calculeQuantiles,       3},
  {"btb_rcppLissage",           (DL_FUNC) &btb_rcppLissage,           13},
  {"btb_rcppLissageMedianSort", (DL_FUNC) &btb_rcppLissageMedianSort,  8},
  {NULL, NULL, 0}
};

void R_init_btb(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
