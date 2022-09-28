#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
//extern SEXP btb_calculeQuantiles(SEXP, SEXP, SEXP);
extern SEXP _btb_constituerMatriceEffectifs(SEXP, SEXP);
extern SEXP _btb_constituerGrappes(SEXP, SEXP);
//extern SEXP btb_coordonneesGrappe(SEXP, SEXP);
extern SEXP _btb_rcppLissage(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _btb_rcppLissageMedian(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _btb_rcppLissageMedianGrappe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
//  {"btb_calculeQuantiles",            (DL_FUNC) &btb_calculeQuantiles,            3},
  {"_btb_constituerMatriceEffectifs",  (DL_FUNC) &_btb_constituerMatriceEffectifs,  2},
  {"_btb_constituerGrappes",           (DL_FUNC) &_btb_constituerGrappes,           2},
//  {"btb_coordonneesGrappe",           (DL_FUNC) &btb_coordonneesGrappe,           2},
  {"_btb_rcppLissage",                 (DL_FUNC) &_btb_rcppLissage,                15},
  {"_btb_rcppLissageMedian",           (DL_FUNC) &_btb_rcppLissageMedian,           7},
  {"_btb_rcppLissageMedianGrappe",     (DL_FUNC) &_btb_rcppLissageMedianGrappe,    13},
  {NULL, NULL, 0}
};

void R_init_btb(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
