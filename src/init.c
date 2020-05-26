#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
//extern SEXP btb_calculeQuantiles(SEXP, SEXP, SEXP);
extern SEXP btb_constituerMatriceEffectifs(SEXP, SEXP);
extern SEXP btb_constituerGrappes(SEXP, SEXP);
//extern SEXP btb_coordonneesGrappe(SEXP, SEXP);
extern SEXP btb_rcppLissage(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP btb_rcppLissageMedian(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP btb_rcppLissageMedianGrappe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
//  {"btb_calculeQuantiles",            (DL_FUNC) &btb_calculeQuantiles,            3},
  {"btb_constituerMatriceEffectifs",  (DL_FUNC) &btb_constituerMatriceEffectifs,  2},
  {"btb_constituerGrappes",           (DL_FUNC) &btb_constituerGrappes,           2},
//  {"btb_coordonneesGrappe",           (DL_FUNC) &btb_coordonneesGrappe,           2},
  {"btb_rcppLissage",                 (DL_FUNC) &btb_rcppLissage,                15},
  {"btb_rcppLissageMedian",           (DL_FUNC) &btb_rcppLissageMedian,           7},
  {"btb_rcppLissageMedianGrappe",     (DL_FUNC) &btb_rcppLissageMedianGrappe,    13},
  {NULL, NULL, 0}
};

void R_init_btb(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
