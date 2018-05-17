// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#ifndef RCPP_LISSAGE_MEDIAN_GRAPPE_H
#define RCPP_LISSAGE_MEDIAN_GRAPPE_H

NumericMatrix rcppLissageMedianGrappe(
    int iMinObsGrappe
  , IntegerVector vXObservation
  , IntegerVector vYObservation
  , IntegerVector vLigneObservation
  , IntegerVector vColonneObservation
  , int iPas
  , int iRayon
  , NumericMatrix mVariables
  , IntegerVector vXCentroide
  , IntegerVector vYCentroide
  , IntegerVector vLigneCentroide
  , IntegerVector vColonneCentroide
  , NumericVector vQuantile
);

#endif
