// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#ifndef RCPP_LISSAGE_H
#define RCPP_LISSAGE_H

NumericMatrix rcppLissage(
    IntegerVector vXObservation
  , IntegerVector vYObservation
  , IntegerVector vLigneObservation
  , IntegerVector vColonneObservation
  , int iPas
  , int iRayon
  , int iNeighbor
  , NumericMatrix mVariables
  , IntegerMatrix mXcentroide
  , IntegerMatrix mYcentroide
  , IntegerMatrix mIcentroide
  , int iNbCentroides
  , Nullable <Function> updateProgress = R_NilValue
)
