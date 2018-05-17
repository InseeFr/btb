// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#ifndef LISSAGE_MEDIAN_H
#define LISSAGE_MEDIAN_H

std::vector<double> calculeQuantiles(std::vector<double>& vModalites, std::vector<double>& vPonderation, const std::vector<double> vQuantiles);

arma::mat rcppLissageMedian(
    std::vector<int> vXobservations
  , std::vector<int> vYobservations
  , int iRayon
  , arma::mat& mVar
  , std::vector<int> vXCentroides
  , std::vector<int> vYCentroides
  , std::vector<double> vQuantiles
);

#endif
