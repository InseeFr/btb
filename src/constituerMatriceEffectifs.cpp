// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>       
#include <RcppParallel.h>
#include <armadillo>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
arma::Mat<int> constituerMatriceEffectifs(IntegerVector vLigneObservation, IntegerVector vColonneObservation)
{
  RVector<int> vLigneObservations(vLigneObservation);
  RVector<int> vColonneObservations(vColonneObservation);
  
  const unsigned int iMaxLigne = *max_element(vLigneObservations.begin(), vLigneObservations.end()) + 1;
  const unsigned int iMaxColonne = *max_element(vColonneObservations.begin(), vColonneObservations.end()) + 1;
  const unsigned int iMaxDimension = std::max(iMaxLigne, iMaxColonne);
  const unsigned int iTailleMatrice = pow(2, std::ceil( log(iMaxDimension) / log(2.0)));
  const unsigned int iNbObs = vLigneObservations.length(); // nombre d'observations total
  
  arma::Mat<int> mEffectifs2n(iTailleMatrice, iTailleMatrice, fill::zeros); // la matrice des effectifs pour chacun des carreaux
  
  for(unsigned int iObsCourante = 0; iObsCourante < iNbObs; ++iObsCourante)
    mEffectifs2n(vLigneObservations[iObsCourante], vColonneObservations[iObsCourante]) += 1;
  
  return mEffectifs2n;
}

