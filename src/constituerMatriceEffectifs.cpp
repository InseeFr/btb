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

//' @title constituerMatriceEffectifs
//' @name constituerMatriceEffectifs
//' @description 
//' Function constituting a membership matrix (Fonction constituant une matrice des effectifs)
//' @usage constituerMatriceEffectifs(vLigneObservation, vColonneObservation)
//'
//' @param 
//'  - vLigneObservation : 
//'      A \code{vector} containing the line number of each observation
//'      (Un \code{vector} contenant le numéro de ligne de chaque observation.)
//'  - vColonneObservation :
//'      A \code{vector} containing the column number of each observation
//'      ()Un \code{vector} contenant le numéro de colonne de chaque observation.)
//' @return 
//'  Returns a \code{matrix} with the number of observations for each cell.
//'  
//'  (Retourne une \code{matrix} avec le nombre d'observations pour chaque cellule.)
//'
//' @author
//'    - Psar Analyse Urbaine Insee 
//'    - Arlindo Dos Santos 
//'    - Francois Semecurbe}
//'
//' @examples : 
//'  dfObservations <- data.frame(x = c(15, 35, 15, 25, 35, 55, 45, 45, 55, 65, 70, 75, 85, 90,
//'                                     65, 75, 85, 65, 70, 75, 85, 90, 65, 70, 75)
//'                                 ,  y = c(10, 10, 30, 30, 35, 35, 45, 55, 55, 65, 65, 65, 65, 65,
//'                                          70, 70, 70, 75, 75, 75, 75, 75, 85, 85, 85))
//'  cellSize <- 20L
//'  # calcul de l'indice des observations
//'  # on prend le rectangle englobant 
//'  # et on positionne le debut de la numérotation sur la première observation
//'  dfObservations$col <- as.integer(floor((dfObservations$x) / cellSize)
//'                                     - floor(min(dfObservations$x / cellSize)) + 1)
//'  dfObservations$row <- as.integer(floor((dfObservations$y) / cellSize) 
//'                                     - floor(min(dfObservations$y / cellSize)) + 1)
//'  
//'  mEffectifs <- constituerMatriceEffectifs(dfObservations$row - 1, dfObservations$col - 1)


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

