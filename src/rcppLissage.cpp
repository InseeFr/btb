// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>       
#include <RcppParallel.h>
#include <vector>
#include <armadillo>
#include <time.h>       // clock_t, clock, CLOCKS_PER_SEC
#include <math.h>       // pow 

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;
using namespace arma;

/*
*    date        version         auteur              commentaire
* 2016/08/09      0.1.3      Arlindo Dos Santos      remplacement de l'appel à clock_gettime  par clock()
*                                                  
*/

/* 
*  arguments
*  vXobservations : vecteur avec les coordonnées x des observations
*  vYobservations : vecteur avec les coordonnées y des observations
*  vColonneObservation      : vecteur avec l'indice colonne des observations 
*  vLigneObservation      : vecteur avec l'indice ligne des observations 
*  iPas           : longueur du côté d'un carreau
*  iRayon         : rayon de lissage         
*  iNeighbor      : voisinage étendu                                 
*  mVar           : matrice avec les variables non lisees
*  mXcentroides   : matrice contenant la valeur de la coordonnée x pour les centroides, stocke en (i;j) ou i est le numero d'ordre en abscisse du carreau et j le numero d'ordre en ordonnée du carreau
*  mYcentroides   : matrice contenant la valeur de la coordonnée y pour les centroides, stocke en (i;j) ou i est le numero d'ordre en abscisse du carreau et j le numero d'ordre en ordonnée du carreau
*  mIcentroides   : matrice contenant le numéro d'ordre du centroide, stocké en (i;j) où i est le numéro d'ordre en abscisse du carreau et j le numéro d'ordre en ordonnée du carreau; la valeur -1 doit figurer si le carreau n'est pas dans la grille
*  iNbCentroides  : nombre de centroides
*  
*  retourne 
*  une matrice "mVariablesLissees" contenant les différentes variables en colonnes et les valeurs lissées en ligne
*  cette matrice ne comporte pas les coordonnées des centroides ni des observations
*/
// [[Rcpp::export]]
NumericMatrix rcppLissage(
    IntegerVector vXObservation
  , IntegerVector vYObservation
  , IntegerVector vLigneObservation
  , IntegerVector vColonneObservation
  , int iPas
  , int iRayon
  , int iNeighbor
  , NumericMatrix mVariables

  , int iNumberCols
  , int iNumberRows
  , int iMinXCentroide
  , int iMinYCentroide
  
  , IntegerMatrix mIcentroide
  , int iNbCentroides
  , Nullable <Function> updateProgress = R_NilValue
)
{
  // debut cast pour threadsafety - cf https://rcppcore.github.io/RcppParallel/#thread_safety
  const RVector<int> vXObservations(vXObservation);
  const RVector<int> vYObservations(vYObservation);
  const RVector<int> vLigneObs(vLigneObservation);
  const RVector<int> vColonneObs(vColonneObservation);
  const RMatrix<double> mVar(mVariables);
  
  const RMatrix<int> mIcentroides(mIcentroide);
  // fin cast pour threadsafety 
  
  // debut benchmark
  clock_t timeBegin = clock();
  double dTempsPasse;
  double dTempsTotal = 0;
  int iTempsRestant = 0;
  int iPourcentageEffectue;
  int iPourcentageEffectuePrecedent = 0;
  std::stringstream message;
  // fin benchmark
  
  const int iNbVoisins = ceil(double(iRayon) / iPas - 0.5);   // fenêtrage: le cercle de rayon iRayon est circonscrit dans le rectangle de 2*iNbVoisins * 2*iNbVoisins
  
  const int iNbCentroidesAbscisse = iNumberCols;
  const int iNbCentroidesOrdonnee = iNumberRows;
  
  const int iNbVars = mVar.ncol();                // nombre de variables a traiter
  const int iNbObs = vXObservations.length();     // nombre d'observations

  int iCol, iLigne;
  int iColMin, iColMax, iLigneMin, iLigneMax;
  int iVarCourante;
  long double dRayonCarre = std::pow((long double)iRayon, 2); // rayon de lissage au carre
  long double dDistanceCarre;               // distance au carré entre une observation et un centroide
  long double dSommePonderation;            // double contenant la somme des ponderations qui sont appliquées depuis l'observation considérée
  
  arma::mat mPonderation(iNbCentroidesOrdonnee, iNbCentroidesAbscisse);   // matrice contenant la ponderation à appliquer à la valeur du centroide pour l'observation considérée
  arma::mat mVariablesLissees(iNbCentroides, iNbVars, fill::zeros);  // matrice contenant le resultat (variables lissees)
  
  /* on parcourt toutes les observations */
  for(int iIndiceObsCourante = 0; iIndiceObsCourante < iNbObs; ++iIndiceObsCourante)
  {
    dSommePonderation = 0;
    
    // fenêtrage
    iColMin = fmax(vColonneObs[iIndiceObsCourante] - iNbVoisins - 1, 0);
    iColMax = fmin(vColonneObs[iIndiceObsCourante] + iNbVoisins, iNbCentroidesAbscisse - 1);
    iLigneMin = fmax(vLigneObs[iIndiceObsCourante] - iNbVoisins - 1, 0);
    iLigneMax = fmin(vLigneObs[iIndiceObsCourante] + iNbVoisins, iNbCentroidesOrdonnee - 1);
    
    /* on parcourt tous les centroides susceptibles d'être intéressé par cette observation */
    for(iCol = iColMin; iCol <= iColMax; ++iCol)
    {
      for(iLigne = iLigneMin; iLigne <= iLigneMax; ++iLigne)
      {
        // vérifier que le centroide est bien dans la grille fournie par l'utilisateur
        if(mIcentroides(iLigne, iCol) == -1)
          continue;
        
        dDistanceCarre = std::pow(long(vXObservations[iIndiceObsCourante]) - (iMinXCentroide + iCol * iPas), 2) + std::pow(long(vYObservations[iIndiceObsCourante]) - (iMinYCentroide + iLigne * iPas), 2);

        if (dDistanceCarre < dRayonCarre)
        {
          mPonderation(iLigne, iCol) = std::pow(1 - (dDistanceCarre / dRayonCarre), 2);
          dSommePonderation += mPonderation(iLigne, iCol);
        }
        else
          mPonderation(iLigne, iCol) = 0; // à conserver impérativement pour ne pas avoir à refaire le calcul de distance lors de la normalisation juste ci-dessous
      }
    }
    
    if(dSommePonderation > 0)
    {
      for(iCol = iColMin; iCol <= iColMax; ++iCol)
      {
        for(iLigne = iLigneMin; iLigne <= iLigneMax; ++iLigne)
        {
          if(mIcentroides(iLigne, iCol) == -1)
            continue;
          
          for (iVarCourante = 0; iVarCourante < iNbVars; ++iVarCourante) /* pour chacune des variables a lisser */
          {
            /* on calcule le lissage : pondération de la valeur de la variable par le poids de lissage afin de normaliser sa valeur */
            mVariablesLissees(mIcentroides(iLigne, iCol), iVarCourante) += mVar(iIndiceObsCourante, iVarCourante) * mPonderation(iLigne, iCol) / dSommePonderation;
          }
          mPonderation(iLigne, iCol) = 0;
        }
      }
    }
    
    // debut benchmark
    dTempsPasse = (clock() - timeBegin) / CLOCKS_PER_SEC;
    iPourcentageEffectue = 100 * iIndiceObsCourante / iNbObs;
    if(iPourcentageEffectuePrecedent != iPourcentageEffectue)  
    {
      dTempsTotal = dTempsPasse * 100 / iPourcentageEffectue;
      iTempsRestant = ceil(dTempsTotal - dTempsPasse);
      iPourcentageEffectuePrecedent = iPourcentageEffectue;
      message.str("");
      message << "Smoothing progress: " << iPourcentageEffectue << "% - minimum remaining time: " << (iTempsRestant / 60) << "m " << (iTempsRestant % 60) << "s";
      if(updateProgress.isNotNull())
        as<Function>(updateProgress)(iPourcentageEffectue, message.str());
      else
        Rcpp::Rcout << "\r" << message.str();
    }
    // fin benchmark
  }
  
  // debut benchmark
  message.str("");
  message << "Smoothing duration: " << floor(dTempsTotal / 60) << "m " << ((int)dTempsTotal % 60) << "s                                                                                           ";
  if(updateProgress.isNotNull())
    as<Function>(updateProgress)(iPourcentageEffectue, message.str());
  else
    Rcpp::Rcout << "\n" << message.str();
  // fin benchmark
  
  return(wrap(mVariablesLissees));
}
