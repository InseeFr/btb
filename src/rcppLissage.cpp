#include <Rcpp.h>
#include <time.h>
#include <iostream>
using namespace Rcpp;

/*
 *    date        version         auteur              commentaire
 * 2016/08/09      0.1.3      Arlindo Dos Santos      remplacement de l'appel a clock_gettime  par clock()
 *                                                  
 */

/* 
*  arguments
*  vXobservations : vecteur avec les coordonnees x des observations
*  vYobservations : vecteur avec les coordonnees y des observations
*  vIndicesX      : vecteur avec l'indice i des observations 
*  vIndicesY      : vecteur avec l'indice j des observations 
*  iPas           : longueur du cote d'un carreau
*  iRayon         : rayon de lissage                                          
*  mVar           : matrice avec les variables non lisees
*  mXcentroides   : matrice contenant la valeur de la coordonnee x pour les centroides, stocke en (i;j) ou i est le numero d'ordre en abscisse du carreau et j le numero d'ordre en ordonnee du carreau
*  mYcentroides   : matrice contenant la valeur de la coordonnee y pour les centroides, stocke en (i;j) ou i est le numero d'ordre en abscisse du carreau et j le numero d'ordre en ordonnee du carreau
*  mIcentroides   : matrice contenant le numero d'ordre du centroide, stocke en (i;j) ou i est le numero d'ordre en abscisse du carreau et j le numero d'ordre en ordonnee du carreau
*  iNbCentroides  : nombre de centroides
*  
*  retourne 
*  une matrice "mVariablesLissees" contenant les differentes variables en colonnes et les valeurs lissees en ligne
*  cette matrice ne comporte pas les coordonnees des centroides ni des observations
*/
// [[Rcpp::export]]
NumericMatrix rcppLissage(
                   NumericVector vXobservations
                  ,NumericVector vYobservations
                  ,NumericVector vIndicesX
                  ,NumericVector vIndicesY
                  ,int iPas
                  ,int iRayon
                  ,NumericMatrix mVar
                  ,NumericMatrix mXcentroides
                  ,NumericMatrix mYcentroides
                  ,NumericMatrix mIcentroides
                  ,int iNbCentroides
)
{
  // debut benchmark
  clock_t timeBegin = clock();
  double dTempsPasse;
  double dTempsTotal = 0;
  int iTempsRestant = 0;
  int iPourcentageEffectue;
  int iPourcentageEffectuePrecedent = 0;
  // fin benchmark

  int i, j;
  int iMin, iMax, jMin, jMax;
  int iVarCourante;
  int iNbVoisins = ceil(double(iRayon) / iPas);   // nombre de cases voisines a etudier pour rechercher des centroides a alimenter
  int iNbCentroidesAbscisse = mXcentroides.ncol();
  int iNbCentroidesOrdonnee = mXcentroides.nrow();
  int iNbVars = mVar.ncol();                // nombre de variables a traiter
  int iNbObs = vXobservations.length();     // nombre d'observations
  long double dRayonCarre = pow((long double)iRayon, 2); // rayon de lissage au carre
  long double dDistanceCarre;               // distance au carre entre une observation et un centroide
  long double dSommePonderation;            // double contenant la somme des ponderations qui sont appliquees depuis l'observation consideree

  NumericMatrix mPonderation(iNbCentroidesOrdonnee, iNbCentroidesAbscisse);   // matrice contenant la ponderation a appliquer a la valeur du centroide pour l'observation consideree
  NumericMatrix mVariablesLissees(iNbCentroides, iNbVars);  // matrice contenant le resultat (variables lissees)

  /* on parcourt toutes les observations */
  for(int iIndiceObsCourante = 0; iIndiceObsCourante < iNbObs; ++iIndiceObsCourante)
  {
    dSommePonderation = 0;

    // i est l'indice en abscisse (numero de colonne) 
    iMin = fmax(vIndicesX[iIndiceObsCourante] - iNbVoisins - 1, 0);
    iMax = fmin(vIndicesX[iIndiceObsCourante] + iNbVoisins, iNbCentroidesAbscisse - 1);
    
    // j est l'indice en ordonnee (numero de ligne)
    jMin = fmax(vIndicesY[iIndiceObsCourante] - iNbVoisins - 1, 0);
    jMax = fmin(vIndicesY[iIndiceObsCourante] + iNbVoisins, iNbCentroidesOrdonnee - 1);

    /* on parcourt tous les centroides suceptibles d'etre interesse par cette observation */
    for(i = iMin; i <= iMax; ++i)
    {
      for(j = jMin; j <= jMax; ++j)
      {
        dDistanceCarre = pow(long(vXobservations[iIndiceObsCourante]) - mXcentroides(j, i), 2) + pow(long(vYobservations[iIndiceObsCourante]) - mYcentroides(j, i), 2);

        if (dDistanceCarre < dRayonCarre)
        {
          mPonderation(j, i) = pow(1 - (dDistanceCarre / dRayonCarre), 2);
          dSommePonderation += mPonderation(j, i);
        }
      }
    }

    if(dSommePonderation > 0)
    {
      for(i = iMin; i <= iMax; ++i)
      {
        for(j = jMin; j <= jMax; ++j)
        {
          for (iVarCourante = 0; iVarCourante < iNbVars; ++iVarCourante) /* pour chacune des variables a lisser */
          {
            /* on calcule le lissage : ponderation de la valeur de la variable par le poids de lissage afin de normaliser sa valeur */
            mVariablesLissees(mIcentroides(j, i), iVarCourante) += mVar(iIndiceObsCourante, iVarCourante) * mPonderation(j, i) / dSommePonderation;
          }
          mPonderation(j, i) = 0;
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
      Rcpp::Rcout << "\rSmoothing progress: " << iPourcentageEffectue << "% - remaining time: " << (iTempsRestant / 60) << "m " << (iTempsRestant % 60) << "s                     ";
    }
    // fin benchmark
  }

  // debut benchmark
  Rcpp::Rcout << "\rElapsed time smoothing: " << floor(dTempsTotal / 60) << "m " << ((int)dTempsTotal % 60) << "s                                                                                           ";
  // fin benchmark
  
  return(mVariablesLissees);
}
