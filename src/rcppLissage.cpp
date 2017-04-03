#include <Rcpp.h>       // Rcpp::Rcout
#include <time.h>       // clock_t, clock, CLOCKS_PER_SEC
#include <math.h>       // pow 

using namespace Rcpp;

/*
 *    date        version         auteur              commentaire
 * 2016/08/09      0.1.3      Arlindo Dos Santos      remplacement de l'appel à clock_gettime  par clock()
 *                                                  
 */

/* 
*  arguments
*  vXobservations : vecteur avec les coordonnées x des observations
*  vYobservations : vecteur avec les coordonnées y des observations
*  vIndicesX      : vecteur avec l'indice i des observations 
*  vIndicesY      : vecteur avec l'indice j des observations 
*  iPas           : longueur du côté d'un carreau
*  iRayon         : rayon de lissage         
*  iNeighbor      : voisinage étendu                                 
*  mVar           : matrice avec les variables non lisees
*  mXcentroides   : matrice contenant la valeur de la coordonnée x pour les centroides, stocke en (i;j) ou i est le numero d'ordre en abscisse du carreau et j le numero d'ordre en ordonnée du carreau
*  mYcentroides   : matrice contenant la valeur de la coordonnée y pour les centroides, stocke en (i;j) ou i est le numero d'ordre en abscisse du carreau et j le numero d'ordre en ordonnée du carreau
*  mIcentroides   : matrice contenant le numéro d'ordre du centroide, stocké en (i;j) où i est le numéro d'ordre en abscisse du carreau et j le numéro d'ordre en ordonnée du carreau
*  iNbCentroides  : nombre de centroides
*  
*  retourne 
*  une matrice "mVariablesLissees" contenant les différentes variables en colonnes et les valeurs lissées en ligne
*  cette matrice ne comporte pas les coordonnées des centroides ni des observations
*/
// [[Rcpp::export]]
NumericMatrix rcppLissage(
                    NumericVector vXobservations
                  , NumericVector vYobservations
                  , NumericVector vIndicesX
                  , NumericVector vIndicesY
                  , int iPas
                  , int iRayon
                  , int iNeighbor
                  , NumericMatrix mVar
                  , NumericMatrix mXcentroides
                  , NumericMatrix mYcentroides
                  , NumericMatrix mIcentroides
                  , int iNbCentroides
                  , Nullable <Function> updateProgress = R_NilValue
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
  int iNbVoisins = ceil(double(iRayon) / iPas - 0.5);   // le cercle de rayon iRayon est circonscrit dans le rectangle de 2*iNbVoisins * 2*iNbVoisins
  int iNbCentroidesAbscisse = mXcentroides.ncol();
  int iNbCentroidesOrdonnee = mXcentroides.nrow();
  int iNbVars = mVar.ncol();                // nombre de variables a traiter
  int iNbObs = vXobservations.length();     // nombre d'observations
  long double dRayonCarre = std::pow((long double)iRayon, 2); // rayon de lissage au carre
  long double dDistanceCarre;               // distance au carré entre une observation et un centroide
  long double dSommePonderation;            // double contenant la somme des ponderations qui sont appliquées depuis l'observation considérée

  NumericMatrix mPonderation(iNbCentroidesOrdonnee, iNbCentroidesAbscisse);   // matrice contenant la ponderation a appliquer a la valeur du centroide pour l'observation consideree
  NumericMatrix mVariablesLissees(iNbCentroides, iNbVars);  // matrice contenant le resultat (variables lissees)
  // NumericMatrix mVariablesLissees(iNbCentroides, iNbVars + 1);  // matrice contenant le resultat (variables lissees) version avec nbObsPondere
  
  /* on parcourt toutes les observations */
  for(int iIndiceObsCourante = 0; iIndiceObsCourante < iNbObs; ++iIndiceObsCourante)
  {
    dSommePonderation = 0;
    
    // i est l'indice en abscisse (numéro de colonne) 
    iMin = fmax(vIndicesX[iIndiceObsCourante] + iNeighbor - iNbVoisins - 1, 0);
    iMax = fmin(vIndicesX[iIndiceObsCourante] + iNeighbor + iNbVoisins, iNbCentroidesAbscisse - 1);
    
    // j est l'indice en ordonnée (numéro de ligne)
    jMin = fmax(vIndicesY[iIndiceObsCourante] + iNeighbor - iNbVoisins - 1, 0);
    jMax = fmin(vIndicesY[iIndiceObsCourante] + iNeighbor + iNbVoisins, iNbCentroidesOrdonnee - 1);
    
    /* on parcourt tous les centroides susceptibles d'être intéréssé par cette observation */
    for(i = iMin; i <= iMax; ++i)
    {
      for(j = jMin; j <= jMax; ++j)
      {
        // vérifier que le centroide est bien dans la grille fournie par l'utilisateur
        if(NumericVector::is_na (mXcentroides(j, i)))
          continue;

        dDistanceCarre = std::pow(long(vXobservations[iIndiceObsCourante]) - mXcentroides(j, i), 2) + std::pow(long(vYobservations[iIndiceObsCourante]) - mYcentroides(j, i), 2);

        if (dDistanceCarre < dRayonCarre)
        {
          mPonderation(j, i) = std::pow(1 - (dDistanceCarre / dRayonCarre), 2);
          dSommePonderation += mPonderation(j, i);
          // mVariablesLissees(mIcentroides(j, i), iNbVars) += mPonderation(j, i); // calcul de la colonne nbObsPondere
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
            /* on calcule le lissage : pondération de la valeur de la variable par le poids de lissage afin de normaliser sa valeur */
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
      std::stringstream message;
      message << "\rSmoothing progress: " << iPourcentageEffectue << "% - minimum remaining time: " << (iTempsRestant / 60) << "m " << (iTempsRestant % 60) << "s";
      // if(updateProgress != NULL)
      //   updateProgress(iPourcentageEffectue, message.str());
      if(updateProgress.isNotNull())
        as<Function>(updateProgress)(iPourcentageEffectue, message.str());
      Rcpp::Rcout << message.str();
    }
    // fin benchmark
  }

  // debut benchmark
  Rcpp::Rcout << "\rElapsed time smoothing: " << floor(dTempsTotal / 60) << "m " << ((int)dTempsTotal % 60) << "s                                                                                           ";
  // fin benchmark
  
  return(mVariablesLissees);
}
