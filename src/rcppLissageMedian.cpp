#include <Rcpp.h>       // Rcpp::Rcout
#include <time.h>       // clock_t, clock, CLOCKS_PER_SEC
#include <math.h>       // pow 
#include <utility>      // std::swap

using namespace Rcpp;

/*
 *    date        version         auteur              commentaire
 * 2016/08/09      0.1.3      Arlindo Dos Santos      remplacement de l'appel a clock_gettime  par clock()
 *                                                    
 */

/*
 * ATTENTION: effet de bord: les 2 vecteurs en paramètre sont modifiés par la fonction
 * http://www.algolist.net/Algorithms/Sorting/Quicksort
 * 
 */
void quickSort(NumericVector vModalites, NumericVector vPonderations, int iLeft, int iRight) 
{
  int i = iLeft, j = iRight;
  double dPivot = vModalites[(iLeft + iRight) / 2];

  /* partition */
  while (i <= j)
  {
    while (vModalites[i] < dPivot)
      i++;
    
    while (vModalites[j] > dPivot)
      j--;
    
    if (i <= j)
    {
      std::swap(vModalites[i], vModalites[j]);
      std::swap(vPonderations[i], vPonderations[j]);
      i++;
      j--;
    }
  }
  
  /* recursion */
  if (iLeft < j)
    quickSort(vModalites, vPonderations, iLeft, j);
  
  if (i < iRight)
    quickSort(vModalites, vPonderations, i, iRight);
}

/*
 * arguments:
 *    - vModalites:    NumericVector contenant la liste des modalites
 *    - vPonderation:  NumericVector contenant la liste des pondérations
 *    - vQuantiles:    NumericVector contenant la liste des quantiles à déterminer
 * 
 * retourne
 *    - NumericVector contenant les quantiles
 */
// [[Rcpp::export]]
NumericVector calculeQuantiles(NumericVector vModalites, NumericVector vPonderation, NumericVector vQuantiles)
{
  int iNbModalites = vModalites.length();
  int iNbPonderations = vPonderation.length();
  int iNbQuantiles = vQuantiles.length();
  
  // debut des controles de validité
  if(iNbModalites != iNbPonderations)
  {
    Rcpp::Rcerr << "Il doit y avoir autant de modalités que de ponderations. NbModalites: " << iNbModalites << " - NbPonderations" << iNbPonderations << "\n";
    throw ""; 
  } 
  
  if(iNbModalites == 0)
  {
    Rcpp::Rcerr << "Il doit y avoir au moins une modalité";
    throw ""; 
  } 
  
  if(iNbQuantiles == 0)
  {
    Rcpp::Rcerr << "Il doit y avoir au moins un quantile";
    throw ""; 
  } 
  
  for(int i = 0; i < iNbQuantiles; ++i)
  {
    if(vQuantiles(i) < 0 || vQuantiles(i) > 1)
    {
      Rcpp::Rcerr << "Valeur de quantile invalide: " << vQuantiles(i);
      throw "";
    }
  }
  // fin des controles de validite
  
  NumericVector vPonderationCumulee(iNbPonderations);
  NumericVector vResultat(iNbQuantiles);
  long double dCible = 0;
  
  quickSort(vModalites, vPonderation, 0, iNbModalites - 1);

  vPonderationCumulee[0] = vPonderation[0];   // calcul des ponderations cumulees
  for(int iPonderation = 1; iPonderation < iNbPonderations; ++iPonderation)
    vPonderationCumulee[iPonderation] = vPonderationCumulee[iPonderation - 1] + vPonderation[iPonderation];

  for(int iQuantile = 0; iQuantile < vQuantiles.length(); ++iQuantile)    // pour tous les quantiles demandes
  {
    dCible = vQuantiles(iQuantile) * vPonderationCumulee(iNbPonderations - 1);
    
    for(int iPonderation = 0; iPonderation < iNbPonderations; ++iPonderation)
    {
      if(std::abs(vPonderationCumulee[iPonderation] - dCible) < std::numeric_limits<double>::epsilon())
      {
        if(iPonderation < iNbPonderations - 1)
        {
          vResultat[iQuantile] = (vModalites[iPonderation] + vModalites[iPonderation + 1]) / 2;
          break; 
        }
        else 
        {// dernier élément du tableau
          vResultat[iQuantile] = vModalites[iPonderation];
          break; 
        }
      }
      
      if(vPonderationCumulee[iPonderation] > dCible &&
         std::abs(vPonderationCumulee[iPonderation] - dCible) > std::numeric_limits<double>::epsilon())
      {
        vResultat[iQuantile] = vModalites[iPonderation];
        break;
      }
    }
  }
  
  return vResultat;
}

/*
 * 
 */
// [[Rcpp::export]]
NumericMatrix rcppLissageMedianSort(
    NumericVector vXobservations
  , NumericVector vYobservations
  , int iRayon
  , NumericMatrix mVar
  , NumericVector vXCentroides
  , NumericVector vYCentroides
  , NumericVector vQuantiles
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
  
  int i;
  int iVarCourante;
  int iIndiceObsCourante;
  int iNbVars = mVar.ncol();                // nombre de variables à traiter
  int iNbObs = vXobservations.length();     // nombre d'observations
  int iNbCentroides = vXCentroides.length();     // nombre de centroides
  int iNbQuantiles = vQuantiles.length();
  long double dRayonCarre = std::pow((long double)iRayon, 2); // rayon de lissage au carre
  long double dDistanceCarre;               // distance au carré entre une observation et un centroide
  
  long double dPonderation;
  List lResultat(iNbCentroides);
  NumericVector vValeursQuantiles(iNbVars * iNbQuantiles);
  NumericMatrix mValeursMedianes(iNbCentroides, iNbVars * iNbQuantiles + 1);
  
  /* on parcourt tous les centroides */
  for(int iIndiceCentroide = 0; iIndiceCentroide < iNbCentroides; ++iIndiceCentroide)
  {
    NumericVector vIndiceObservations; 
    NumericVector vPonderations; 
    dPonderation = 0;
    
    /* on parcourt toutes les observations */
    for(iIndiceObsCourante = 0; iIndiceObsCourante < iNbObs; ++iIndiceObsCourante)
    {
      dDistanceCarre = std::pow(long(vXobservations[iIndiceObsCourante]) - vXCentroides(iIndiceCentroide), 2) + std::pow(long(vYobservations[iIndiceObsCourante]) - vYCentroides(iIndiceCentroide), 2);
      
      if (dDistanceCarre < dRayonCarre)
      {
        dPonderation = std::pow(1 - (dDistanceCarre / dRayonCarre), 2);
        vPonderations.push_back(dPonderation);
        vIndiceObservations.push_back(iIndiceObsCourante);
      }
    }

    if(vIndiceObservations.length() > 0)
    {
      mValeursMedianes(iIndiceCentroide, 0) = vIndiceObservations.length();
      for (iVarCourante = 0; iVarCourante < iNbVars; ++iVarCourante) /* pour chacune des variables a lisser */
      {
        NumericVector vModalites;
        for(i = 0; i < vIndiceObservations.length(); ++i)
          vModalites.push_back(mVar(vIndiceObservations(i), iVarCourante));

        NumericVector vValeursQuantilesVarCourante = calculeQuantiles(vModalites, vPonderations, vQuantiles);
        for(i = 0; i < iNbQuantiles; ++i)
          mValeursMedianes(iIndiceCentroide, i + iVarCourante * iNbQuantiles + 1) = vValeursQuantilesVarCourante(i);
      }
    }
    // debut benchmark
    dTempsPasse = (clock() - timeBegin) / CLOCKS_PER_SEC;
    iPourcentageEffectue = 100 * iIndiceCentroide / iNbCentroides;
    if(iPourcentageEffectuePrecedent != iPourcentageEffectue)  
    {
      dTempsTotal = dTempsPasse * 100 / iPourcentageEffectue;
      iTempsRestant = ceil(dTempsTotal - dTempsPasse);
      iPourcentageEffectuePrecedent = iPourcentageEffectue;
      std::stringstream message;
      message << "\rMedian smoothing progress: " << iPourcentageEffectue << "% - minimum remaining time: " << (iTempsRestant / 60) << "m " << (iTempsRestant % 60) << "s";
      // if(updateProgress != NULL)
      //   updateProgress(iPourcentageEffectue, message.str());
      if(updateProgress.isNotNull())
        as<Function>(updateProgress)(iPourcentageEffectue, message.str());
      Rcpp::Rcout << message.str();
    }
    // fin benchmark
  }

  // debut benchmark
  Rcpp::Rcout << "\rElapsed time median smoothing: " << floor(dTempsTotal / 60) << "m " << ((int)dTempsTotal % 60)<< "s                                                                 ";
  // fin benchmark
  
  return(mValeursMedianes);
}
