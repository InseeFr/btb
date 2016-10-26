#include <Rcpp.h>
#include <time.h>
#include <iostream>
using namespace Rcpp;

// double dAccumulateurTempsTri = 0;
  
// /*
//  * Fonction qui retourne l'indice de l'element le plus petit dans un intervalle sur un vecteur de pairs
//  * La comparaison se fait sur sur le premier element de la paire 
//  * 
//  * arguments
//  *    vVecteurPair: vecteur de pairs ou effectuer la recherche
//  *    iDebut: index du premier element du vecteur a prendre en consideration
//  *    iFin: index du dernier element du vecteur a prendre en consideration
//  * 
//  * resultat
//  *    indice de l'element ayant la plus petite valeur
//  */
// int getIndiceElementMin(std::vector <std::pair <double, double> > vVecteurPair, int iDebut, int iFin)
// {
//   int iIndiceResultat = iDebut;
//   double dMinVal = vVecteurPair[iDebut].first;
//   for(int i = iDebut + 1; i <= iFin; ++i)
//   {
//     if(vVecteurPair[i].first < dMinVal)
//     {
//       dMinVal = vVecteurPair[i].first;
//       iIndiceResultat = i;
//     }
//   }
//   return iIndiceResultat;
// }
// 
// /*
//  * Recherche de la mediane pour des frequences decimales 
//  * 
//  * Application de la recherche quickselect pour une complexite moyenne en O(n)
//  *    https://fr.wikipedia.org/wiki/Quickselect
//  *    http://stackoverflow.com/questions/34271792/rcpp-function-to-find-the-median-given-a-vector-of-values-and-their-frequencies
//  * 
//  * arguments
//  *    vValues: vecteur de valeurs
//  *    vFreqs: vecteur de frequences
//  * 
//  * resultat
//  *    double repr?sentant la valeur mediane
//  * 
//  */
// double fast_median_freq(std::vector<double> vValues, std::vector<double> vFreqs) 
// {
//   const int iLength = vFreqs.size();
//   std::vector < std::pair <double, double> > vDonnees;
//   double dFreqSum = 0;
//   
//   for (int i = 0; i < iLength; ++i) 
//   {
//     vDonnees.push_back(std::pair<double, double>(vValues[i], vFreqs[i]));
//     dFreqSum += vFreqs[i];
//   }
//   
//   double dTarget = dFreqSum / 2;
//   int iLow = 0;
//   int iHigh = iLength - 1;
//   
//   while (true) 
//   {
//     // Random pivot; move to the end
//     int iRandomIndex = iLow + (rand() % (iHigh - iLow + 1));
//     std::swap(vDonnees[iRandomIndex], vDonnees[iHigh]);
//     
//     // In-place pivot
//     int iHighPos = iLow;  // Start of values higher than pivot
//     double dLowSum = 0;  // Sum of frequencies of elements below pivot
// 
//     for (int iPos = iLow; iPos < iHigh; ++iPos) 
//     {
//       if (vDonnees[iPos].first <= vDonnees[iHigh].first) 
//       {
//         dLowSum += vDonnees[iPos].second;
//         std::swap(vDonnees[iHighPos], vDonnees[iPos]);
//         ++iHighPos;
//       }
//     }
//     std::swap(vDonnees[iHighPos], vDonnees[iHigh]);  // Move pivot to "iHighPos"
// 
//     // for(int i = 0; i < iLength; ++i)
//     //   Rcpp::Rcout << "\n vDonnees[" << i << "]: " << vDonnees[i].first << " => " << vDonnees[i].second;
//     // 
//     // Rcpp::Rcout << "\n std::numeric_limits<double>::epsilon(): " << std::numeric_limits<double>::epsilon();
//     // 
//     // Rcpp::Rcout << "\n iHighPos: " << iHighPos << " - iLow: " << iLow << " - iHigh: " << iHigh;
//     // Rcpp::Rcout << "\n dLowSum: " << dLowSum << " - vDonnees[iHighPos].second: " << vDonnees[iHighPos].second << " - dTarget: " << dTarget;
//     // Rcpp::Rcout << "\n (dLowSum + vDonnees[iHighPos].second): " << (dLowSum + vDonnees[iHighPos].second);
//     // Rcpp::Rcout << "\n (dLowSum == dTarget): " << (dLowSum == dTarget);
//     // Rcpp::Rcout << "\n (dLowSum < dTarget): " << (dLowSum < dTarget);
//     // Rcpp::Rcout << "\n (dLowSum + vDonnees[iHighPos].second) gt dTarget: " << ((dLowSum + vDonnees[iHighPos].second) > dTarget);
//     // Rcpp::Rcout << "\n (dLowSum + vDonnees[iHighPos].second) - dTarget: " << (dLowSum + vDonnees[iHighPos].second) - dTarget;
//     // Rcpp::Rcout << "\n std::fabs(dLowSum + vDonnees[iHighPos].second - dTarget): " << std::fabs(dLowSum + vDonnees[iHighPos].second - dTarget);
//     // Rcpp::Rcout << "\n std::fabs(dLowSum - dTarget): " << std::fabs(dLowSum - dTarget);
//     // Rcpp::Rcout << "\n std::fabs(dLowSum + vDonnees[iHighPos].second - dTarget) gt std::numeric_limits<double>::epsilon(): " << (std::fabs(dLowSum + vDonnees[iHighPos].second - dTarget) > std::numeric_limits<double>::epsilon());
//     
//     if(std::fabs(dLowSum - dTarget) < std::numeric_limits<double>::epsilon())  // le vecteur de gauche a atteinds la cible
//     {
//       Rcpp::Rcout << "\n cAS 3\n";
//       double dNextHighest = std::min_element(vDonnees.begin() + iHighPos, vDonnees.begin() + iLength - 1)->first;
//       return (vDonnees[iHighPos - 1].first + dNextHighest) / 2;
//     }
//     
//     if( dLowSum < dTarget && ((dLowSum + vDonnees[iHighPos].second) > dTarget) )
//     {
//       Rcpp::Rcout << "\n cAS 1\n";
//       double dNextHighest = std::min_element(vDonnees.begin() + iHighPos, vDonnees.begin() + iLength - 1)->first;
//       return (vDonnees[iHighPos - 1].first + dNextHighest) / 2;
//     }
//     
//     if( dLowSum > dTarget )
//     {
//       iHigh = iHighPos - 1;
//       Rcpp::Rcout << "\n cAS 2\n";
//     }
//     else if( dLowSum < dTarget && dLowSum + vDonnees[iHighPos].second < dTarget)
//     {
//       iLow = iHighPos + 1;
//       dTarget -= (dLowSum + vDonnees[iHighPos].second);
//       Rcpp::Rcout << "\n cAS 4\n";
//     }
//     else
//     {
//       Rcpp::Rcout << "\n cAS 5\n";
//     }
//   }
// }

/*
 * ATTENTION: effet de bord: les 2 vecteurs en parametre sont modifies par la fonction
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
 *    - vPonderation:  NumericVector contenant la liste des ponderations
 *    - vQuantiles:    NumericVector contenant la liste des quantiles a determiner
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
  
  // debut des controles de validite
  if(iNbModalites != iNbPonderations)
  {
    Rcpp::Rcerr << "Il doit y avoir autant de modalites que de ponderations. NbModalites: " << iNbModalites << " - NbPonderations" << iNbPonderations << "\n";
    throw ""; 
  } 
  
  if(iNbModalites == 0)
  {
    Rcpp::Rcerr << "Il doit y avoir au moins une modalite";
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
  
  // timespec ts_beg, ts_end;
  // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_beg);
  quickSort(vModalites, vPonderation, 0, iNbModalites - 1);
  // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);
  // dAccumulateurTempsTri += (ts_end.tv_sec - ts_beg.tv_sec) + (ts_end.tv_nsec - ts_beg.tv_nsec) / 1e9;
  
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
        {// dernier element du tableau
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
  ,NumericVector vYobservations
  ,int iRayon
  ,NumericMatrix mVar
  ,NumericVector vXCentroides
  ,NumericVector vYCentroides
  ,NumericVector vQuantiles
)
{
  // debut benchmark
  timespec ts_beg, ts_end;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_beg);
  double dTempsPasse;
  double dTempsTotal = 0;
  int iTempsRestant = 0;
  int iPourcentageEffectue;
  int iPourcentageEffectuePrecedent = 0;
  // fin benchmark
  
  int i;
  int iVarCourante;
  int iIndiceObsCourante;
  int iNbVars = mVar.ncol();                // nombre de variables a traiter
  int iNbObs = vXobservations.length();     // nombre d'observations
  int iNbCentroides = vXCentroides.length();     // nombre de centroides
  int iNbQuantiles = vQuantiles.length();
  long double dRayonCarre = pow(iRayon, 2); // rayon de lissage au carre
  long double dDistanceCarre;               // distance au carre entre une observation et un centroide
  
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
      dDistanceCarre = pow(long(vXobservations[iIndiceObsCourante]) - vXCentroides(iIndiceCentroide), 2) + pow(long(vYobservations[iIndiceObsCourante]) - vYCentroides(iIndiceCentroide), 2);
      
      if (dDistanceCarre < dRayonCarre)
      {
        dPonderation = pow(1 - (dDistanceCarre / dRayonCarre), 2);
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
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);
    dTempsPasse = (ts_end.tv_sec - ts_beg.tv_sec) + (ts_end.tv_nsec - ts_beg.tv_nsec) / 1e9;
    iPourcentageEffectue = 100 * iIndiceCentroide / iNbCentroides;
    if(iPourcentageEffectuePrecedent != iPourcentageEffectue)  
    {
      dTempsTotal = dTempsPasse * 100 / iPourcentageEffectue;
      iTempsRestant = ceil(dTempsTotal - dTempsPasse);
      iPourcentageEffectuePrecedent = iPourcentageEffectue;
      Rcpp::Rcout << "\rMedian smoothing progress: " << iPourcentageEffectue << "% - remaining time: " << floor(iTempsRestant / 60) << "m " << (iTempsRestant % 60) << "s                                                                     ";
    }
    // fin benchmark
  }

  // debut benchmark
  Rcpp::Rcout << "\rElapsed time: " << floor(dTempsTotal / 60) << "m " << ((int)dTempsTotal % 60)<< "s                                                                 ";
  // Rcpp::Rcoutt << "\naccumulateur tri: " << dAccumulateurTempsTri;
  // fin benchmark
  
  return(mValeursMedianes);
}
