/*
 *    date        version         auteur              commentaire
 * 2017/08/23      0.0.1      Arlindo Dos Santos      version initiale
 *                                                    
 */

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>       
#include <vector>
#include <armadillo>
#include "utils.h"
#include "lissageMedian.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

/**
* arguments:
*    - vModalites:    std::vector<double> contenant la liste des modalites
*    - vPonderation:  std::vector<double> contenant la liste des pondérations
*    - vQuantiles:    std::vector<double> contenant la liste des quantiles à déterminer
* 
* retourne
*    - std::vector<double> contenant les quantiles
*/
// [[Rcpp::export]]
std::vector<double> calculeQuantiles(std::vector<double>& vModalites, std::vector<double>& vPonderation, const std::vector<double> vQuantiles)
{
  int iNbModalites = vModalites.size();
  int iNbPonderations = vPonderation.size();
  int iNbQuantiles = vQuantiles.size();
  
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
    if(vQuantiles[i] < 0 || vQuantiles[i] > 1)
    {
      Rcpp::Rcerr << "Valeur de quantile invalide: " << vQuantiles[i];
      throw "";
    }
  }
  // fin des controles de validite
  
  std::vector<double> vPonderationCumulee(iNbPonderations);
  std::vector<double> vResultat(iNbQuantiles);
  long double dCible = 0;
  
  quickSort(vModalites, vPonderation, 0, iNbModalites - 1);
  
  // Rcout << "vModalites";
  // for(unsigned int i = 0; i < vModalites.size(); ++i)
  //   Rcout << " " << i << ": "<< vModalites[i] << " ";
  // 
  // Rcout << std::endl << "vPonderations";
  // for(unsigned int i = 0; i < vPonderation.size(); ++i)
  //   Rcout << " " << i << ": "<< vPonderation[i] << " ";
  // 
  // Rcout << std::endl << "vQuantiles";
  // for(unsigned int i = 0; i < vQuantiles.size(); ++i)
  //   Rcout << " " << i << ": "<< vQuantiles[i] << " ";
  
  vPonderationCumulee[0] = vPonderation[0];   // calcul des ponderations cumulees
  for(int iPonderation = 1; iPonderation < iNbPonderations; ++iPonderation)
    vPonderationCumulee[iPonderation] = vPonderationCumulee[iPonderation - 1] + vPonderation[iPonderation];
  
  for(unsigned int iQuantile = 0; iQuantile < vQuantiles.size(); ++iQuantile)    // pour tous les quantiles demandes
  {
    dCible = vQuantiles[iQuantile] * vPonderationCumulee[iNbPonderations - 1];
    
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

/**
* 
*/
// [[Rcpp::export]]
arma::mat rcppLissageMedian(
    std::vector<int> vXobservations
    , std::vector<int> vYobservations
    , int iRayon
    , arma::mat& mVar
    , std::vector<int> vXCentroides
    , std::vector<int> vYCentroides
    , std::vector<double> vQuantiles
)
{
  // // debut benchmark
  // clock_t timeBegin = clock();
  // double dTempsPasse;
  // double dTempsTotal = 0;
  // int iTempsRestant = 0;
  // int iPourcentageEffectue;
  // int iPourcentageEffectuePrecedent = 0;
  // std::stringstream message;
  // // fin benchmark
  
  unsigned int i;
  unsigned int iVarCourante;
  unsigned int iIndiceObsCourante;
  unsigned int iNbVars = mVar.n_cols;                // nombre de variables à traiter
  unsigned int iNbObs = vXobservations.size();     // nombre d'observations
  unsigned int iNbCentroides = vXCentroides.size();     // nombre de centroides
  unsigned int iNbQuantiles = vQuantiles.size();
  long double dRayonCarre = std::pow((long double)iRayon, 2); // rayon de lissage au carre
  long double dDistanceCarre;               // distance au carré entre une observation et un centroide
  long double dPonderation;
  
  arma::mat mValeursMedianes(iNbCentroides, iNbVars * iNbQuantiles + 1 + 2, fill::zeros); // +2 pour remettre à la fin les colonnes x et y des centroides associés
  
  // for(i = 0; i < iNbQuantiles; ++i)
  //   Rcout << "vQuantiles[" << i << "]: "<< vQuantiles[i] << std::endl;
  
  /* on parcourt tous les centroides */
  for(unsigned int iIndiceCentroide = 0; iIndiceCentroide < iNbCentroides; ++iIndiceCentroide)
  {
    std::vector<double> vIndiceObservations; 
    std::vector<double> vPonderations; 
    dPonderation = 0;
    
    // Rcout << std::endl << "iIndiceCentroide: " << iIndiceCentroide << std::endl;
    
    /* on parcourt toutes les observations */
    for(iIndiceObsCourante = 0; iIndiceObsCourante < iNbObs; ++iIndiceObsCourante)
    {
      dDistanceCarre = std::pow(long(vXobservations[iIndiceObsCourante]) - vXCentroides[iIndiceCentroide], 2) + std::pow(long(vYobservations[iIndiceObsCourante]) - vYCentroides[iIndiceCentroide], 2);
      // Rcout << "dDistanceCarre: " << dDistanceCarre << std::endl;
      
      if(dDistanceCarre < dRayonCarre)
      {
        dPonderation = std::pow(1 - (dDistanceCarre / dRayonCarre), 2);
        vPonderations.push_back(dPonderation);
        vIndiceObservations.push_back(iIndiceObsCourante);
        // Rcout << "iIndiceObsCourante: " << iIndiceObsCourante << " - dPonderation: " << dPonderation << std::endl;
      }
    }
    
    if(vIndiceObservations.size() > 0)
    {
      mValeursMedianes(iIndiceCentroide, 0) = vIndiceObservations.size(); /* nb obs ayant servi à calculer les quantiles */
      // Rcout << "mValeursMedianes(iIndiceCentroide, 0): " << mValeursMedianes(iIndiceCentroide, 0) << std::endl;
      for (iVarCourante = 0; iVarCourante < iNbVars; ++iVarCourante) /* pour chacune des variables a lisser */
      {
        std::vector<double> vModalites;
        for(i = 0; i < vIndiceObservations.size(); ++i)
          vModalites.push_back(mVar(vIndiceObservations[i], iVarCourante));
        
        std::vector<double> vValeursQuantilesVarCourante = calculeQuantiles(vModalites, vPonderations, vQuantiles);
        for(i = 0; i < iNbQuantiles; ++i)
          mValeursMedianes(iIndiceCentroide, i + iVarCourante * iNbQuantiles + 1) = vValeursQuantilesVarCourante[i];
      }
    }
    // // debut benchmark
    // dTempsPasse = (clock() - timeBegin) / CLOCKS_PER_SEC;
    // iPourcentageEffectue = 100 * iIndiceCentroide / iNbCentroides;
    // if(iPourcentageEffectuePrecedent != iPourcentageEffectue)  
    // {
    //   dTempsTotal = dTempsPasse * 100 / iPourcentageEffectue;
    //   iTempsRestant = ceil(dTempsTotal - dTempsPasse);
    //   iPourcentageEffectuePrecedent = iPourcentageEffectue;
    //   message.str("");
    //   message << "Median smoothing progress: " << iPourcentageEffectue << "% - minimum remaining time: " << (iTempsRestant / 60) << "m " << (iTempsRestant % 60) << "s";
    //   if(updateProgress.isNotNull())
    //     as<Function>(updateProgress)(iPourcentageEffectue, message.str());
    //   else
    //     Rcpp::Rcout << "\r" << message.str();
    // }
    // // fin benchmark
  }
  
  // // debut benchmark
  // message.str("");
  // message << "Elapsed time median smoothing: " << floor(dTempsTotal / 60) << "m " << ((int)dTempsTotal % 60)<< "s                                                                 ";
  // if(updateProgress.isNotNull())
  //   as<Function>(updateProgress)(iPourcentageEffectue, message.str());
  // else
  //   Rcpp::Rcout << "\n" << message.str();
  // // fin benchmark
  
  unsigned char iColonneX = mVar.n_cols * vQuantiles.size() + 1 + 0;
  unsigned char iColonneY = iColonneX + 1;
  
  int xSize = vXCentroides.size();
  for (int i = 0; i < xSize; ++i) 
  {
    mValeursMedianes(i, iColonneX) = vXCentroides[i];
    mValeursMedianes(i, iColonneY) = vYCentroides[i];
  }
  
  return(mValeursMedianes);
}
