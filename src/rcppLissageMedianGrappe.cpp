/*
 *    date        version         auteur              commentaire
 * 2017/08/23      0.0.1      Arlindo Dos Santos      version initiale
 *                                                    
 */

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>       
#include <RcppParallel.h>
#include <vector>
#include <armadillo>
#include <math.h>
#include "lissageMedian.h"
#include "constituerGrappes.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;
using namespace arma;

// Lissage median par grappe
struct LissageMedianGrappe : public Worker
{
  // attributs de la classe
  const std::vector<int> vNoGrappes;
  const arma::Mat<int> mGrappes;
  const arma::Mat<int> mEffectifs2n;
  const RVector<int> vXObservations;
  const RVector<int> vYObservations;
  const RVector<int> vLigneObservations;
  const RVector<int> vColonneObservations;
  const int iPas;
  const int iRayon;
  const RMatrix<double> mVar;
  const RVector<int> vXCentroides;
  const RVector<int> vYCentroides;
  const RVector<int> vLigneCentroides;
  const RVector<int> vColonneCentroides;
  const std::vector<double> vQuantile;
  int iCompteurAvancement;
  clock_t timeBegin;
  arma::mat mResultatFinal;
  
  // constructeur par defaut
  LissageMedianGrappe(
      const std::vector<int> vNoGrappes
    , const arma::Mat<int> mGrappes
    , const arma::Mat<int> mEffectifs2n
    , const RVector<int> vXObservations
    , const RVector<int> vYObservations
    , const RVector<int> vLigneObservations
    , const RVector<int> vColonneObservations
    , const int iPas
    , const int iRayon
    , const RMatrix<double> mVar
    , const RVector<int> vXCentroides
    , const RVector<int> vYCentroides
    , const RVector<int> vLigneCentroides
    , const RVector<int> vColonneCentroides
    , const std::vector<double> vQuantile
  , arma::mat mResultatFinal
  )
    :
     vNoGrappes(vNoGrappes)
    , mGrappes(mGrappes)
    , mEffectifs2n(mEffectifs2n)
    , vXObservations(vXObservations)
    , vYObservations(vYObservations)
    , vLigneObservations(vLigneObservations)
    , vColonneObservations(vColonneObservations)
    , iPas(iPas)
    , iRayon(iRayon)
    , mVar(mVar)
    , vXCentroides(vXCentroides)
    , vYCentroides(vYCentroides)
    , vLigneCentroides(vLigneCentroides)
    , vColonneCentroides(vColonneCentroides)
    , vQuantile(vQuantile)
    , iCompteurAvancement(0)
    , timeBegin(clock())
    , mResultatFinal(mResultatFinal)
  {}
  
  void operator()(std::size_t begin, std::size_t end)
  {
    // debut benchmark
    std::stringstream message;
    int iPourcentageEffectue = 0;
    double dTempsPasse;
    double dTempsTotal = 0;
    int iTempsRestant;
    // fin benchmark
    
    const unsigned int iNbCentroides = vXCentroides.length(); // nombre de centroides total
    const unsigned int iNbObs = vXObservations.length(); // nombre d'observations total
    const unsigned int iNbVar = mVar.ncol(); // nombre de variables
    const unsigned int iNbQuantiles = vQuantile.size(); // nombre de quantiles
    
    const unsigned int iMaxLigne = *max_element(vLigneObservations.begin(), vLigneObservations.end()) + 1;
    const unsigned int iMaxColonne = *max_element(vColonneObservations.begin(), vColonneObservations.end()) + 1;
    const unsigned int iMaxDimension = std::max(iMaxLigne, iMaxColonne);
    const unsigned int iTailleMatrice = pow(2, std::ceil( log(iMaxDimension) / log(2.0)));
    
    unsigned int iCompteur;
    unsigned int iLigneCourante;
    unsigned int iColonneCourante;
    unsigned int iVarCourante;
    unsigned int iNoGrappe;
    unsigned int iObsCourante;
    
    unsigned int iLigneMinCentroideGrappe;
    unsigned int iColonneMinCentroideGrappe;
    unsigned int iLigneMaxCentroideGrappe;
    unsigned int iColonneMaxCentroideGrappe;
    unsigned int iNbCentroidesGrappe;
    unsigned int iNbMaxCentroidesGrappe;
    
    unsigned int iLigneMinCentroideGrappeEtendue;
    unsigned int iColonneMinCentroideGrappeEtendue;
    unsigned int iLigneMaxCentroideGrappeEtendue;
    unsigned int iColonneMaxCentroideGrappeEtendue;
    unsigned int iNbObsGrappeEtendue;
    
    for(std::size_t i = begin; i < end; ++i)
    {
      iNoGrappe = vNoGrappes[i];
      // Rcout << std::endl << "iNoGrappe: " << iNoGrappe << std::endl;
      
      vector<int> vCoordGrappe = coordonneesGrappe(log2(iTailleMatrice), iNoGrappe);
      iLigneMinCentroideGrappe = vCoordGrappe[0];
      iColonneMinCentroideGrappe = vCoordGrappe[1];
      
      // Rcout << std::endl << "vCoordGrappe: " << iLigneMinCentroideGrappe << " - " << iColonneMinCentroideGrappe << std::endl;
      
      // coordonnées extremes de la grappe considérée
      iCompteur = 0;
      while(   iLigneMinCentroideGrappe + pow(2, iCompteur) < iTailleMatrice
                 && iNoGrappe == (unsigned int)mGrappes(iLigneMinCentroideGrappe + pow(2, iCompteur), iColonneMinCentroideGrappe)
      )
      {
        ++iCompteur;
      }
      // Rcout << "iCompteur " << iCompteur << std::endl;
      
      iLigneMaxCentroideGrappe = iLigneMinCentroideGrappe + pow(2, iCompteur) - 1;
      iColonneMaxCentroideGrappe = iColonneMinCentroideGrappe + pow(2, iCompteur) - 1;
      
      // Rcout << "grappe " << iNoGrappe << " row: " << iLigneMinCentroideGrappe << "->" << iLigneMaxCentroideGrappe << " - col: " << iColonneMinCentroideGrappe << "->" << iColonneMaxCentroideGrappe << std::endl;
      
      // compter et récupérer les centroides de la grappe
      
      // créer les vecteurs de stockage des centroides de la grappe
      // Remarque: il se peut que nous ne remplissions pas tout le vecteur car certains centroides peuvent être manquants
      iNbMaxCentroidesGrappe = (iLigneMaxCentroideGrappe - iLigneMinCentroideGrappe + 1) * (iColonneMaxCentroideGrappe - iColonneMinCentroideGrappe + 1);
      
      // Rcout << "iNbMaxCentroidesGrappe: " << iNbMaxCentroidesGrappe << std::endl;
      
      std::vector<int> vXCentroidesGrappe(iNbMaxCentroidesGrappe);
      std::vector<int> vYCentroidesGrappe(iNbMaxCentroidesGrappe);
      
      iNbCentroidesGrappe = 0;
      for(iCompteur = 0; iCompteur < iNbCentroides; ++iCompteur)
        if(  (unsigned int)vLigneCentroides[iCompteur] >= iLigneMinCentroideGrappe   
          && (unsigned int)vLigneCentroides[iCompteur] <= iLigneMaxCentroideGrappe
          && (unsigned int)vColonneCentroides[iCompteur] >= iColonneMinCentroideGrappe 
          && (unsigned int)vColonneCentroides[iCompteur] <= iColonneMaxCentroideGrappe)
        {
          vXCentroidesGrappe[iNbCentroidesGrappe] = vXCentroides[iCompteur];
          vYCentroidesGrappe[iNbCentroidesGrappe] = vYCentroides[iCompteur];
          ++iNbCentroidesGrappe;
          // Rcout << "iCompteur: " << iCompteur << " - iNbCentroidesGrappe: " << iNbCentroidesGrappe << std::endl;
        }
        
      // réduire la taille des vecteurs 
      vXCentroidesGrappe.resize(iNbCentroidesGrappe);
      vYCentroidesGrappe.resize(iNbCentroidesGrappe);
      
      // for(unsigned int i = 0; i < iNbMaxCentroidesGrappe; ++i)
      //   Rcout << "vXCentroidesGrappe[" << i << "]: " << vXCentroidesGrappe[i] << std::endl;
      
      // Rcout << "iNbCentroidesGrappe: " << iNbCentroidesGrappe << std::endl;
      
      // étendre la grappe à la zone couverte par le rayon
      iLigneMinCentroideGrappeEtendue = (unsigned int)max(0, (int)iLigneMinCentroideGrappe - (int)ceil((double)iRayon / iPas - 0.5));
      iLigneMaxCentroideGrappeEtendue = min(iTailleMatrice - 1, iLigneMaxCentroideGrappe + (int)ceil((double)iRayon / iPas - 0.5));
      iColonneMinCentroideGrappeEtendue = (unsigned int)max(0, (int)iColonneMinCentroideGrappe - (int)ceil((double)iRayon / iPas - 0.5));
      iColonneMaxCentroideGrappeEtendue = min(iTailleMatrice - 1, iColonneMaxCentroideGrappe + (int)ceil((double)iRayon / iPas - 0.5));
      
      // Rcout << "grappe etendue i: " << iLigneMinCentroideGrappeEtendue << "->" << iLigneMaxCentroideGrappeEtendue << " - j: " << iColonneMinCentroideGrappeEtendue << "->" << iColonneMaxCentroideGrappeEtendue << std::endl;
      // Rcout << "avant compter le nombre d'observations de la grappe" << std::endl;
      
      // compter le nombre d'observations de la grappe étendue
      iNbObsGrappeEtendue = 0;
      for(iLigneCourante = iLigneMinCentroideGrappeEtendue; iLigneCourante <= iLigneMaxCentroideGrappeEtendue; ++iLigneCourante)
        for(iColonneCourante = iColonneMinCentroideGrappeEtendue; iColonneCourante <= iColonneMaxCentroideGrappeEtendue; ++iColonneCourante)
          iNbObsGrappeEtendue += mEffectifs2n(iLigneCourante, iColonneCourante);
      
      // Rcout << "iNbObsGrappeEtendue: " << iNbObsGrappeEtendue << std::endl;
      
      // créer les vecteurs de stockage des observations de la grappe étendue
      std::vector<int> vXObservationsGrappe(iNbObsGrappeEtendue);
      std::vector<int> vYObservationsGrappe(iNbObsGrappeEtendue);
      arma::mat mVarGrappe(iNbObsGrappeEtendue, iNbVar);
      
      // Rcout << "avant récupération des observations de la grappe étendue" << std::endl;
      // récupération des observations de la grappe étendue
      iObsCourante = 0;
      for(iCompteur = 0; iCompteur < iNbObs; ++iCompteur)
        if(  (unsigned int)vLigneObservations[iCompteur] >= iLigneMinCentroideGrappeEtendue
          && (unsigned int)vLigneObservations[iCompteur] <= iLigneMaxCentroideGrappeEtendue
          && (unsigned int)vColonneObservations[iCompteur] >= iColonneMinCentroideGrappeEtendue 
          && (unsigned int)vColonneObservations[iCompteur] <= iColonneMaxCentroideGrappeEtendue)
        {
          // Rcout << "iCompteur: " << iCompteur << " - iObsCourante: " << iObsCourante << std::endl;
          vXObservationsGrappe[iObsCourante] = vXObservations[iCompteur];
          vYObservationsGrappe[iObsCourante] = vYObservations[iCompteur];
          for(iVarCourante = 0; iVarCourante < iNbVar; ++iVarCourante)
            mVarGrappe(iObsCourante, iVarCourante) = mVar(iCompteur, iVarCourante);
          ++iObsCourante;
        }
        
      // Rcout << "vXObservationsGrappe ";
      // for(unsigned int i = 0; i < iNbObsGrappeEtendue; ++i)
      //   Rcout << vXObservationsGrappe[i] << " ";
      // 
      // Rcout << std::endl << "vYObservationsGrappe ";
      // for(unsigned int i = 0; i < iNbObsGrappeEtendue; ++i)
      //   Rcout << vYObservationsGrappe[i] << " ";
      // 
      // Rcout << std::endl << "vXCentroidesGrappe ";
      // for(unsigned int i = 0; i < iNbCentroidesGrappe; ++i)
      //   Rcout << vXCentroidesGrappe[i] << " ";
      // 
      // Rcout << std::endl << "vYCentroidesGrappe ";
      // for(unsigned int i = 0; i < iNbCentroidesGrappe; ++i)
      //   Rcout << vYCentroidesGrappe[i] << " ";
      // 
      // Rcout << std::endl << "mVarGrappe ";
      // for(unsigned int i = 0; i < mVarGrappe.n_rows; ++i)
      //   Rcout << mVarGrappe[i] << " ";
      
      // appel calcul médiane
      arma::mat mResultat = rcppLissageMedian(vXObservationsGrappe, vYObservationsGrappe, iRayon, mVarGrappe, vXCentroidesGrappe, vYCentroidesGrappe, vQuantile);
      
      // Rcout << mResultat;
      
      unsigned int iColonneX = iNbVar * iNbQuantiles + 1; 
      unsigned int iColonneY = iColonneX + 1; 
      
      // écrire le résultat partiel dans la matrice finale
      for(unsigned int iLigneResultat = 0; iLigneResultat < (unsigned int)mResultat.n_rows; ++iLigneResultat)
      {
        for(unsigned int iLigneResultatFinal = 0; iLigneResultatFinal < (unsigned int)mResultatFinal.n_rows; ++iLigneResultatFinal)
        {
          if(mResultat(iLigneResultat, iColonneX) == mResultatFinal(iLigneResultatFinal, iColonneX) &&
             mResultat(iLigneResultat, iColonneY) == mResultatFinal(iLigneResultatFinal, iColonneY) )
          {
            for(iColonneCourante = 0; iColonneCourante < iColonneX; ++iColonneCourante)
            {             
              mResultatFinal(iLigneResultatFinal, iColonneCourante) = mResultat(iLigneResultat, iColonneCourante);
            }
            break;
          }
        }
      }
      // debut benchmark
      ++iCompteurAvancement;
      message.str("");
      iPourcentageEffectue = 100 * iCompteurAvancement / vNoGrappes.size();
      dTempsPasse = (clock() - timeBegin) / CLOCKS_PER_SEC;
      dTempsTotal = dTempsPasse * 100 / iPourcentageEffectue;
      iTempsRestant = ceil(dTempsTotal - dTempsPasse);
      message << "\rMedian smoothing progress (parallel clusters): " << iPourcentageEffectue << "% - minimum remaining time: " << (iTempsRestant / 60) << "m " << (iTempsRestant % 60) << "s           ";
      Rcout << message.str();
      // fin benchmark
    }
    // debut benchmark
    if(iPourcentageEffectue == 100)
    {   
      message.str("");
      message << "Elapsed time median smoothing: " << floor(dTempsTotal / 60) << "m " << ((int)dTempsTotal % 60)<< "s                                                                 ";
      Rcpp::Rcout << "\n" << message.str();
    }
    // fin benchmark
  }
};


// arma::mat rcppLissageMedianGrappe(
/**
 * arguments
 * 
 * vLigneObservations : indices doivent commencer à 0
 * vColonneObservations : indices doivent commencer à 0
 * 
 * vLigneCentroides : indices doivent commencer à 0
 * vColonneCentroides : indices doivent commencer à 0
 * 
 * retourne: 
 * matrice 
 *  - lignes : centroides
 *  - colonnes: nbObs, V1q1 V1q2.. V1qn, ..., Vnq1, Vnq2, ...Vnqn, x, y
 */
// [[Rcpp::export]]
NumericMatrix rcppLissageMedianGrappe(
    int iMinObsGrappe
  , IntegerVector vXObservation
  , IntegerVector vYObservation
  , IntegerVector vLigneObservation
  , IntegerVector vColonneObservation
  , int iPas
  , int iRayon
  , NumericMatrix mVariables
  , IntegerVector vXCentroide
  , IntegerVector vYCentroide
  , IntegerVector vLigneCentroide
  , IntegerVector vColonneCentroide
  , NumericVector vQuantile
  )
{ 
    // debut cast pour threadsafety - cf https://rcppcore.github.io/RcppParallel/#thread_safety
    RVector<int> vXObservations(vXObservation);
    RVector<int> vYObservations(vYObservation);
    RVector<int> vLigneObservations(vLigneObservation);
    RVector<int> vColonneObservations(vColonneObservation);
    RVector<int> vXCentroides(vXCentroide);
    RVector<int> vYCentroides(vYCentroide);
    RVector<int> vLigneCentroides(vLigneCentroide);
    RVector<int> vColonneCentroides(vColonneCentroide);
    RMatrix<double> mVar(mVariables);
    
    std::vector<double> vQuantiles(vQuantile.length());
    for(int i = 0; i < vQuantile.length(); ++i)
      vQuantiles[i] = vQuantile[i];
    // fin cast pour threadsafety 
    
    const unsigned int iNbCentroides = vXCentroides.length(); // nombre de centroides total
    const unsigned int iNbObs = vXObservations.length(); // nombre d'observations total
    const unsigned int iNbVar = mVar.ncol(); // nombre de variables
    const unsigned int iNbQuantiles = vQuantiles.size(); // nombre de quantiles
    
    const unsigned int iMaxLigne = *max_element(vLigneObservations.begin(),vLigneObservations.end()) + 1;
    const unsigned int iMaxColonne = *max_element(vColonneObservations.begin(), vColonneObservations.end()) + 1;
    const unsigned int iMaxDimension = std::max(iMaxLigne, iMaxColonne);
    const unsigned int iTailleMatrice = pow(2, std::ceil( log(iMaxDimension) / log(2.0)));
    
    unsigned int iObsCourante;

    // Rcout << "iMaxDimension: " << iMaxDimension << " - iTailleMatrice: " << iTailleMatrice << " - iMaxLigne: " << iMaxLigne << " - iMaxColonne: " << iMaxColonne << std::endl;

    arma::Mat<int> mEffectifs2n(iTailleMatrice, iTailleMatrice, fill::zeros);// la matrice des effectifs pour chacun des carreaux
   
    // remplir la matrice des effectifs
    for(iObsCourante = 0; iObsCourante < iNbObs; ++iObsCourante)
      mEffectifs2n(vLigneObservations[iObsCourante], vColonneObservations[iObsCourante]) += 1;
    
    // for(unsigned int iLigne = 0; iLigne < iTailleMatrice; ++iLigne)
    // {
    //   Rcout << std::endl;
    //   for(unsigned int iColonne = 0; iColonne < iTailleMatrice; ++iColonne)
    //     Rcout << mEffectifs2n(iLigne, iColonne) << " ";
    // }        
    
    std::vector<int> vNoGrappe;
    arma::Mat<int> mGrappes = constituerGrappes(iMinObsGrappe, mEffectifs2n, vNoGrappe);

    Rcout << "number of clusters: " << vNoGrappe.size() << std::endl;
    // Rcout << mGrappes;
            
    // NumericMatrix mResultatFinal(iNbCentroides, iNbVar * iNbQuantiles + 2 + 1); // liste des colonnes: nbObs, V1q1 V1q2.. V1qn, ..., Vnq1, Vnq2, ...Vnqn, x, y
    arma::mat mResultatFinal(iNbCentroides, iNbVar * iNbQuantiles + 2 + 1, fill::zeros); // liste des colonnes: nbObs, V1q1 V1q2.. V1qn, ..., Vnq1, Vnq2, ...Vnqn, x, y
        
    // Rcout << "iNbObs:" << iNbObs << " - iNbCentroides:" << iNbCentroides << " - iNbVar:" << iNbVar << " - iNbQuantiles:" << iNbQuantiles << std::endl;
    
    // initialiser le tableau résultat avec les coordonnées des centroides 
    unsigned int iColonneX = iNbVar * iNbQuantiles + 1; 
    unsigned int iColonneY = iColonneX + 1; 
    for(unsigned int iLigneCourante = 0; iLigneCourante < iNbCentroides; ++iLigneCourante)
    {
      mResultatFinal(iLigneCourante, iColonneX) = vXCentroides[iLigneCourante];
      mResultatFinal(iLigneCourante, iColonneY) = vYCentroides[iLigneCourante];
    }    
    
    LissageMedianGrappe lmg(vNoGrappe, mGrappes, mEffectifs2n, vXObservations, vYObservations, vLigneObservations, vColonneObservations, iPas, iRayon, mVar, vXCentroides, vYCentroides, vLigneCentroides, vColonneCentroides, vQuantiles, mResultatFinal);
    parallelFor(0, vNoGrappe.size(), lmg);

    // // tracer résultat final
    // Rcout << std::endl << mResultatFinal;
    
    return wrap(lmg.mResultatFinal);
}
