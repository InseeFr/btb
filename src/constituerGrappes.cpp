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

using namespace Rcpp;
using namespace std;
using namespace arma;

/** Fonction récursive mettant à jour la matrice mGrappe qui contient l'identifiant de grappe pour chaque case 
 *      Permet, entre autres, de constituer des grappes respectant le secret statistique
 *  Remarque: les grappes avec 0 observation sont autorisées
 *
 * algorithme top-down - variante du quadTree: 
 *  - Si les 4 sous-carrés contiennent 
 *        -    au moins iNbObsMin observations 
 *        - ou 0 observation
 *      => alors découper en 4
 *  - les zones avec 0 observation portent le numéro de grappe "-1"
 *
 *  arguments
 *    - iNbObsMin          : nombre minimum d'observation que doit contenir chaque grappe. 
 *    - profondeurMax      : nombre de fois maximum que peut être découpé la matrice mEffectifs2n en 4
 *    - mEffectifs2n       : matrice des effectifs - doit être carrée et la longueur de son côté doit être une puissance de 2
 *    - mGrappe            : matrice destiné à recevoir le numéro de grappe pour chaque case - cette matrice est de la meme taille que mEffectifs2n
 *    - vNoGrappe          : vecteur destiné à recevoir la liste des numéros de grappes présents dans la matrice résultat (utile pour un traitement en parallèle ultérieur)
 *    - profondeurCourante : nombre de fois que la sous-zone a été découpée
 *    - iTaille            : longueur en nombre de case du carré considéré
 *    - iRowReference      : numéro d'indice ligne correspondant à la case Nord-ouest du sous-carré considéré
 *    - iColReference      : numéro d'indice colonne correspondant à la case Nord-ouest du sous-carré considéré
 *  
 *  retourne : void (le résultat est renseigné par effet de bord 
 *    - dans la matrice mGrappe 
 *    - dans le vecteur vNoGrappe
 *      
 */
void quadTree2(  const unsigned int iNbObsMin
               , const unsigned short profondeurMax
               , const arma::Mat<int>& mEffectifs2n
               , arma::Mat<int>& mGrappe
               , std::vector<int> &vNoGrappe
               , const unsigned short profondeurCourante
               , const unsigned int iTaille
               , const unsigned int iRowReference
               , const unsigned int iColReference
)
{
  // Rcout << std::endl << "profondeurMax: " << profondeurMax  << " - profondeurCourante: " << profondeurCourante << " - iTaille: " << iTaille << " - iRowRef: " << iRowReference << " - iColRef: " << iColReference << " - Grappe.size(): " << vNoGrappe.size() << std::endl;
  // 
  // Rcout << "NoGrappes: ";
  // for(std::vector<int>::iterator it = vNoGrappe.begin(); it < vNoGrappe.end(); ++it)
  //   Rcout << " " << *it;
  // 
  // Rcout << std::endl;
  
  if (iTaille == 1)
    return;
  
  const unsigned int iDemiTaille = iTaille / 2;
  unsigned long int iNbObsClusterNO = 0;
  unsigned long int iNbObsClusterNE = 0;
  unsigned long int iNbObsClusterSO = 0;
  unsigned long int iNbObsClusterSE = 0;
  unsigned long int iRow, iCol, iRowMax, iColMax, iColMin;
  
  // Etape 1: vérifier qu'il y a suffisamment d'observations dans chaque sous-zone pour découper en 4
  
  // Nord - ouest
  for(iRow = iRowReference; iNbObsClusterNO < iNbObsMin && iRow < iRowReference + iDemiTaille; ++iRow)
    for(iCol = iColReference; iNbObsClusterNO < iNbObsMin && iCol < iColReference + iDemiTaille; ++iCol)
      iNbObsClusterNO += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterNO < iNbObsMin && iNbObsClusterNO != 0)
    return;
  
  // Nord - est
  for(iRow = iRowReference; iNbObsClusterNE < iNbObsMin && iRow < iRowReference + iDemiTaille; ++iRow)
    for(iCol = iColReference + iDemiTaille; iNbObsClusterNE < iNbObsMin && iCol < iColReference + iTaille; ++iCol)
      iNbObsClusterNE += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterNE < iNbObsMin && iNbObsClusterNE != 0)
    return;
  
  // Sud - ouest
  for(iRow = iRowReference + iDemiTaille; iNbObsClusterSO < iNbObsMin && iRow < iRowReference + iTaille; ++iRow)
    for(iCol = iColReference; iNbObsClusterSO < iNbObsMin && iCol < iColReference + iDemiTaille; ++iCol)
      iNbObsClusterSO += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterSO < iNbObsMin && iNbObsClusterSO != 0)
    return;
  
  // Sud - est
  for(iRow = iRowReference + iDemiTaille; iNbObsClusterSE < iNbObsMin && iRow < iRowReference + iTaille; ++iRow)
    for(iCol = iColReference + iDemiTaille; iNbObsClusterSE < iNbObsMin && iCol < iColReference + iTaille; ++iCol)
      iNbObsClusterSE += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterSE < iNbObsMin && iNbObsClusterSE != 0)
    return;
  
  // Rcout << "NO: " << iNbObsClusterNO << " - NE: " << iNbObsClusterNE << " - SO: " << iNbObsClusterSO << " - SE: " << iNbObsClusterSE << std::endl;
  
  // Etape 2: On supprime l'ancien numero de grappe
  if(std::find(vNoGrappe.begin(), vNoGrappe.end(), mGrappe(iRowReference, iColReference)) != vNoGrappe.end())
    vNoGrappe.erase(std::find(vNoGrappe.begin(), vNoGrappe.end(), mGrappe(iRowReference, iColReference)));
  
  // Etape 3: On enregistre les numeros de grappe
  long int lValNordOuest = std::pow(4, profondeurMax - profondeurCourante);
  long int lValNordEst = mGrappe(iRowReference, iColReference) + lValNordOuest;
  long int lValSudOuest = mGrappe(iRowReference + iDemiTaille, iColReference) + 2 * lValNordOuest;
  long int lValSudEst = mGrappe(iRowReference + iDemiTaille, iColReference + iDemiTaille) + 3 * lValNordOuest;
  
  if(iNbObsClusterNO != 0)
    if (std::find(vNoGrappe.begin(), vNoGrappe.end(), mGrappe(iRowReference, iColReference)) == vNoGrappe.end())
      vNoGrappe.push_back(mGrappe(iRowReference, iColReference));
    
  if(iNbObsClusterNE != 0)
    if (std::find(vNoGrappe.begin(), vNoGrappe.end(), lValNordEst) == vNoGrappe.end())
      vNoGrappe.push_back(lValNordEst);

  if(iNbObsClusterSO != 0)
    if (std::find(vNoGrappe.begin(), vNoGrappe.end(), lValSudOuest) == vNoGrappe.end())
      vNoGrappe.push_back(lValSudOuest);
  
  if(iNbObsClusterSE != 0)
    if (std::find(vNoGrappe.begin(), vNoGrappe.end(), lValSudEst) == vNoGrappe.end())
      vNoGrappe.push_back(lValSudEst);
    
  // Rcout << "vNoGrappe.size(): " << vNoGrappe.size() << " (" << mGrappe(iRowReference, iColReference) << "," << lValNordEst << "," << lValSudOuest << "," << lValSudEst << ")" << std::endl;
  
  // Etape 4: On renumérote les 3 zones
  iRowMax = iRowReference + iTaille;
  iColMax = iColReference + iTaille;
  iColMin = iColReference + iDemiTaille;
  for(iRow = iRowReference + iDemiTaille; iRow < iRowMax; ++iRow)
  {
    for(iCol = iColMin; iCol < iColMax; ++iCol)
    {
      if(iNbObsClusterNO == 0)
        mGrappe(iRow - iDemiTaille, iCol - iDemiTaille) = -1;
      
      if(iNbObsClusterNE == 0)
        mGrappe(iRow - iDemiTaille, iCol              ) = -1;
      else
        mGrappe(iRow - iDemiTaille, iCol              ) = lValNordEst;
      
      if(iNbObsClusterSO == 0)
        mGrappe(iRow              , iCol - iDemiTaille) = -1;
      else
        mGrappe(iRow              , iCol - iDemiTaille) = lValSudOuest;
      
      if(iNbObsClusterSE == 0)
        mGrappe(iRow              , iCol              ) = -1;
      else
        mGrappe(iRow              , iCol              ) = lValSudEst;
    }  
  }    
  
  // Rcout << mGrappe;
  
  // Etape 5: appels récursifs sur les 4 sous-zones
  if(iNbObsClusterNO != 0)
    quadTree2(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference              , iColReference);
  
  if(iNbObsClusterNE != 0)
    quadTree2(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference              , iColReference + iDemiTaille);
  
  if(iNbObsClusterSO != 0)
    quadTree2(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference + iDemiTaille, iColReference);
  
  if(iNbObsClusterSE != 0)
    quadTree2(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference + iDemiTaille, iColReference + iDemiTaille);
}

/** Fonction permettant de constituer des grappes de carreaux
*  Remarques
*    - version de la fonction qui créée une grappe unique, portant le numéro -1, pour toutes les cases sans observation
*    - la taille de la matrice en entrée qui doit être carrée et la longueur du côté doit être une puissance de 2 
*    - la taille de la matrice en entrée est limitée à 32768 x 32768 pour éviter l'overflow dans les integer en R 2^32
*    - utiliser une StringMatrix à la place d'une IntegerMatrix pour mGrappe n'est pas une bonne idée car le traitement est 24 fois plus long
* 
*  arguments
*    - iNbObsMin    : nombre minimum d'observation que doit contenir chaque grappe. 
*                     Permet de constituer des grappes respectant le secret statistique
*                     Remarque: les grappes avec 0 observation sont autorisées
*    - mEffectifs2n : matrice des effectifs - doit être carrée et la longueur de son côté doit être une puissance de 2
*    - vNoGrappe    : vecteur destiné à recevoir la liste des numéros de grappes présents dans la matrice résultat (utile pour un traitement en parallèle ultérieur)
*  
*  retourne  : une matrice carree de la meme taille que mEffectifs2n dont chaque case contient l'identifiant de la grappe à laquelle elle a été affectée
*              par effet de bord, vNoGrappe est également retourné
*  
*/
arma::Mat<int> constituerGrappes2(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs, std::vector<int>& vNoGrappe) 
{
  // Rcout << mEffectifs;
  
  const unsigned int iLongueur = mEffectifs.n_cols;
  const unsigned short profondeurMax = std::ceil(log( iLongueur ) / log(2));
  arma::Mat<int> mGrappe(iLongueur, iLongueur, fill::zeros);
  
  if(profondeurMax >= 16)
  {
    Rcpp::Rcerr << "makeClusterObsMin does not allow matrix greater than 32768 x 32768 to prevent overflow in R integer 2^32"; 
    return mGrappe;
  }
  
  try
  {
    quadTree2(iNbObsMin, profondeurMax, mEffectifs, mGrappe, vNoGrappe, 1, iLongueur, 0, 0);
    if(vNoGrappe.size() == 0)
      vNoGrappe.push_back(0);
  }
  catch(const std::exception & e)
  {
    Rcpp::Rcerr << "Exception in makeClusterObsMin: " << e.what(); 
  }
  
  return mGrappe;
}

/** Fonction permettant de constituer des grappes de carreaux
 * 
 * version :
 *  - destinée à une utilisation depuis R; exposée par le package; on ne peut pas récupérer la liste des numéros de grappe via sa référence
 *  - toutes les cases sans observation ont le même numéro de grappe
 * 
 * surcharge makeClusterObsMin1 const unsigned int iNbObsMin const arma::Mat int & mEffectifs std::vector int & vNoGrappe 
 * 
 */
// // [[Rcpp::export]]
arma::Mat<int> constituerGrappes2(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs)
{
  vector<int> vInt (1,0);
  return constituerGrappes2(iNbObsMin, mEffectifs, vInt);
}

/** Fonction récursive mettant à jour la matrice mGrappe qui contient l'identifiant de grappe pour chaque case 
 *      Permet, entre autres, de constituer des grappes respectant le secret statistique
 *  Remarque: les grappes avec 0 observation sont autorisées
 *
 * algorithme top-down - variante du quadTree: 
 *  - Si les 4 sous-carrés contiennent 
 *        -    au moins iNbObsMin observations 
 *        - ou 0 observation
 *      => alors découper en 4
 *  - On ne redécoupe pas les sous-carrés  
 *        -    contenant 0 observation
 *        - ou dont un carreau contient toutes les observations
 *
 *  arguments
 *    - iNbObsMin          : nombre minimum d'observation que doit contenir chaque grappe. 
 *    - profondeurMax      : nombre de fois maximum que peut être découpé la matrice mEffectifs2n en 4
 *    - mEffectifs2n       : matrice des effectifs - doit être carrée et la longueur de son côté doit être une puissance de 2
 *    - mGrappe            : matrice destiné à recevoir le numéro de grappe pour chaque case - cette matrice est de la meme taille que mEffectifs2n
 *    - vNoGrappe          : vecteur destiné à recevoir la liste des numéros de grappes présents dans la matrice résultat (utile pour un traitement en parallèle ultérieur)
 *    - profondeurCourante : nombre de fois que la sous-zone a été découpée
 *    - iTaille            : longueur en nombre de case du carré considéré
 *    - iRowReference      : numéro d'indice ligne correspondant à la case Nord-ouest du sous-carré considéré
 *    - iColReference      : numéro d'indice colonne correspondant à la case Nord-ouest du sous-carré considéré
 *  
 *  retourne : void (le résultat est renseigné par effet de bord 
 *    - dans la matrice mGrappe 
 *    - dans le vecteur vNoGrappe
 *      
 */
void quadTree( const unsigned int iNbObsMin
             , const unsigned short profondeurMax
             , const arma::Mat<int>& mEffectifs2n
             , arma::Mat<int>& mGrappe
             , std::vector<int> &vNoGrappe
             , const unsigned short profondeurCourante
             , const unsigned int iTaille
             , const unsigned int iRowReference
             , const unsigned int iColReference
)
{
  // Rcout << std::endl << "profondeurMax: " << profondeurMax  << " - profondeurCourante: " << profondeurCourante << " - iTaille: " << iTaille << " - iRowRef: " << iRowReference << " - iColRef: " << iColReference << " - Grappe.size(): " << vNoGrappe.size() << std::endl;
  // 
  // Rcout << "NoGrappes: ";
  // for(std::vector<int>::iterator it = vNoGrappe.begin(); it < vNoGrappe.end(); ++it)
  //   Rcout << " " << *it;
  // 
  // Rcout << std::endl;
  
  if (iTaille == 1)
    return;
  
  const unsigned int iDemiTaille = iTaille / 2;
  unsigned long int iNbObsClusterNO = 0;
  unsigned long int iNbObsClusterNE = 0;
  unsigned long int iNbObsClusterSO = 0;
  unsigned long int iNbObsClusterSE = 0;
  unsigned long int iRow, iCol, iRowMax, iColMax, iColMin;
  
  // Etape 1: vérifier qu'il y a suffisamment d'observations dans chaque sous-zone pour découper en 4
  
  // Nord - ouest
  for(iRow = iRowReference; iNbObsClusterNO < iNbObsMin && iRow < iRowReference + iDemiTaille; ++iRow)
    for(iCol = iColReference; iNbObsClusterNO < iNbObsMin && iCol < iColReference + iDemiTaille; ++iCol)
      iNbObsClusterNO += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterNO < iNbObsMin && iNbObsClusterNO != 0)
    return;
  
  // Nord - est
  for(iRow = iRowReference; iNbObsClusterNE < iNbObsMin && iRow < iRowReference + iDemiTaille; ++iRow)
    for(iCol = iColReference + iDemiTaille; iNbObsClusterNE < iNbObsMin && iCol < iColReference + iTaille; ++iCol)
      iNbObsClusterNE += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterNE < iNbObsMin && iNbObsClusterNE != 0)
    return;
  
  // Sud - ouest
  for(iRow = iRowReference + iDemiTaille; iNbObsClusterSO < iNbObsMin && iRow < iRowReference + iTaille; ++iRow)
    for(iCol = iColReference; iNbObsClusterSO < iNbObsMin && iCol < iColReference + iDemiTaille; ++iCol)
      iNbObsClusterSO += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterSO < iNbObsMin && iNbObsClusterSO != 0)
    return;
  
  // Sud - est
  for(iRow = iRowReference + iDemiTaille; iNbObsClusterSE < iNbObsMin && iRow < iRowReference + iTaille; ++iRow)
    for(iCol = iColReference + iDemiTaille; iNbObsClusterSE < iNbObsMin && iCol < iColReference + iTaille; ++iCol)
      iNbObsClusterSE += mEffectifs2n(iRow, iCol);
  
  if(iNbObsClusterSE < iNbObsMin && iNbObsClusterSE != 0)
    return;
  
  // On ne découpe pas les sous-zones contenant 0 obs ou dont un des 4 quartiers contient toutes les observations
  const unsigned long int iTotalGrappe = iNbObsClusterNO + iNbObsClusterNE + iNbObsClusterSO + iNbObsClusterSE;
  if(iNbObsClusterNO == iTotalGrappe || iNbObsClusterNE == iTotalGrappe || iNbObsClusterSO == iTotalGrappe || iNbObsClusterSE == iTotalGrappe)
    return;
  
  // Rcout << "NO: " << iNbObsClusterNO << " - NE: " << iNbObsClusterNE << " - SO: " << iNbObsClusterSO << " - SE: " << iNbObsClusterSE << std::endl;
  
  // Etape 2: On supprime l'ancien numero de grappe
  if(std::find(vNoGrappe.begin(), vNoGrappe.end(), mGrappe(iRowReference, iColReference)) != vNoGrappe.end())
    vNoGrappe.erase(std::find(vNoGrappe.begin(), vNoGrappe.end(), mGrappe(iRowReference, iColReference)));
  
  // Etape 3: On enregistre les numeros de grappe
  long int lValNordOuest = std::pow(4, profondeurMax - profondeurCourante);
  long int lValNordEst = mGrappe(iRowReference, iColReference) + lValNordOuest;
  long int lValSudOuest = mGrappe(iRowReference + iDemiTaille, iColReference) + 2 * lValNordOuest;
  long int lValSudEst = mGrappe(iRowReference + iDemiTaille, iColReference + iDemiTaille) + 3 * lValNordOuest;
  
  vNoGrappe.push_back(mGrappe(iRowReference, iColReference));
  vNoGrappe.push_back(lValNordEst);
  vNoGrappe.push_back(lValSudOuest);
  vNoGrappe.push_back(lValSudEst);
  
  // Rcout << "vNoGrappe.size(): " << vNoGrappe.size() << " (" << mGrappe(iRowReference, iColReference) << "," << lValNordEst << "," << lValSudOuest << "," << lValSudEst << ")" << std::endl;
  
  // Etape 4: On renumérote les 3 zones
  iRowMax = iRowReference + iTaille;
  iColMax = iColReference + iTaille;
  iColMin = iColReference + iDemiTaille;
  for(iRow = iRowReference + iDemiTaille; iRow < iRowMax; ++iRow)
    for(iCol = iColMin; iCol < iColMax; ++iCol)
    {
      mGrappe(iRow - iDemiTaille, iCol              ) = lValNordEst;
      mGrappe(iRow              , iCol - iDemiTaille) = lValSudOuest;
      mGrappe(iRow              , iCol              ) = lValSudEst;
    }  
    
  // Etape 5: appels récursifs sur les 4 sous-zones
  quadTree(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference              , iColReference);
  quadTree(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference              , iColReference + iDemiTaille);
  quadTree(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference + iDemiTaille, iColReference);
  quadTree(iNbObsMin, profondeurMax, mEffectifs2n, mGrappe, vNoGrappe, profondeurCourante + 1, iDemiTaille, iRowReference + iDemiTaille, iColReference + iDemiTaille);
}

/** Fonction permettant de constituer des grappes de carreaux
*  la taille de la matrice en entrée qui doit être carrée et dont la longueur du côté doit être une puissance de 2 
*  est limitée à 32768 x 32768 pour éviter l'overflow dans les integer en R 2^32
* 
*  arguments
*  iNbObsMin    : nombre minimum d'observation que doit contenir chaque grappe. 
*                  Permet de constituer des grappes respectant le secret statistique
*                  Remarque: les grappes avec 0 observation sont autorisées
*  mEffectifs2n : matrice des effectifs - doit être carrée et la longueur de son côté doit être une puissance de 2
*  
*  retourne  : une matrice carree de la meme taille que mEffectifs2n 
*              dont chaque case contient l'identifiant de la grappe à laquelle elle a été affectée
*              
*  Remarque: utiliser une StringMatrix à la place d'une IntegerMatrix pour mGrappe n'est pas une bonne idée 
*  car le traitement est 24 fois plus long
*/
arma::Mat<int> constituerGrappes(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs, std::vector<int> &vNoGrappe) 
{
  // Rcout << mEffectifs;
  
  const unsigned int iLongueur = mEffectifs.n_cols;
  const unsigned short profondeurMax = std::ceil(log( iLongueur ) / log(2));
  arma::Mat<int> mGrappe(iLongueur, iLongueur, fill::zeros);
  
  if(profondeurMax >= 16)
  {
    Rcpp::Rcerr << "makeClusterObsMin does not allow matrix greater than 32768 x 32768 to prevent overflow in R integer 2^32"; 
    return mGrappe;
  }
  
  try
  {
    quadTree(iNbObsMin, profondeurMax, mEffectifs, mGrappe, vNoGrappe, 1, iLongueur, 0, 0);
    if(vNoGrappe.size() == 0)
      vNoGrappe.push_back(0);
  }
  catch(const std::exception & e)
  {
    Rcpp::Rcerr << "Exception in makeClusterObsMin: " << e.what(); 
  }
  
  // Rcout << mGrappe;
  return mGrappe;
}

/** Fonction permettant de constituer des grappes de carreaux
 * 
 * version destinée à une utilisation depuis R; exposée par le package; on ne peut pas récupérer la liste des numéros de grappe via sa référence  
 * 
 * surcharge constituerGrappes const unsigned int iNbObsMin const arma::Mat int & mEffectifs std::vector int & vNoGrappe 
 * 
 */
// [[Rcpp::export]]
arma::Mat<int> constituerGrappes(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs) 
{
  vector<int> vInt (1,0);
  return constituerGrappes(iNbObsMin, mEffectifs, vInt);
}

/** Fonction récursive permettant de retrouver les coordonnées du point en haut à gauche (Nord-est) d'un identifiant de grappe en fonction de son numéro
 * 
 * iNiveauMax : entier représentant l'exposant de la puissance de 2 pour obtenir la taille du côté de la matrice contenant les numeros de grappe 
 *              exemple: si la matrice est 128 * 128, alors iNiveauMax = log2(128) = 7
 * iNiveaucourant : niveau de décomposition (valeur comprise entre 1 et iNiveauMax)
 * iNombre : numéro de grappe à décomposer
 * vIndice : vecteur d'entiers de taille 2 qui contient le résultat 
 *          - Le premier élément est la coordonnée numéro de ligne - commence à 0
 *          - Le deuxième élément est la coordonnée numéro de colonne - commence à 0
 * 
 * retourne : void (résultat dans vIndice modifié par effet de bord)
 * 
 */
void decomposer(const int iNiveauMax, const int iNiveaucourant, const int iNombre, std::vector<int> &vIndice)
{
  // Rcout << "iNiveauMax: " << iNiveauMax  << " - iNiveaucourant: " << iNiveaucourant  << " - iNombre: " << iNombre << std::endl;
  
  if(iNiveaucourant > iNiveauMax)
    return;
  
  int iPositionNiveau = floor(iNombre / pow(4, iNiveauMax - iNiveaucourant));
  
  if(iPositionNiveau == 1)
    vIndice[1] += pow(2, iNiveauMax - iNiveaucourant);
  else if(iPositionNiveau == 2)
    vIndice[0] += pow(2, iNiveauMax - iNiveaucourant);
  else if(iPositionNiveau == 3)
  {
    vIndice[0] += pow(2, iNiveauMax - iNiveaucourant);
    vIndice[1] += pow(2, iNiveauMax - iNiveaucourant);
  }
  decomposer(iNiveauMax, iNiveaucourant + 1, iNombre - iPositionNiveau * pow(4, iNiveauMax - iNiveaucourant), vIndice);
}

/** Fonction permettant de retrouver les coordonnées du point en haut à gauche (Nord-est) d'un identifiant de grappe en fonction de son numéro
 * 
 * arguments
 * iNiveauMax : entier représentant l'exposant de la puissance de 2 pour obtenir la taille du côté de la matrice contenant les numeros de grappe 
 *              exemple: si la matrice est 128 * 128, alors iNiveauMax = log2(128) = 7
 * iNoGrappe  : numéro de grappe à décomposer
 *  
 * retourne: un vecteur d'entiers de taille 2. 
 *          - Le premier élément est la coordonnée numéro de ligne - commence à 0
 *          - Le deuxième élément est la coordonnée numéro de colonne - commence à 0
 *
 */
// [[Rcpp::export]]
std::vector<int> coordonneesGrappe(int iNiveauMax, int iNoGrappe)
{
  std::vector<int> vIndice(2);
  decomposer(iNiveauMax, 1, iNoGrappe, vIndice);
  
  // Rcout << "coordonneesGrappe: (" << vIndice[0] << ";"<< vIndice[1] << ")" << std::endl;  return vIndice;
  return vIndice;
}
