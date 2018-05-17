/*
 *    date        version         auteur              commentaire
 * 2017/08/23      0.0.1      Arlindo Dos Santos      version initiale
 *                                                    
 */

// [[Rcpp::plugins(cpp11)]]
#include <vector>

using namespace std;

/* 
 *    date        version         auteur              commentaire
 * 2017/08/02      0.0.1      Arlindo Dos Santos      version initiale
 */

/**
 * ATTENTION: effet de bord: les 2 vecteurs en paramètre sont modifiés par la fonction
 * http://www.algolist.net/Algorithms/Sorting/Quicksort
 * 
 */
void quickSort(std::vector<double>& vModalites, std::vector<double>& vPonderations, int iLeft, int iRight) 
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
