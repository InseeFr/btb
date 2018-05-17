#ifndef CONSTITUER_GRAPPE_H
#define CONSTITUER_GRAPPE_H

void quadTree2( const unsigned int iNbObsMin
                   , const unsigned short profondeurMax
                   , const arma::Mat<int>& mEffectifs2n
                   , arma::Mat<int>& mGrappe
                   , std::vector<int> &vNoGrappe
                   , const unsigned short profondeurCourante
                   , const unsigned int iTaille
                   , const unsigned int iRowReference
                   , const unsigned int iColReference
);
  
arma::Mat<int> constituerGrappes2(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs, std::vector<int> &vNoGrappe);
arma::Mat<int> constituerGrappes2(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs);

void quadTree( const unsigned int iNbObsMin
                   , const unsigned short profondeurMax
                   , const arma::Mat<int>& mEffectifs2n
                   , arma::Mat<int>& mGrappe
                   , std::vector<int> &vNoGrappe
                   , const unsigned short profondeurCourante
                   , const unsigned int iTaille
                   , const unsigned int iRowReference
                   , const unsigned int iColReference
);

arma::Mat<int> constituerGrappes(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs, std::vector<int> &vNoGrappe);
arma::Mat<int> constituerGrappes(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs);

void decomposer(const int iNiveauMax, const int iNiveaucourant, const int iNombre, std::vector<int> &vIndice);
std::vector<int> coordonneesGrappe(int iNiveauMax, int iNoGrappe);
  
void test(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs, std::vector<int>& vNoGrappe);
void test(const unsigned int iNbObsMin, const arma::Mat<int>& mEffectifs);
  
#endif
