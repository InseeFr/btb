#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix rcpp_lissage(NumericVector x_estim,NumericVector y_estim,NumericVector x_pt,NumericVector y_pt,NumericVector rayon,NumericMatrix var, NumericVector taille_int) 
  {
    int col=var.ncol();
    int bandwith2=pow(rayon[0],2); /*la bande passante*/
    int taille_intermediaire=taille_int[0];
    int estim_length=x_estim.length(); /*la taille du vecteur des points d estimation*/
    int pt_length=x_pt.length();
    int i,j,index,t;
    float x_pt_temp,y_pt_temp,dist_carree;
    NumericVector sommeponderation(pt_length, 0.0);
    NumericMatrix liss(estim_length,col);
    NumericVector ponderation(taille_intermediaire,0.0);
    NumericVector temp(taille_intermediaire,0.0); 
    NumericVector coord_estimation_ponderation(taille_intermediaire,0.0);
   
    for (i=0;i<pt_length;i++) /*on parcourt le vecteur des points*/
     {
      x_pt_temp=x_pt[i]; 
      y_pt_temp=y_pt[i];
      //ponderation=temp;
      index=0;
       for (j=0;j<estim_length;j++) /*on parcourt le vecteur des points d estimation*/
        {
          /*calcul de la distance entre le point courant et le point d estimation*/
          dist_carree=pow(x_pt_temp-x_estim[j],2)+pow(y_pt_temp-y_estim[j],2);
          /*on affecte un poids uniquement aux points situes a l interieur de la bande passante*/
          /*le poids correspond a la fonction bisquare*/
          if (dist_carree<bandwith2) 
          {
            coord_estimation_ponderation[index]=j;
            ponderation[index]=pow(1-(dist_carree/bandwith2),2);
            sommeponderation[i]+=ponderation[index]; 
            index+=1;
          }
        }
      if (sommeponderation[i]>0)
        {
          for (j=0;j<index;j++)
          {
            for (t=0;t<col;t++) /*pour chacune des variables a lisser*/
              {
                /*on calcule le lissage : ponderation de la valeur de la variable par le poids de lissage*/
               liss(coord_estimation_ponderation[j],t)+=var(i,t)*ponderation[j]/sommeponderation[i];
              }
          
        }
      }
   }
  return(liss);
}  
