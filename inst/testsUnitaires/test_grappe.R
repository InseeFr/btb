# cf inst/testsunitaires/makeCluster.ods

library(RUnit)
library(dplyr)
library(DescTools)
Rcpp::sourceCpp('btb/src/constituerGrappes.cpp')
Rcpp::sourceCpp('btb/src/constituerMatriceEffectifs.cpp')

test.constituerGrappes <- function()
{
  dfObservations <- data.frame(x = c(15, 35, 15, 25, 35, 55, 45, 45, 55, 65, 70, 75, 85, 90, 65, 75, 85, 65, 70, 75, 85, 90, 65, 70, 75)
                            ,  y = c(10, 10, 30, 30, 35, 35, 45, 55, 55, 65, 65, 65, 65, 65, 70, 70, 70, 75, 75, 75, 75, 75, 85, 85, 85)
  )

  cellSize <- 20L

  # calcul de l'indice des observations - on prend le rectangle englobant et on positionne le debut de la numérotation sur la première observation
  dfObservations$col <- as.integer(floor((dfObservations$x) / cellSize) - floor(min(dfObservations$x / cellSize)) + 1)
  dfObservations$row <- as.integer(floor((dfObservations$y) / cellSize) - floor(min(dfObservations$y / cellSize)) + 1)
  
  #### matrice des effectifs
  iLongueur <- 2 ^ ceiling(log( max(dfObservations$col, dfObservations$row) ) / log(2))
  dfEffectifs <- dfObservations %>% group_by(col, row) %>% summarise(nbObs = n())
  mEffectifsSummarise <- DescTools::as.matrix.xtabs(xtabs(nbObs~row+col, data=dfEffectifs))

  mEffectifs <- constituerMatriceEffectifs(dfObservations$row - 1, dfObservations$col - 1)

  # vérifier qu'on obtient la même matrice des effectifs par deux méthodes différentes
  checkEquals(mEffectifs[1:5, 1:5], matrix(mEffectifsSummarise[1:5, 1:5], nrow = 5))
  
  # vérifier qu'on retrouve bien le même nombre d'observations 
  checkEquals(nrow(dfObservations), sum(mEffectifs))
  
  #### matrice des grappes
  mGrappes <- constituerGrappes(1, mEffectifs)

  mResultatAttendu = matrix(data = c(0,1,4,4,16,16,16,16,2,3,4,4,16,16,16,16,8,8,12,13,16,16,16,16,8,8,14,15,16,16,16,16,32,32,32,32,48,48,48,48,32,32,32,32,48,48,48,48,32,32,32,32,48,48,48,48,32,32,32,32,48,48,48,48), nrow = 8, ncol = 8)
  mResultatAttendu = as.integer(mResultatAttendu)
  mResultatAttendu = t(matrix(mResultatAttendu, ncol = 8))
  
  checkEquals(mResultatAttendu, mGrappes)
}

