library(RUnit)
Rcpp::sourceCpp('btb/src/rcppLissage.cpp')
Rcpp::sourceCpp('btb/src/rcppLissageMedianGrappe.cpp')
source("btb/R/kernelSmoothing.R")

test.kernelSmoothing <- function()
{
  # cf btb/tests/kernel_smoothing_2_obs.ods onglet obs1 & obs2
  lResultatAttendu <- list()
  lResultatAttendu[[1]] <- as.integer(c(10, 30, 50, 10, 30, 50, 10, 30, 50, 10, 30, 50))
  lResultatAttendu[[2]] <- as.integer(c(30, 30, 30, 50, 50, 50, 70, 70, 70, 90, 90, 90))
  lResultatAttendu[[3]] <- c(0.1978939, 1.1326320, 0.7356645, 2.2257595, 4.0865530, 2.3904194, 2.7336609, 4.3929076, 2.2948416, 1.0933719, 1.4505392, 0.2657565)
  
  iPas <- 20L
  iRayon <- 41L
  dfObservations <- data.frame(x = c(22, 35), y = c(70, 55), V1 = c(10, 13))
  
  # automatic grid
  # Call mode 1
  dtLisse <- kernelSmoothing(dfObservations, "2154", iPas, iRayon)
  checkEqualsNumeric(lResultatAttendu[[3]], dtLisse$V1, tolerance = 2e-8)
  checkEqualsNumeric(sum(dtLisse$V1), sum(dfObservations$V1))
  
  # Call mode 3
  vXCentroides <- rep(seq(from = 10, to = 50, by = iPas), 4)
  vYCentroides <- rep(seq(from = 30, to = 90, by = iPas), each = 3)
  dtCentroides <- data.frame(cbind(x = vXCentroides, y = vYCentroides))
  dtLisse <- kernelSmoothing(dfObservations, "2154", iPas, iRayon, dfCentroids = dtCentroides)
  checkEqualsNumeric(lResultatAttendu[[3]], dtLisse$V1, tolerance = 2e-8)
  checkEqualsNumeric(sum(dtLisse$V1), sum(dfObservations$V1))
  
  # user grid
  lResultatAttendu <- list()
  lResultatAttendu[[1]] = c(1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2)                         # nbObs
  lResultatAttendu[[2]] = c(13, 13, 13, 10, 10, 10, 10, 10, 10, 10, 10, 10)             # q 0.1
  lResultatAttendu[[3]] = c(13, 13, 13, 10, 13, 13, 10, 10, 13, 10, 10, 10)             # q 0.5
  lResultatAttendu[[4]] = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 10, 13, 13)             # q 0.9
  lResultatAttendu[[5]] = as.integer(c(10, 30, 50, 10, 30, 50, 10, 30, 50, 10, 30, 50)) # x
  lResultatAttendu[[6]] = as.integer(c(30, 30, 30, 50, 50, 50, 70, 70, 70, 90, 90, 90)) # y

  # Call mode 2
  dtLisse <- kernelSmoothing(dfObservations, "2154", iPas, iRayon, vQuantiles = c(0.1, 0.5, 0.9))
  checkEquals(lResultatAttendu[[1]], dtLisse$nbObs)
  checkEquals(lResultatAttendu[[2]], dtLisse$V1_01)
  checkEquals(lResultatAttendu[[3]], dtLisse$V1_05)
  checkEquals(lResultatAttendu[[4]], dtLisse$V1_09)
  checkEquals(lResultatAttendu[[5]], dtLisse$x)
  checkEquals(lResultatAttendu[[6]], dtLisse$y)
  
  # Call mode 4
  dtLisse <- kernelSmoothing(dfObservations, "2154", iPas, iRayon, dfCentroids = dtCentroides, vQuantiles = c(0.1, 0.5, 0.9))
  checkEquals(lResultatAttendu[[1]], dtLisse$nbObs)
  checkEquals(lResultatAttendu[[2]], dtLisse$V1_01)
  checkEquals(lResultatAttendu[[3]], dtLisse$V1_05)
  checkEquals(lResultatAttendu[[4]], dtLisse$V1_09)
  checkEquals(lResultatAttendu[[5]], dtLisse$x)
  checkEquals(lResultatAttendu[[6]], dtLisse$y)
  
  # tester les controles de validité
  dfObservations[1, "x"] <- NA
  checkException(kernelSmoothing(dfObservations, "2154", iPas, iRayon))
  
  # Vérifier que le lissage est conservatif
  load("btb/data/reunion.RData")
  reunionSmooth <- kernelSmoothing(dfObservations = reunion, "32740", iCellSize = 200L, iBandwidth = 400L)
  checkEquals(sum(reunion$houhold), sum(reunionSmooth$houhold))
}
