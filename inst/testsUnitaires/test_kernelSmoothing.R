library(RUnit)
Rcpp::sourceCpp('btb/src/rcppLissage.cpp')

test.kernelSmoothing <- function()
{
  # cf btb/tests/kernel_smoothing_2_obs.ods onglet obs1 & obs2
  lResultatAttendu <- list()
  lResultatAttendu[[1]] <- as.integer(c(10, 30, 50, 10, 30, 50, 10, 30, 50, 10, 30, 50))
  lResultatAttendu[[2]] <- as.integer(c(30, 30, 30, 50, 50, 50, 70, 70, 70, 90, 90, 90))
  lResultatAttendu[[3]] <- c(0.1978939, 1.1326320, 0.7356645, 2.2257595, 4.0865530, 2.3904194, 2.7336609, 4.3929076, 2.2948416, 1.0933719, 1.4505392, 0.2657565)
  
  iPas <- 20
  iRayon <- 41
  dfObservations <- data.frame(x = c(22, 35), y = c(70, 55), V1 = c(10, 13))
  
  # automatic grid
  # Call mode 1
  dfLisse <- kernelSmoothing(dfObservations, iPas, iRayon)
  checkEqualsNumeric(lResultatAttendu[[3]], dfLisse@.Data[[3]], tolerance = 2e-8)
  
  # Call mode 3
  vXCentroides <- rep(seq(from = 10, to = 50, by = iPas), 4)
  vYCentroides <- rep(seq(from = 30, to = 90, by = iPas), each = 3)
  dfCentroides <- data.frame(cbind(x = vXCentroides, y = vYCentroides))
  dfLisse <- kernelSmoothing(dfObservations, iPas, iRayon, dfCentroids = dfCentroides)
  dfLisse <- dfLisse[dfLisse$V1 > 0, ]
  checkEqualsNumeric(lResultatAttendu[[3]], dfLisse@.Data[[3]], tolerance = 2e-8)
  
  # user grid
  lResultatAttendu <- list()
  lResultatAttendu[[1]] = as.integer(c(10, 30, 50, 10, 30, 50, 10, 30, 50, 10, 30, 50))
  lResultatAttendu[[2]] = as.integer(c(30, 30, 30, 50, 50, 50, 70, 70, 70, 90, 90, 90))
  lResultatAttendu[[3]] = c(1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2)
  lResultatAttendu[[4]] = c(13, 13, 13, 10, 10, 10, 10, 10, 10, 10, 10, 10)
  lResultatAttendu[[5]] = c(13, 13, 13, 10, 13, 13, 10, 10, 13, 10, 10, 10)
  lResultatAttendu[[6]] = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 10, 13, 13)
  # lResultatAttendu[[7]] = c(17, 17, 17, 15, 15, 15, 15, 15, 15, 15, 15, 15)   # V2
  # lResultatAttendu[[8]] = c(17, 17, 17, 15, 17, 17, 15, 15, 17, 15, 15, 15)   # V2
  # lResultatAttendu[[9]] = c(17, 17, 17, 17, 17, 17, 17, 17, 17, 15, 17, 17)   # V2
  
  # Call mode 2
  dfLisse <- kernelSmoothing(dfObservations, iPas, iRayon, vQuantiles = c(0.1, 0.5, 0.9))
  checkEquals(lResultatAttendu, dfLisse@.Data)
  
  # Call mode 4
  dfLisse <- kernelSmoothing(dfObservations, iPas, iRayon, dfCentroids = dfCentroides, vQuantiles = c(0.1, 0.5, 0.9))
  checkEquals(lResultatAttendu, dfLisse@.Data)
  
}