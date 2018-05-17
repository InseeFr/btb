library(RUnit)
Rcpp::sourceCpp('btb/src/lissageMedian.cpp')

test.calculeQuantiles <- function()
{
  ######## elements tries - ponderations entieres - nb elements impair ########
  # test 1
  vQuantilesResult <- calculeQuantiles(c(2, 5, 6), c(1, 1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(2, 5, 6))
  
  # test 2
  vQuantilesResult <- calculeQuantiles(c(6, 7, 8, 9, 10), c(1, 1, 1, 1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(6, 8, 10)) 
  
  # test 3
  vQuantilesResult <- calculeQuantiles(c(1, 2, 3, 8, 9), c(8, 2, 1, 6, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(1, 2, 8))
  
  ######## elements tries - ponderations entieres - nb elements pair
  # test 4
  vQuantilesResult <- calculeQuantiles(c(8, 9), c(1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(8, 8.5, 9))
  
  # test 5
  vQuantilesResult <- calculeQuantiles(c(7, 8, 9, 10), c(1, 1, 1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8.5, 10))
  
  # test 6
  vQuantilesResult <- calculeQuantiles(c(7, 8, 9, 10), c(12, 3, 6, 8), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8, 10))
  
  ######## elements tries - ponderations decimales - nb elements impair
  # test 7
  vQuantilesResult <- calculeQuantiles(c(7, 8, 9, 10), c(12, 3, 6, 8), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8, 10)) 
  
  # test 8
  vQuantilesResult <- calculeQuantiles(c(6, 7, 8, 9, 10), c(0.3, 0.2, 0.1, 0.2, 0.3), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(6, 8, 10)) 
  
  ######## elements tries - ponderations d?cimales - nb elements pair
  # test 9
  vQuantilesResult <- calculeQuantiles(c(7, 8, 9, 10), c(0.5, 0.4, 0.3, 0.6), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8.5, 10)) 
  
  # test 10
  vQuantilesResult <- calculeQuantiles(c(7, 8, 9, 10), c(0.11, 0.1, 0.11, 0.09), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8, 10)) 
  
  # test 11
  vQuantilesResult <- calculeQuantiles(c(7, 8, 9, 10), c(0.11, 0.09, 0.11, 0.09), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8.5, 10)) 

    ######## elements non tries - ponderations entieres - nb elements impair
  # test 12
  vQuantilesResult <- calculeQuantiles(c(9, 8, 7), c(1, 1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8, 9)) 
  
  # test 13
  vQuantilesResult <- calculeQuantiles(c(10, 7, 6, 9, 8), c(1, 1, 1, 1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(6, 8, 10)) 
  
  ######## elements non tries - ponderations entieres - nb elements pair
  # test 14
  vQuantilesResult <- calculeQuantiles(c(9, 8), c(1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(8, 8.5, 9)) 
  
  # test 15
  vQuantilesResult <- calculeQuantiles(c(7, 10, 9, 8), c(1, 1, 1, 1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8.5, 10)) 
  
  # test 16
  vQuantilesResult <- calculeQuantiles(c(10, 8, 9, 7), c(8, 2, 6, 12), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8.5, 10)) 
  
  ######## elements non tries - ponderations d?cimales - nb elements impair
  # test 17
  vQuantilesResult <- calculeQuantiles(c(7, 8, 9), c(0.1, 0.1, 0.1), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8, 9)) 
  
  # test 18
  vQuantilesResult <- calculeQuantiles(c(10, 7, 8, 9, 6), c(0.3, 0.2, 0.1, 0.2, 0.3), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(6, 8, 10)) 
  
  ######## elements non tries - ponderations d?cimales - nb elements pair
  # test 19
  vQuantilesResult <- calculeQuantiles(c(7, 8, 10, 9), c(0.5, 0.4, 0.6, 0.3), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8.5, 10)) 
  
  # test 20
  vQuantilesResult <- calculeQuantiles(c(7, 8, 10, 9), c(0.11, 0.1, 0.09, 0.11), c(0.1, 0.5, 0.9))
  checkEquals(vQuantilesResult, c(7, 8, 10)) 
  
  ############# tests validite arguments ###############
  checkException(calculeQuantiles(c(), c(), c()))
  checkException(calculeQuantiles(c(), c(1), c(1)))
  checkException(calculeQuantiles(c(1), c(), c(1)))
  checkException(calculeQuantiles(c(1), c(1), c()))
  checkException(calculeQuantiles(c(1), c(1, 2), c(1)))
  checkException(calculeQuantiles(c(1), c(1), c(2)))
  checkException(calculeQuantiles(c(1), c(1), c(-0.5)))
}

test.rcppLissageMedian <- function()
{
vXObservations <- c(22, 35)
vYObservations <- c(70, 55)
iPas <- 20
iRayon <- 41
mVar <- matrix(c(10, 13, 15, 17), nrow = 2, ncol = 2, byrow = FALSE)
dimnames(mVar) = list( c("row1", "row2"),  c("V1", "V2"))
vXCentroides <- rep(seq(from = 10, to = 90, by = 20), 5)
vYCentroides <- rep(seq(from = 10, to = 90, by = 20), each = 5)
vQuantiles <- c(0.1, 0.5, 0.9)
  lissageMedian <- rcppLissageMedian(vXObservations, vYObservations, iRayon, mVar, vXCentroides, vYCentroides, vQuantiles )

  mResultatAttendu <- matrix(0, nrow = 25, ncol = 7)
  mResultatAttendu[ 6, ] <- c(1, 13, 13, 13, 17, 17, 17)
  mResultatAttendu[ 7, ] <- c(2, 13, 13, 13, 17, 17, 17)
  mResultatAttendu[ 8, ] <- c(1, 13, 13, 13, 17, 17, 17)
  mResultatAttendu[11, ] <- c(2, 10, 10, 13, 15, 15, 17)
  mResultatAttendu[12, ] <- c(2, 10, 13, 13, 15, 17, 17)
  mResultatAttendu[13, ] <- c(2, 10, 13, 13, 15, 17, 17)
  mResultatAttendu[14, ] <- c(1, 13, 13, 13, 17, 17, 17)
  mResultatAttendu[16, ] <- c(2, 10, 10, 13, 15, 15, 17)
  mResultatAttendu[17, ] <- c(2, 10, 10, 13, 15, 15, 17)
  mResultatAttendu[18, ] <- c(2, 10, 13, 13, 15, 17, 17)
  mResultatAttendu[19, ] <- c(1, 13, 13, 13, 17, 17, 17)
  mResultatAttendu[21, ] <- c(1, 10, 10, 10, 15, 15, 15)
  mResultatAttendu[22, ] <- c(2, 10, 10, 13, 15, 15, 17)
  mResultatAttendu[23, ] <- c(2, 10, 10, 13, 15, 15, 17)

  checkEquals(lissageMedian[, 1:7], mResultatAttendu)
}
