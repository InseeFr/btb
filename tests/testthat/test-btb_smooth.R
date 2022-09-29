test_that("btb_smooth works", {
  
  lResultatAttendu <- list()
  lResultatAttendu[[1]] <- as.integer(c(10, 30, 50, 10, 30, 50, 10, 30, 50, 10, 30, 50))
  lResultatAttendu[[2]] <- as.integer(c(30, 30, 30, 50, 50, 50, 70, 70, 70, 90, 90, 90))
  lResultatAttendu[[3]] <- c(0.19789385715986,1.13263202097131,0.73566453611993,2.22575946257372,4.08655298545109,2.39041940585168,2.73366093149227,4.39290763076285,2.29484163589690,1.09337186614331,1.45053918192735,0.26575648564970)
  
  iPas <- 20L
  iRayon <- 41L
  dfObservations <- data.frame(x = c(22, 35), y = c(70, 55), V1 = c(10, 13))
  
  # automatic grid
  # Call mode 1
  dtLisse <- btb::btb_smooth(dfObservations, "2154", iPas, iRayon)
  testthat::expect_equal(lResultatAttendu[[3]],dtLisse$V1, tolerance = 2e-8)
  testthat::expect_equal(sum(dtLisse$V1), sum(dfObservations$V1))
  
  # Call mode 3
  vXCentroides <- rep(seq(from = 10, to = 50, by = iPas), 4)
  vYCentroides <- rep(seq(from = 30, to = 90, by = iPas), each = 3)
  dtCentroides <- data.frame(cbind(x = vXCentroides, y = vYCentroides))
  dtLisse <- btb::btb_smooth(dfObservations, "2154", iPas, iRayon, dfCentroids = dtCentroides)
  testthat::expect_equal(lResultatAttendu[[3]], dtLisse$V1, tolerance = 2e-8)
  testthat::expect_equal(sum(dtLisse$V1), sum(dfObservations$V1))
  
  
  # user grid
  lResultatAttendu <- list()
  lResultatAttendu[[1]] = c(1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2)                         # nbObs
  lResultatAttendu[[2]] = c(13, 13, 13, 10, 10, 10, 10, 10, 10, 10, 10, 10)             # q 0.1
  lResultatAttendu[[3]] = c(13, 13, 13, 10, 13, 13, 10, 10, 13, 10, 10, 10)             # q 0.5
  lResultatAttendu[[4]] = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 10, 13, 13)             # q 0.9
  lResultatAttendu[[5]] = as.integer(c(10, 30, 50, 10, 30, 50, 10, 30, 50, 10, 30, 50)) # x
  lResultatAttendu[[6]] = as.integer(c(30, 30, 30, 50, 50, 50, 70, 70, 70, 90, 90, 90)) # y
  
  # Call mode 2
  dtLisse <- btb::btb_smooth(dfObservations, "2154", iPas, iRayon, vQuantiles = c(0.1, 0.5, 0.9))
  testthat::expect_equal(lResultatAttendu[[1]], dtLisse$nbObs)
  testthat::expect_equal(lResultatAttendu[[2]], dtLisse$V1_01)
  testthat::expect_equal(lResultatAttendu[[3]], dtLisse$V1_05)
  testthat::expect_equal(lResultatAttendu[[4]], dtLisse$V1_09)
  testthat::expect_equal(lResultatAttendu[[5]], dtLisse$x)
  testthat::expect_equal(lResultatAttendu[[6]], dtLisse$y)
  
  # Call mode 4
  dtLisse <- btb::btb_smooth(dfObservations, "2154", iPas, iRayon, dfCentroids = dtCentroides, vQuantiles = c(0.1, 0.5, 0.9))
  testthat::expect_equal(lResultatAttendu[[1]], dtLisse$nbObs)
  testthat::expect_equal(lResultatAttendu[[2]], dtLisse$V1_01)
  testthat::expect_equal(lResultatAttendu[[3]], dtLisse$V1_05)
  testthat::expect_equal(lResultatAttendu[[4]], dtLisse$V1_09)
  testthat::expect_equal(lResultatAttendu[[5]], dtLisse$x)
  testthat::expect_equal(lResultatAttendu[[6]], dtLisse$y)
  
  # tester les controles de validité
  dfObservations[1, "x"] <- NA
  testthat::expect_error(btb::btb_smooth(dfObservations, "2154", iPas, iRayon))
  
  # Vérifier que le lissage est conservatif
  #load("./data/reunion.RData")
  reunionSmooth <- btb::btb_smooth(dfObservations = btb::reunion, "3274", iCellSize = 200L, iBandwidth = 400L)
  testthat::expect_equal(sum(reunion$houhold), sum(reunionSmooth$houhold))
  
})


