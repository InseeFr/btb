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
  
})


