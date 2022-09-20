test_that("btb_ptsToGrid works", {
  library(sf)
  
  pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
  dfResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
  
  # Le résultat est un objet de classe sf
  testthat::expect_equal(class(dfResult),c("sf","data.frame"))
})
