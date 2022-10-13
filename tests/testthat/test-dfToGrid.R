test_that("dfToGrid works", {
  
  pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
  sEPSG=2154
  
  res_new <- btb::btb_ptsToGrid(pts = pts, iCellSize = 200, sEPSG=sEPSG,names_centro = c("x","y"))
  res_old <- btb::dfToGrid(df = pts, iCellSize = 200, sEPSG=sEPSG)
  
  
  
  testthat::expect_identical(res_new,res_old)
  
})