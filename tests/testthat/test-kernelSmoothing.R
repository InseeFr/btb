test_that("kernelSmoothing works", {
  
  iPas <- 20L
  iRayon <- 41L
  pts <- data.frame(x = c(22, 35), y = c(70, 55), V1 = c(10, 13))
  
  res_new <- btb::btb_smooth(pts, "2154", iPas, iRayon)
  res_old <- btb::kernelSmoothing(pts, "2154", iPas, iRayon)
  
  testthat::expect_identical(res_new,res_old)
  
})
