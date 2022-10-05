test_that("btb_ptsToGrid works", {

  pts1 <- data.frame(x_centroide = c(100, 100, 300, 300, 500), y_centroide = c(100, 300, 100, 300, 100))
  sfpts1 <- sf::st_as_sf(pts1,coords=c("x_centroide","y_centroide"),crs=2154)
  
  result0 <- btb::btb_ptsToGrid(pts = pts1, iCellSize = 200)
  result1 <- btb_ptsToGrid(pts = pts1, sEPSG = "2154", iCellSize = 200)
  result2 <- btb_ptsToGrid(pts = sfpts1, iCellSize = 200)
  
  # The result is an sf object
  testthat::expect_equal(class(result1),c("sf","data.frame"))
  
  # Not working if no iCellSize
  testthat::expect_error(btb_ptsToGrid(pts = pts1, sEPSG = "2154"))
  
  # Same result with df and sf 
  testthat::expect_identical(result1,result2)
  
  # Take into account epsg of sf object
  testthat::expect_false(identical(result0,result2))
})
