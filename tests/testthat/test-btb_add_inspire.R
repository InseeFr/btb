test_that("btb_add_inspire works !", {
  
  
  # Test normal behaviour
  pts_centro1 <- pixel_france[1:1000,] 
  
  testthat::expect_true(
    is.data.frame(
      btb::btb_add_inspire(pts_centro1, iCellSize = 1000, sEPSG = "2154",names_centro = c("x","y"))
    )
  )
  
  pts_centro2 <- data.frame(x_centroide=c(150,350,550,150,350,550),y_centroide=c(175,175,175,375,375,375))
  testthat::expect_true(
    is.data.frame(
      btb::btb_add_inspire(pts_centro2, iCellSize = 200, sEPSG = "2154")
    )
  )
  
  # Test error if points don't represent centroids
  pts1 <- btb::dfRestaurantParis
  testthat::expect_error(
    btb::btb_add_inspire(pts1, iCellSize = 200,sEPSG = "2154", names_centro = c("x","y"))
  )
  
  
})
