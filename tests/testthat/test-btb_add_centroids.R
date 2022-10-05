test_that("btb_add_centroids works", {
  
  # Test : works properly *************************************************
  
  pts1 <- data.frame( x = c(656913.1 , 348296.3 , 842276.3 , 716750.0 , 667418.2),
                     y = c(6855995 , 6788073 , 6385680 , 7003984 , 6585793), 
                     val=1:5)
  
  res1 <- btb::btb_add_centroids(pts1,100,names_centro = c("centroX","centroY"))
  classe_res1 <- class(res1)
  
  pts2 <- sf::st_as_sf(pts1,coords=c("x","y"),crs=2154)
  res2 <- btb_add_centroids(pts2,50)
  classe_res2 <- class(res2)
  
  testthat::expect_true(is.data.frame(res1))
  testthat::expect_true("sf" %in% classe_res2)
  
  # Test errors *************************************************
  
  testthat::expect_error(btb::btb_add_centroids(pts1)) # no iCellsize
  testthat::expect_error(btb::btb_add_centroids(pts1,iCellSize = -1)) # iCellsize <= 0
  testthat::expect_error(btb::btb_add_centroids(pts1,iCellSize = "carre")) # iCellsize not numeric
  
  # If pts is an sf object with degrees coordinates 
  pts3 <- sf::st_transform(pts2,4326) # GPS system
  testthat::expect_error(
    btb::btb_add_centroids(pts3,iCellSize = 75)
  )
  
  # Test warnings : duplicated column names
  testthat::expect_warning(
    btb::btb_add_centroids(pts1, iCellSize = 10 , names_centro =  c("x","y"),add = T)
  )

  testthat::expect_silent(
    btb::btb_add_centroids(pts1, iCellSize = 10 , names_centro = c("x","y"),add = F)
  )
  
  # Only works with points !
  squares <- btb::pixel_france[1:100,] %>% btb::btb_ptsToGrid(iCellSize = 200, names_centro = c('x','y'), sEPSG = 2154)
  testthat::expect_error(btb::btb_add_centroids(squares,iCellSize = 200))
  
  
  
  
})
