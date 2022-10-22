#' @title Beyond the Border - Kernel Density Estimation for Urban Geography
#' @name btb
#' @description 
#' 
#' The kernelSmoothing() function allows you to square and smooth geolocated data. It calculates a classical kernel smoothing (conservative) or a geographically weighted median. There are four major call modes of the function. 
#' The first call mode is kernelSmoothing(obs, epsg, cellsize, bandwidth) for a classical kernel smoothing and automatic grid.
#' The second call mode is kernelSmoothing(obs, epsg, cellsize, bandwidth, quantiles) for a geographically weighted median and automatic grid.
#' The third call mode is kernelSmoothing(obs, epsg, cellsize, bandwidth, centroids) for a classical kernel smoothing and user grid.
#' The fourth call mode is kernelSmoothing(obs, epsg, cellsize, bandwidth, quantiles, centroids) for a geographically weighted median and user grid.
#' @import methods sf RcppParallel mapsf dplyr
#' @importFrom Rcpp evalCpp
#' @useDynLib btb
NULL