#' @title Beyond the Border - Kernel Density Estimation for Urban Geography
#' @name btb
#' @description 
#' 
#' The kernelSmoothing() function allows you to square and smooth geolocated data. It calculates a classical kernel smoothing (conservative) or a geographically weighted median. There are four major call modes of the function. 
#' The first call mode is kernelSmoothing(obs, epsg, cellsize, bandwith) for a classical kernel smoothing and automatic grid.
#' The second call mode is kernelSmoothing(obs, epsg, cellsize, bandwith, quantiles) for a geographically weighted median and automatic grid.
#' The third call mode is kernelSmoothing(obs, epsg, cellsize, bandwith, centroids) for a classical kernel smoothing and user grid.
#' The fourth call mode is kernelSmoothing(obs, epsg, cellsize, bandwith, quantiles, centroids) for a geographically weighted median and user grid.
#' Geographically weighted summary statistics : a framework for localised exploratory data analysis, C.Brunsdon & al., in Computers, Environment and Urban Systems C.Brunsdon & al. (2002) <doi:10.1016/S0198-9715(01)00009-6>, 
#' Statistical Analysis of Spatial and Spatio-Temporal Point Patterns, Third Edition, Diggle, pp. 83-86, (2003) <doi:10.1080/13658816.2014.937718>.
#' @import methods sf RcppParallel mapsf dplyr
#' @importFrom Rcpp evalCpp
#' @useDynLib btb
NULL