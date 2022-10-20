#' kernelSmoothing
#' @description Function kernelSmoothing has been replaced by `btb_smooth`
#' 
#' @param dfObservations points
#' @param sEPSG epsg code
#' @param iCellSize cell size
#' @param iBandwidth smoothing bandwidth
#' @param vQuantiles quantiles
#' @param dfCentroids user grid
#' @param fUpdateProgress message parameter
#' @param iNeighbor number of neighbors pixels
#' @param iNbObsMin technical parameter
#'
#' @return a warning message
#' @export

kernelSmoothing <- function(dfObservations,sEPSG, iCellSize, iBandwidth, vQuantiles = NULL, dfCentroids = NULL, fUpdateProgress = NULL, iNeighbor = NULL, iNbObsMin = 250){
  
  warning("`btb::kernelSmoothing()` was deprecated in btb 0.2.0. Please use `btb::btb_smooth())` instead.",
          immediate. = TRUE) 
  return(btb::btb_smooth(pts = dfObservations, 
                         sEPSG = sEPSG,
                         iCellSize = iCellSize ,
                         iBandwidth = iBandwidth,
                         vQuantiles = vQuantiles ,
                         dfCentroids = dfCentroids ,
                         iNeighbor = iNeighbor,
                         inspire = F,
                         iNbObsMin = iNbObsMin))
}
