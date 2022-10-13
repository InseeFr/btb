#' dfToGrid
#' 
#' @description Function dfToGrid has been replaced by btb_ptsToGrid
#' @param df centroid points
#' @param sEPSG epsg code
#' @param iCellSize cell size
#'
#' @return a grid
#' @export
#'

dfToGrid <- function(df, sEPSG, iCellSize = NULL){
  warning("Function dfToGrid is deprecated : please use btb_ptsToGrid instead")
  return(
    btb::btb_ptsToGrid(pts = df,
                       sEPSG = sEPSG,
                       iCellSize = iCellSize,
                       names_centro = c("x","y"),
                       inspire = F)
    )
}
