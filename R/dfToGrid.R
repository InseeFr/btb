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
  
  warning("`btb::dfToGrid()` was deprecated in btb 0.2.0. Please use `btb::btb_ptsToGrid())` instead.",
          immediate. = TRUE) 
  
  return(
    btb::btb_ptsToGrid(pts = df,
                       sEPSG = sEPSG,
                       iCellSize = iCellSize,
                       names_centro = c("x","y"),
                       inspire = FALSE)
    )
}
