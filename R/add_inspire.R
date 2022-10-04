#' btb_add_inspire
#' 
#' @description 
#' 
#' Creates Inpire names for a grid defined with :
#'   
#'   - centroids coordinates
#'   - squares size(s)
#'   - Projection system 
#' 
#' (Ajoute les noms des identifiants en norme Inspire des carreaux définis à partir de :
#'   
#'   - des coordonnées de leurs centroides
#'   - la taille de leur côté en mètres
#'   - un système de projection)
#'
#' @param pts (**df object**) : table of centroids.  
#' @param sEPSG (**integer or character**) : epsg code.
#' @param iCellSize (**integer**) : cells size(s) in meters. Can be a vector for irregular grids
#' @param names_centro (**character vector**) : vector of names for longitude/latitude variables
#'
#' @return (**df**) pts table with additional *idInspire* variable
#' @export
#' @examples
#' pts <- data.frame(x_centroide = c(100, 100, 300, 300, 500), y_centroide = c(100, 300, 100, 300, 100))
#' btb_add_inspire(pts, sEPSG = 2154, iCellSize = 200)

btb_add_inspire <- function(pts, sEPSG, iCellSize, names_centro = c("x_centroide","y_centroide")){
  
  # Checks : 
  stopifnot("Problem with centroid df" = nrow(df)>0) 
  stopifnot("Problem with centroids coordinates names" = sum(names_centro %in% colnames(pts)) == 2) 
  stopifnot("sEPSG invalid" = identical(nchar(as.character(sEPSG)),4L)) 
  stopifnot("iCellSize invalid" = length(iCellSize) %in% c(1,nrow(pts))) 
  
  # Code : 
  pts$idInspire <- paste0(
    "CRS",sEPSG,
    "RES",iCellSize,"m",
    "N",format(pts[[names_centro[2]]] - iCellSize / 2, scientific = F, trim = T),
    "E",format(pts[[names_centro[1]]] - iCellSize / 2, scientific = F, trim = T)
  )
  return(pts)
}

