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
#' @param pts_centro (**df object**) : table of centroids.  
#' @param sEPSG (**integer or character**) : epsg code.
#' @param iCellSize (**integer**) : cells size(s) in meters. Can be a vector for irregular grids
#' @param names_centro (**character vector**) : vector of names for longitude/latitude variables. Default c("x_centroide","y_centroide").
#'
#' @return (**df**) pts_centro table with additional *idInspire* variable
#' @export
#' @examples
#' pts_centro <- data.frame(x_centroide = c(100, 100, 300, 300, 500), 
#' y_centroide = c(100, 300, 100, 300, 100))
#' btb_add_inspire(pts_centro, sEPSG = 2154, iCellSize = 200)

btb_add_inspire <- function(pts_centro, sEPSG, iCellSize, names_centro = c("x_centroide","y_centroide")){
  
  # Checks : 
  stopifnot("Problem with centroid df" = nrow(pts_centro)>0) 
  stopifnot("Problem with centroids coordinates names" = sum(names_centro %in% colnames(pts_centro)) == 2) 
  stopifnot("sEPSG invalid" = nchar(as.character(sEPSG)) >= 4L)
  stopifnot("iCellSize invalid" = length(iCellSize) %in% c(1,nrow(pts_centro))) 
  
  nb_lines_checked <- min(15,nrow(pts_centro))
  res_modulo_x <- pts_centro[1:nb_lines_checked,names_centro[1]] %% iCellSize
  res_modulo_y <- pts_centro[1:nb_lines_checked,names_centro[2]] %% iCellSize
  stopifnot("Mismatch beetween centroids coordinates and iCellSize"=identical(length(unique(res_modulo_x)),1L) ) 
  stopifnot("Mismatch beetween centroids coordinates and iCellSize"=identical(length(unique(res_modulo_y)),1L) ) 
  
  # Code : 
  pts_centro$idInspire <- paste0(
    "CRS",sEPSG,
    "RES",iCellSize,"m",
    "N",format(pts_centro[[names_centro[2]]] - iCellSize / 2, scientific = F, trim = T),
    "E",format(pts_centro[[names_centro[1]]] - iCellSize / 2, scientific = F, trim = T)
  )
  return(pts_centro)
}

