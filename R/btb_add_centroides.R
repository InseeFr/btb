
#' btb_add_centroides
#' 
#' @description 
#' 
#' - Link some points to their centroids in a grid segmentation
#' - Relie des points aux centroides des carreaux auxquels ces points appartiennent (dans un découpage de l'espace en grille carroyée)
#'
#' @param pts : points (`df` of `sf` object)
#' @param iCellSize : 
#' - size of the square cells (**meters**)
#' - Taille des côtés des carreaux (**mètres**)
#' @param coords :
#' - Names of the latitude-longitude variables (**character vector**)
#' - Noms des variables de latitude-longitude (**vecteur character**) 
#'
#' @return 
#' - `pts` table with additional centroids coordinates `x_centro` and `y_centro` (`df` of `sf` object)
#' - Table `pts` avec les coordonnées des centroïdes `x_centro` and `y_centro` (objet `df` of `sf` )
#' @export
#'
#' @details 
#'  Works with sf points but only with coordinates in meters (and not degrees !). Do not use sf points with GPS coordinates for example.
#' @examples
#' pts <- data.frame(
#' x = c(656913.1 , 348296.3 , 842276.3 , 716750.0 , 667418.2),
#' y = c(6855995 , 6788073 , 6385680 , 7003984 , 6585793),
#' val=1:5)
#' btb_add_centroides(pts,100)
#' pts2 <- sf::st_as_sf(pts,coords=c("x","y"),crs=2154)
#' btb_add_centroides(pts2,50)



btb_add_centroides <- function(pts,iCellSize,coords=c("x","y")){
  
  # Checks *************
  
  stopifnot("Incorrect cells size"=is.numeric(iCellSize))
  
  if("sf" %in% class(pts)){
    proj_units <- sf::st_crs(pts , parameters = TRUE)$units_gdal
    stopifnot("Coordintates unit must be meters (not degrees)"=identical(proj_units,"metre"))
  }else{
    stopifnot("Coordinates names not found"=coords %in% colnames(pts))
    stopifnot("NA values in coordinates"=identical(sum(is.na(pts[,coords])),0L))
  }
  if(sum(c("x_centro","y_centro") %in% colnames(pts))>0){
    warning("Variables x_centro and y_centro have been replaced !")
  }
  
  # code ****************
  
  if("sf" %in% class(pts)){
    coord_numeric <- sf::st_coordinates(pts)
    pts <- pts %>% cbind(coord_numeric)
    coords <- colnames(coord_numeric)
  }
  
  x_centro <-  pts[[coords[1]]] - (pts[[coords[1]]] %% iCellSize) + (iCellSize / 2)
  y_centro <-  pts[[coords[2]]] - (pts[[coords[2]]] %% iCellSize) + (iCellSize / 2)
  pts <- pts %>% cbind(x_centro=x_centro,y_centro=y_centro)
  return(pts)
}

