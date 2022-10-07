#' @title  btb_ptsToGrid
#' 
#' @encoding UTF-8
#' 
#' @description 
#' Function to compute a grid (regular or not) from centroid points. 
#' 
#' (Fonction permettant de générer une grille (régulière ou non) à partir de centroïdes)
#' 
#' @param pts 
#' A simple `data.frame` with the centroids coordinates of the squares to draw, or a `sf` object of centroides. To generate an irregular grid, a column with each cell size must be provided and named `iCellSize`.
#' 
#' (Un simple \code{data.frame} comportant les coordonnées des carrés à dessiner, ou un objet `sf` des centroides. Pour obtenir une grille irrégulière, il faut fournir une colonne indiquant la taille de chaque carreau, et nommée `iCellSize`.
#' 
#' @param sEPSG 
#' EPSG code of projection (\code{character}). For example, the RGF93 / Lambert-93 projection has "2154" code.
#' 
#' (code EPSG de la projection (\code{character}). Par exemple, la projection RGF93 / Lambert-93 a pour code "2154".)
#' 
#' @param iCellSize 
#' 
#' Cell size of the grid. If this argument is provided, the grid is regular.
#' 
#' (Taille des carreaux de la grille. Si cet argument est fourni, la grille est régulière.)
#' @param names_centro  (**character vector**)
#'  - vector of names for longitude/latitude variables. Default c("x_centroide","y_centroide").
#'  - vecteur des noms des variables de longitude/latitude. Par défaut :  c("x_centroide","y_centroide")
#' @param inspire (boolean) : if TRUE, returns a column for Inspire grid names.
#'
#' @return 
#' Returns an object of class \code{sf} and \code{data.frame}. 
#' 
#' (Retourne un objet de classe \code{sf} et \code{data.frame}.)
#' 
#' @examples  
#' # example 1 - regular grid
#' pts <- data.frame(x_centroide = c(100, 100, 300, 300, 500), 
#' y_centroide = c(100, 300, 100, 300, 100))
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
#' # write_sf(obj = carResult, dsn = "regularGrid.shp", delete_layer = TRUE)
#' 
#' # example 2 - irregular grid
#' pts <- data.frame(x = c(50, 50, 150, 150, 300)
#'                  , y = c(50, 150, 50, 150, 100)
#'                  , iCellSize = c(50, 50, 50, 50, 100))
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154",names_centro=c("x","y"))
#' # write_sf(obj = carResult, dsn = "irregularGrid.shp", delete_layer = TRUE)
#' # Exemple 3 : sf points (no epsg)
#' pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
#' pts <- sf::st_as_sf(pts,coords=c("x","y"))
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
#' # Exemple 3 : sf points (no epsg)
#' pts <- data.frame(x = c(100, 100, 300, 300, 500), 
#' y = c(100, 300, 100, 300, 100))
#' pts <- sf::st_as_sf(pts,coords=c("x","y"),crs=2154)
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
#' @export






btb_ptsToGrid <- function(pts, sEPSG=NA, iCellSize = NULL, names_centro = c("x_centroide","y_centroide"), inspire = F)
{
  # Test of parameters
  stopifnot("pts must be a df object"= is.data.frame(pts))
  stopifnot("No size for cells"= !is.null(iCellSize) | !is.null(pts[["iCellSize"]]))
  stopifnot("sEPSG not valid "= is.na(sEPSG) | identical(nchar(as.character(sEPSG)),4L))
  
  # If pts is a sf objet
  if("sf" %in% class(pts)){
    stopifnot("sf pts : geometry type must be POINT"= identical(as.vector(sf::st_geometry_type(pts,by_geometry=F)),"POINT") )
    
    # Retrieve epsg code from sf points (or NA)
    if(is.na(sEPSG)) sEPSG <- (sf::st_crs(pts))$epsg
    
    # As df with coordinates
    pts[[names_centro[1]]] <- sf::st_coordinates(pts)[,1]
    pts[[names_centro[2]]] <- sf::st_coordinates(pts)[,2]
    pts <- sf::st_drop_geometry(pts) 
  }
  
  stopifnot("Centroids coordinates names not found"=names_centro %in% colnames(pts))
  
  # Withdraw iCellSize as vector for iregular grid
  if(is.null(iCellSize)) iCellSize <- pts$iCellSize
  
  
  # Add Inpire id
  
  if(inspire){
    pts <- btb::btb_add_inspire(pts, sEPSG, iCellSize, names_centro = names_centro)
  }

  
  if (!is.null(iCellSize)) # if regular grid
  {
    r = iCellSize / 2
    pts$geometry <-
      sprintf(
        "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
        pts[[names_centro[1]]] - r,
        pts[[names_centro[2]]] + r,
        pts[[names_centro[1]]] + r,
        pts[[names_centro[2]]] + r,
        pts[[names_centro[1]]] + r,
        pts[[names_centro[2]]] - r,
        pts[[names_centro[1]]] - r,
        pts[[names_centro[2]]] - r,
        pts[[names_centro[1]]] - r,
        pts[[names_centro[2]]] + r
      )
  }
  else # If irregular grid
    pts$geometry <-
      sprintf(
        "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
        pts[[names_centro[1]]] - pts$iCellSize / 2,
        pts[[names_centro[2]]] + pts$iCellSize / 2,
        pts[[names_centro[1]]] + pts$iCellSize / 2,
        pts[[names_centro[2]]] + pts$iCellSize / 2,
        pts[[names_centro[1]]] + pts$iCellSize / 2,
        pts[[names_centro[2]]] - pts$iCellSize / 2,
        pts[[names_centro[1]]] - pts$iCellSize / 2,
        pts[[names_centro[2]]] - pts$iCellSize / 2,
        pts[[names_centro[1]]] - pts$iCellSize / 2,
        pts[[names_centro[2]]] + pts$iCellSize / 2
      )
  
  sfpts <- sf::st_as_sf(pts, wkt = "geometry", crs = as.integer(sEPSG))
  
  return(sfpts)
}
