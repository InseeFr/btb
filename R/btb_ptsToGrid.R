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
#'
#' @return 
#' Returns an object of class \code{sf} and \code{data.frame}. 
#' 
#' (Retourne un objet de classe \code{sf} et \code{data.frame}.)
#' 
#' @examples  
#' # example 1 - regular grid
#' pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
#' # write_sf(obj = carResult, dsn = "regularGrid.shp", delete_layer = TRUE)
#' 
#' # example 2 - irregular grid
#' pts <- data.frame(x = c(50, 50, 150, 150, 300)
#'                  , y = c(50, 150, 50, 150, 100)
#'                  , iCellSize = c(50, 50, 50, 50, 100))
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154")
#' # write_sf(obj = carResult, dsn = "irregularGrid.shp", delete_layer = TRUE)
#' # Exemple 3 : sf points (no epsg)
#' pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
#' pts <- st_as_sf(pts,coords=c("x","y"))
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
#' #' Exemple 3 : sf points (no epsg)
#' pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
#' pts <- st_as_sf(pts,coords=c("x","y"),crs=2154)
#' carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
#' @export






btb_ptsToGrid <- function(pts, sEPSG=NA, iCellSize = NULL, inspire = F)
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
    pts$x <- sf::st_coordinates(pts)[,1]
    pts$y <- sf::st_coordinates(pts)[,2]
    pts <- sf::st_drop_geometry(pts) 
  }
  
  # Withdraw iCellSize as vector for iregular grid
  if(is.null(iCellSize)) iCellSize <- pts$iCellSize
  
  
  # Add Inpire id
  
  if(inspire){
    pts <- add_inspire(pts, sEPSG, iCellSize)
  }

  
  if (!is.null(iCellSize))
  {
    r = iCellSize / 2
    pts$geometry <-
      sprintf(
        "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
        pts$x - r,
        pts$y + r,
        pts$x + r,
        pts$y + r,
        pts$x + r,
        pts$y - r,
        pts$x - r,
        pts$y - r,
        pts$x - r,
        pts$y + r
      )
  }
  else
    pts$geometry <-
      sprintf(
        "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
        pts$x - pts$iCellSize / 2,
        pts$y + pts$iCellSize / 2,
        pts$x + pts$iCellSize / 2,
        pts$y + pts$iCellSize / 2,
        pts$x + pts$iCellSize / 2,
        pts$y - pts$iCellSize / 2,
        pts$x - pts$iCellSize / 2,
        pts$y - pts$iCellSize / 2,
        pts$x - pts$iCellSize / 2,
        pts$y + pts$iCellSize / 2
      )
  
  sfpts <- sf::st_as_sf(pts, wkt = "geometry", crs = as.integer(sEPSG))
  
  
  
  
  return(sfpts)
}
