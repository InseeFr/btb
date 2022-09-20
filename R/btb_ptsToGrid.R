#' @title  btb_ptsToGrid
#' 
#' @encoding UTF-8
#' 
#' @description 
#' Function to compute a grid (regular or not) from a data.frame. 
#' 
#' (Fonction permettant de générer une grille (régulière ou non) à partir d'un data.frame.)
#' 
#' @param pts 
#' A `data.frame` with the centroids coordinates of the squares to draw. To generate an irregular grid, a third column wiht each cell size must be provided. (x, y, iCellSize)
#' 
#' (Un \code{data.frame} comportant les coordonnées des carrés à dessiner. Pour obtenir une grille irrégulière, il faut fournir une troisième colonne indiquant la taille de chaque carreau. (x, y, iCellSize).)
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
#' library(sf)
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
#' @export


btb_ptsToGrid <- function(pts, sEPSG, iCellSize = NULL)
{
  # Test of parameters
  stopifnot("pts must be a df object"= is.data.frame(pts))
  stopifnot("No size for cells"= !is.null(iCellSize) | !is.null(pts[["iCellSize"]]))
  
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
  
  sfpts <- st_as_sf(pts, wkt = "geometry", crs = as.integer(sEPSG))
  return(sfpts)
}