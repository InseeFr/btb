\name{dfToGrid}
\alias{dfToGrid}
\docType{package}

\title{Grid function (Fonction de carroyage)}

\description{
Function to compute a grid (regular or not) from a data.frame. 

(Fonction permettant de générer une grille (régulière ou non) à partir d'un data.frame.)
}

\usage{dfToGrid(df, sEPSG, iCellSize = NULL)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{A \code{data.frame} with the centroids coordinates of the squares to draw. To generate an irregular grid, a third column wiht each cell size must be provided. (x, y, iCellSize)
  
  (Un \code{data.frame} comportant les coordonnées des carrés à dessiner. Pour obtenir une grille irrégulière, il faut fournir une troisième colonne indiquant la taille de chaque carreau. (x, y, iCellSize).)}
  
  \item{sEPSG}{EPSG code of projection (\code{character}). For example, the RGF93 / Lambert-93 projection has "2154" code.
  
  (code EPSG de la projection (\code{character}). Par exemple, la projection RGF93 / Lambert-93 a pour code "2154".)}
  
  \item{iCellSize}{Cell size of the grid. If this argument is provided, the grid is regular.
  
  (Taille des carreaux de la grille. Si cet argument est fourni, la grille est régulière.)}
  
}

\value{
Returns an object of class \code{sf} and \code{data.frame}. 

(Retourne un objet de classe \code{sf} et \code{data.frame}.)
} 
 
\author{Psar Analyse Urbaine Insee - Thierry Cornely, Laure Genebes, Arlindo Dos Santos, Cynthia Faivre, Auriane Renaud and Francois Semecurbe
}

\examples{
library(sf)
# example 1 - regular grid
df <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
dfResult <- dfToGrid(df = df, sEPSG = "2154", iCellSize = 200)
# write_sf(obj = dfResult, dsn = "regularGrid.shp", delete_layer = TRUE)

# example 2 - irregular grid
df <- data.frame(x = c(50, 50, 150, 150, 300)
               , y = c(50, 150, 50, 150, 100)
               , iCellSize = c(50, 50, 50, 50, 100))
dfResult <- dfToGrid(df = df, sEPSG = "2154")
# write_sf(obj = dfResult, dsn = "irregularGrid.shp", delete_layer = TRUE)
}