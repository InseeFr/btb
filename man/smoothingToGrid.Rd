\name{smoothingToGrid}
\alias{smoothingToGrid}
\docType{package}

\title{
Building of a basemap from smoothed data
(Fonction de transformation d'un lissage en un fond de carte carroye)
}

\description{
Transforms smoothed data into grid data for export as cartographical data.
(Transforme un lissage en un carroyage pour export cartographique.)
}

\usage{
  smoothingToGrid(grid, epsg, cellSize = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid}{the \code{grid} obtained from the smoothing (le \code{grid} obtenu lors du lissage)}
  \item{epsg}{EPSG code (code EPSG de projection)}
  \item{cellSize}{leave empty (ne pas remplir normalement)}
 }
\details{
Returns a \code{spatialPolygonsdataframe}.
(Retourne un \code{spatialPolygonsdataframe}.)
}  

\author{ Psar Analyse Urbaine Insee - Thierry Cornely, Laure Genebes, Arlindo Dos Santos, Cynthia Faivre, Auriane Renaud and Francois Semecurbe }

\examples{
\dontrun{
data(reunion)
# Smoothing all variables for Reunion (Lissage de toutes les variables pour la Reunion)

# Call mode 1: classic smoothing - automatic grid
reunionSmooth <- kernelSmoothing(dfObservations = reunion, cellSize = 200
                                , bandwidth = 400)

# Call mode 2: median smoothing - automatic grid
reunionSmooth <- kernelSmoothing(dfObservations = reunion, cellSize = 200
                                , bandwidth = 400, vQuantiles = c(0.1, 0.5, 0.9))

# Call mode 3: classic smoothing - user grid
dfCentroidsUser <- data.frame(x = seq(from = 314508, to = 378708, by = 200)
                            , y = seq(from = 7634103, to = 7691103, by = 200))
reunionSmooth <- kernelSmoothing(dfObservations = reunion, cellSize = 200
                                , bandwidth = 400, dfCentroids = dfCentroidsUser)

# Call mode 4: median smoothing - user grid
reunionSmooth <- kernelSmoothing(dfObservations = reunion, cellSize = 200
                                , bandwidth = 400, vQuantiles = c(0.1, 0.5, 0.9)
                                , dfCentroids = dfCentroidsUser)

# Building of the associated basemap (Creation du fond de carte associe)
grid <- smoothingToGrid(reunionSmooth, "32740")

# preview (Apercu)
library(sp)
library(cartography)
cartography::choroLayer(grid, var = "houhold", nclass = 5
                      , method = "fisher-jenks", border = NA)

# Export of the basemap in shapefile format (Export du fond de carte au format shapefile)
library(rgdal)
rgdal::writeOGR(grid, "reunion.shp", "reunion", driver = "ESRI Shapefile")
}}