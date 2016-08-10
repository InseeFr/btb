\name{smoothing_to_grid}
\alias{smoothing_to_grid}
\docType{package}

\title{
Building of a basemap from smoothed data
(Fonction de transformation d un lissage en un fond de carte carroye)
}

\description{
Transforms smoothed data into grid data for export as cartographical data.
(Transforme un lissage en un carroyage pour export cartographique.)
}

\usage{
  smoothing_to_grid(df,epsg,cell_size=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{the \code{data.frame} obtained from the smoothing (le \code{data.frame} obtenu lors du lissage)}
  \item{epsg}{EPSG code (code EPSG de projection)}
  \item{cell_size}{leave empty (ne pas remplir normalement)}
 }
\details{
Returns a \code{spatialPolygonsdataframe}.
(Retourne un \code{spatialPolygonsdataframe}.)
}  

\author{ Psar Analyse Urbaine Insee - Thierry Cornely, Laure Genebes, Arlindo Dos Santos, Cynthia Faivre, Auriane Renaud and Francois Semecurbe }

\examples{
\dontrun{
data(reunion)
# Smoothing of the variable houhold (households) and phouhold (poor households) for Reunion
# Lissage de la variable houhold(menages) et phouhold(menages pauvres) pour La Reunion
reunion_smooth=kernel_smoothing(reunion,200,800,c("houhold","phouhold"))

# Calculating the poverty rate (Calcul du taux de pauvrete)
reunion_smooth$ratio=reunion_smooth$phouhold/reunion_smooth$houhold*100

# Building of the associated basemap (Creation du fond de carte associe)
grid=smoothing_to_grid(reunion_smooth,"32740")

#library(rgdal)
# Export of the basemap in shapefile format (Export du fond de carte au format shapefile)
#writeOGR(grid,"reunion.shp","reunion",driver="ESRI Shapefile")
}
}
