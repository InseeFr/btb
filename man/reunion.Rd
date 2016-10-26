\name{reunion}
\alias{reunion}
\docType{package}
\title{
Households of Reunion (Menages de La Reunion) 
}

\description{
Gridded database with a grid cell resolution of 200 meters. Source : Insee, Localized data : Fiscal Revenue - 31/12/2010 and Housing Tax - 01/01/2011. 

(Donnees carroyees a 200 metres. Source : Insee, Revenus Fiscaux Localises (RFL) au 31 decembre 2010 et Taxe d'habitation (TH) au 1er janvier 2011)
}

\value{
  \item{x}{longitude}
  \item{y}{latitude}
  \item{houhold}{number of households - nombre de menages}
  \item{phouhold}{number of poor households - nombre de menages pauvres}
 }

 
\examples{
\dontrun{
data(reunion)
# Smoothing of the variable houhold (households) and phouhold (poor households) for Reunion
# Lissage de la variable houhold(menages) et phouhold(menages pauvres) pour la Reunion
reunionSmooth <- kernelSmoothing(reunion, 200, 800, c("houhold", "phouhold"))

# Calculating the poverty rate (Calcul du taux de pauvrete)
reunionSmooth$ratio <- reunionSmooth$phouhold / reunionSmooth$houhold * 100

# Building of the associated basemap (Creation du fond de carte associe)
grid <- smoothingToGrid(reunionSmooth, "32740")

#library(rgdal)
# Export of the basemap in shapefile format (Export du fond de carte au format shapefile)
#writeOGR(grid, "reunion.shp", "reunion", driver = "ESRI Shapefile")
}
}
