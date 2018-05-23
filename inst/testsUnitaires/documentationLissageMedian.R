
vModalites <- c(12000, 15000, 32000, 44000, 29000, 18000)
vPonderation <- c(0, 0.2344, 0, 0.6094, 0.8594, 0.9844)
vQuantiles <- c(0.1, 0.5, 0.9)

calculeQuantiles(vModalites = vModalites, vPonderation = vPonderation, vQuantiles = vQuantiles)




# revenus
dfObservations <- data.frame(x = c(15, 35, 15, 25, 35, 55, 45, 45, 55, 65, 70, 75, 85, 90, 65, 75, 85, 65, 70, 75, 85, 90, 65, 70, 75)
                             ,  y = c(10, 10, 30, 30, 35, 35, 45, 55, 55, 65, 65, 65, 65, 65, 70, 70, 70, 75, 75, 75, 75, 75, 85, 85, 85)
                             ,  valeur = c(10, 10, 30, 30, 35, 35, 45, 55, 55, 65, 65, 65, 65, 65, 70, 70, 70, 75, 75, 75, 75, 75, 85, 85, 85)
)

library(btb)
library(sf)
r <- cellSize / 2
df <- dfCentroidesUniques
df$geom <- sprintf("POLYGON ((%i %i, %i %i, %i %i, %i %i, %i %i))", df$x-r, df$y+r, df$x+r, df$y+r, df$x+r, df$y-r, df$x-r, df$y-r, df$x-r, df$y+r)
sfdf <- st_as_sf(df, wkt = "geom", crs = 2154)
rgdal::writeOGR(sfdf, paste0("D:/S3QCEA/temp/temp.shp"), "nomCouche", driver = "ESRI Shapefile", overwrite_layer = TRUE)
sf::write_sf(obj = sfdf, dsn = "D:/S3QCEA/temp/sample.shp")

xOffset = 0
yOffset = 0
cellSize <- 20
bandwidth <- 41
vQuantiles <- NULL
dfCentroids = NULL
fUpdateProgress = NULL
neighbor = max(0, ceiling(bandwidth / cellSize / 2L) - 1L) 
