# Link points to their centroids

- Link some points to their centroids in a grid segmentation

- Relie des points aux centroides des carreaux auxquels ces points
  appartiennent (dans un découpage de l'espace en grille carroyée)

## Usage

``` r
btb_add_centroids(
  pts,
  iCellSize,
  offset = c(0L, 0L),
  names_coords = c("x", "y"),
  names_centro = c("x_centro", "y_centro"),
  add = TRUE
)
```

## Arguments

- pts:

  : points (`df` of `sf` object)

- iCellSize:

  :

  - Size of the square cells (**meters**)

  - Taille des côtés des carreaux (**mètres**)

- offset:

  (**numeric vector of size 2**)

  - Offset for a grid non centered on the geographical referential
    origin

  - Décalage si utilisation d'une grille non centrée sur l'origine du
    référentiel géographique

- names_coords:

  :

  - Names of the latitude-longitude variables (**character vector**)

  - Noms des variables de latitude-longitude (**vecteur character**)

- names_centro:

  - Names of the latitude-longitude variables for added centroids
    (**character vector**)

    - Noms des variables de latitude-longitude pour les centroïdes
      ajoutés (**vecteur character**)

- add:

  (**boolean**)

  - If TRUE : returns pts + centroids coordinates

  - Si TRUE : retourne pts + les coordonnées des centroïdes

## Value

- `pts` table with additional centroids coordinates `x_centro` and
  `y_centro` (`df` of `sf` object)

- Table `pts` avec les coordonnées des centroïdes `x_centro` and
  `y_centro` (objet `df` of `sf` )

## Details

Works with sf points but only with coordinates in meters (and not
degrees !). Do not use sf points with GPS coordinates for example.

## Examples

``` r
pts <- data.frame(
x = c(656913.1 , 348296.3 , 842276.3 , 716750.0 , 667418.2),
y = c(6855995 , 6788073 , 6385680 , 7003984 , 6585793),
val = 1:5)
btb_add_centroids(pts, 100, names_centro = c("centroX", "centroY"))
#>          x       y val centroX centroY
#> 1 656913.1 6855995   1  656950 6855950
#> 2 348296.3 6788073   2  348250 6788050
#> 3 842276.3 6385680   3  842250 6385650
#> 4 716750.0 7003984   4  716750 7003950
#> 5 667418.2 6585793   5  667450 6585750
btb_add_centroids(pts, 100, offset = c(50, 50), names_centro = c("centroX", "centroY"))
#>          x       y val centroX centroY
#> 1 656913.1 6855995   1  656900 6856000
#> 2 348296.3 6788073   2  348300 6788100
#> 3 842276.3 6385680   3  842300 6385700
#> 4 716750.0 7003984   4  716800 7004000
#> 5 667418.2 6585793   5  667400 6585800
pts2 <- sf::st_as_sf(pts, coords = c("x","y"), crs = 2154)
btb_add_centroids(pts2, 50)
#> Simple feature collection with 5 features and 5 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 348296.3 ymin: 6385680 xmax: 842276.3 ymax: 7003984
#> Projected CRS: RGF93 v1 / Lambert-93
#>   val        X       Y x_centro y_centro                 geometry
#> 1   1 656913.1 6855995   656925  6855975 POINT (656913.1 6855995)
#> 2   2 348296.3 6788073   348275  6788075 POINT (348296.3 6788073)
#> 3   3 842276.3 6385680   842275  6385675 POINT (842276.3 6385680)
#> 4   4 716750.0 7003984   716775  7003975   POINT (716750 7003984)
#> 5   5 667418.2 6585793   667425  6585775 POINT (667418.2 6585793)
```
