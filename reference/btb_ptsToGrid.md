# Compute a grid from centroid points

Function to compute a grid (regular or not) from centroid points.

(Fonction permettant de générer une grille (régulière ou non) à partir
de centroïdes)

## Usage

``` r
btb_ptsToGrid(
  pts,
  sEPSG = NA,
  iCellSize = NULL,
  names_centro = c("x_centro", "y_centro"),
  inspire = FALSE
)
```

## Arguments

- pts:

  A simple `data.frame` with the centroids coordinates of the squares to
  draw, or a `sf` object of centroides. To generate an irregular grid, a
  column with each cell size must be provided and named `iCellSize`.

  (Un simple `data.frame` comportant les coordonnées des carrés à
  dessiner, ou un objet `sf` des centroides. Pour obtenir une grille
  irrégulière, il faut fournir une colonne indiquant la taille de chaque
  carreau, et nommée `iCellSize`.

- sEPSG:

  EPSG code of projection (`character`). For example, the RGF93 /
  Lambert-93 projection has "2154" code.

  (code EPSG de la projection (`character`). Par exemple, la projection
  RGF93 / Lambert-93 a pour code "2154".)

- iCellSize:

  Cell size of the grid. If this argument is provided, the grid is
  regular.

  (Taille des carreaux de la grille. Si cet argument est fourni, la
  grille est régulière.)

- names_centro:

  (**character vector**)

  - vector of names for longitude/latitude variables. Default
    c("x_centro","y_centro").

  - vecteur des noms des variables de longitude/latitude. Par défaut :
    c("x_centro","y_centro")

- inspire:

  (boolean) : if TRUE, returns a column for Inspire grid names.

## Value

Returns an object of class `sf` and `data.frame`.

(Retourne un objet de classe `sf` et `data.frame`.)

## Examples

``` r
 
# example 1 - regular grid
pts <- data.frame(x_centro = c(100, 100, 300, 300, 500), 
y_centro = c(100, 300, 100, 300, 100))
carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
# write_sf(obj = carResult, dsn = "regularGrid.shp", delete_layer = TRUE)

# example 2 - irregular grid
pts <- data.frame(x = c(50, 50, 150, 150, 300)
                 , y = c(50, 150, 50, 150, 100)
                 , iCellSize = c(50, 50, 50, 50, 100))
carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154",names_centro=c("x","y"))
# write_sf(obj = carResult, dsn = "irregularGrid.shp", delete_layer = TRUE)
# Exemple 3 : sf points (no epsg)
pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
pts <- sf::st_as_sf(pts,coords=c("x","y"))
carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
# Exemple 3 : sf points (no epsg)
pts <- data.frame(x = c(100, 100, 300, 300, 500), 
y = c(100, 300, 100, 300, 100))
pts <- sf::st_as_sf(pts,coords=c("x","y"),crs=2154)
carResult <- btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)
```
