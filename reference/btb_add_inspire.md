# Creates Inpire names for a grid

Creates Inpire names for a grid defined with :

- centroids coordinates

- Squares size(s)

- Projection system

Ajoute les noms des identifiants en norme Inspire des carreaux définis à
partir :

- des coordonnées de leurs centroides

- de la taille de leur côté en mètres

- d'un système de projection)

## Usage

``` r
btb_add_inspire(
  pts_centro,
  sEPSG,
  iCellSize,
  names_centro = c("x_centro", "y_centro")
)
```

## Arguments

- pts_centro:

  (**df object**) : table of centroids.

- sEPSG:

  (**integer or character**) : epsg code.

- iCellSize:

  (**integer**) : cells size(s) in meters. Can be a vector for irregular
  grids

- names_centro:

  (**character vector**) : vector of names for longitude/latitude
  variables. Default c("x_centro","y_centro").

## Value

(**df**) pts_centro table with additional *idInspire* variable

## Examples

``` r
pts_centro <- data.frame(x_centro = c(100, 100, 300, 300, 500), 
y_centro = c(100, 300, 100, 300, 100))
btb_add_inspire(pts_centro, sEPSG = 2154, iCellSize = 200)
#>   x_centro y_centro              idInspire
#> 1      100      100     CRS2154RES200mN0E0
#> 2      100      300   CRS2154RES200mN200E0
#> 3      300      100   CRS2154RES200mN0E200
#> 4      300      300 CRS2154RES200mN200E200
#> 5      500      100   CRS2154RES200mN0E400
```
