# Depreciated function of btb_smooth

Function kernelSmoothing has been replaced by `btb_smooth`

## Usage

``` r
kernelSmoothing(
  dfObservations,
  sEPSG,
  iCellSize,
  iBandwidth,
  vQuantiles = NULL,
  dfCentroids = NULL,
  fUpdateProgress = NULL,
  iNeighbor = NULL,
  iNbObsMin = 250
)
```

## Arguments

- dfObservations:

  points

- sEPSG:

  epsg code

- iCellSize:

  cell size

- iBandwidth:

  smoothing bandwidth

- vQuantiles:

  quantiles

- dfCentroids:

  user grid

- fUpdateProgress:

  message parameter

- iNeighbor:

  number of neighbors pixels

- iNbObsMin:

  technical parameter

## Value

a warning message
