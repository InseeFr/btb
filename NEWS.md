# btb 0.2.0

* Change of maintainer, new authors (november 2022)
* Change `dfToGrid` name (to `btb_ptsToGrid`)
* Allow `btb_ptsToGrid` to use simple feature points as input
* Change `kernelSmoothing` name (to `btb_smooth`)
* Allow `btb_ptsToGrid` to use simple feature points as input
* New function to link points to their centroid in a grid (`btb_add_centroides`)
* Updates messages when smoothing removed (in RcppLissage.cpp)
* Old functions `dfToGrid` and `kernelSmoothing` replaced by warnings
* Possibility to use Inspire naming for grid cells (with `btb_add_inspire` function) 
* Roxygen documentation
* testthat compatible tests
* New tests
* Add vignette
* Add pkgdown website

<!---
# btb 0.1.8
* add `iNbObsMin` for classical smoothing
* `smoothingToGrid` integrated in `kernelSmoothing`

# btb 0.1.7

* function `updateProgress` for `smoothingToGrid`

# btb 0.1.6

* `iNeighbor` parameter
-->

# btb 0.1.30.3 
* New CRAN version (last version submitted by Arlindo)

# btb 0.1.19

* improvement: ckecking if NA exists in coordinates and variables with a loops (faster than anyNA)
* improvement: two matrix suppressed in parameters sended to rcppLissage
* new function : constituerGrappes
* new function : constituerMatriceEffectifs

# btb 0.1.18

* improvement: kernelSmoothing uses armadillo matrix instead of NumericMatrix

# btb 0.1.17

* improvement: kernelSmoothing with quantiles splits computation in parallel clusters (with RcppParallel) 

# 0.1.16

* improvement: smoothingToGrid uses sf package instead of sp package

# btb 0.1.15

* improvement: verbose output to console only if fUpdateProgress is not provided. (if provided, the output must be managed by the calling function)
* improvement: use of left-value [, "var"] instead of $var

# btb 0.1.13

* bug fixed: upper boundary of loop for was underestimated when bandwidth > neighbor o cellSize
* bug fixed: the offset was wrong when dfCendroids was provided
* improving window when looking for cells to send smoothed value
* new argument for kernelSmoothing and smoothingToGrid: fUpdateProgress
* neighbor argument is 0 if dfCendroids is provided
* adding this NEWS file
* adding unit tests
* encoding in UTF-8 (description file)

# btb 0.1.3

* change c++ library call for solaris compatibility 

# btb 0.1.2

* bug fixed: manage missing values
* bug fixed: doubles comparison
* rename bandwith in bandwidth
* improved documentation
* improved memory improved (double to int when possible)
* kernelSmoothing accept a new argument: dfCentroids
* new smoothing mode: mobile quantile
* log when na value is found
* check if bandwidth > sqrt(cellSize) / 2
* check if x and y are not null
* encoding source files in UTF-8  

# btb 0.1.1

* small bugs correction

# btb 0.1.0

* `kernelSmoothingMedian` function

# btb 0.0.7

* initial release of the package




