# Beyond the Border - Kernel Density Estimation for Urban Geography 
[![CRAN_Status](http://www.r-pkg.org/badges/version/btb)](https://cran.r-project.org/package=btb)

`btb` ("Beyond the Border") is an R package which provides functions dedicated to urban analysis and density estimation using the KDE (kernel density estimator) method. 

A partial transposition of the package in Python is also available: [btbpy](https://github.com/InseeFrLab/btbpy).

## Description

 The kernelSmoothing() function allows you to square and smooth geolocated data. It calculates a classical kernel smoothing (conservative) or a geographically weighted median. There are four major call modes of the function. 

- The first call mode is kernelSmoothing(obs, epsg, cellsize, bandwith) for a classical kernel smoothing and automatic grid.
- The second call mode is kernelSmoothing(obs, epsg, cellsize, bandwith, quantiles) for a geographically weighted median and automatic grid.
- The third call mode is kernelSmoothing(obs, epsg, cellsize, bandwith, centroids) for a classical kernel smoothing and user grid.
- The fourth call mode is kernelSmoothing(obs, epsg, cellsize, bandwith, quantiles, centroids) for a geographically weighted median and user grid.

## Installation

`btb` is available on CRAN and can therefore be readily installed
```
install.packages("btb")
```

It is also possible to install it from the GitHub repository using `devtools`
```
install.packages("devtools")
devtools::install_github("InseeFrLab/btb")
```

## Usage 

Details on how to use the package can be found in its documentation. Some applications for spatial smoothing are presented in [chapter 8](https://www.insee.fr/en/statistiques/fichier/3635545/imet131-l-chapitre-8.pdf) of the [Handbook of Spatial Analysis](https://www.insee.fr/en/information/3635545) published by Insee.

## Contributions

Maintainer: Arlindo Dos Santos <Arlindo.Dos-Santos@insee.fr>

Authors and contributors:
- Arlindo Dos Santos [cre],
- François Sémécurbe [drt, aut],
- Auriane Renaud [ctb],
- Farida Marouchi [ctb]
- Joachim Timoteo [ctb]

## References

- Geographically weighted summary statistics : a framework for localised exploratory data analysis, C.Brunsdon & al., in Computers, Environment and Urban Systems C.Brunsdon & al. (2002) <doi:10.1016/S0198-9715(01)00009-6>, 
- Statistical Analysis of Spatial and Spatio-Temporal Point Patterns, Third Edition, Diggle, pp. 83-86, (2003) <doi:10.1080/13658816.2014.937718>.
