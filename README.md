# Beyond the Border <img src="man/figures/logo.png" width=200 align="right" />

<!-- badges: start -->
  [![R-CMD-check](https://github.com/InseeFr/btb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/InseeFr/btb/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

[![CRAN_Status](http://www.r-pkg.org/badges/version/btb)](https://cran.r-project.org/package=btb)

`btb` ("Beyond the Border - Kernel Density Estimation for Urban Geography") is an R package which provides functions dedicated to urban analysis and density estimation using the KDE (kernel density estimator) method. 

A partial transposition of the package in Python is also available: [btbpy](https://github.com/InseeFrLab/btbpy).

## Description


The `btb_smooth()` function allows you to square and smooth geolocated data. It calculates a classical kernel smoothing (conservative) or a geographically weighted median. There are four major call modes of the function. 

- The first call mode is `btb_smooth(obs, epsg, cellsize, bandwidth)` for a classical kernel smoothing and automatic grid.
- The second call mode is `btb_smooth(obs, epsg, cellsize, bandwidth, quantiles)` for a geographically weighted median and automatic grid.
- The third call mode is `btb_smooth(obs, epsg, cellsize, bandwidth, centroids)` for a classical kernel smoothing and user grid.
- The fourth call mode is `btb_smooth(obs, epsg, cellsize, bandwidth, quantiles, centroids)` for a geographically weighted median and user grid.

## Installation

`btb` is available on CRAN and can therefore be readily installed
```
install.packages("btb")
```

To get a bug fix or to use a feature from the development version, you can install the development version of from GitHub :

```
install.packages("devtools")
devtools::install_github("InseeFr/btb")
```

## Usage 

Details on how to use the package can be found in its [documentation](man). Some applications for spatial smoothing are presented in [chapter 8](https://www.insee.fr/en/statistiques/fichier/3635545/imet131-l-chapitre-8.pdf) of the [Handbook of Spatial Analysis](https://www.insee.fr/en/information/3635545) published by Insee. You advise you to start by consulting the [vignette](https://inseefr.github.io/btb/articles/btb.html) of the package

## Contributions

Maintainer: Solène Colin <solene.colin@insee.fr>

Creators, authors and contributors: 

- Arlindo DOS SANTOS [aut],
- François SEMECURBE [aut],
- Julien PRAMIL [aut]
- Solène COLIN [cre, ctb],
- Kim ANTUNEZ [ctb],
- Auriane RENAUD [ctb],
- Farida MAROUCHI [ctb]
- Joachim TIMOTEO [ctb]


## References

- Geographically weighted summary statistics : a framework for localised exploratory data analysis, C.Brunsdon & al., in Computers, Environment and Urban Systems C.Brunsdon & al. (2002) [doi:10.1016/S0198-9715(01)00009-6](https://doi.org/10.1016/S0198-9715(01)00009-6) 
- Statistical Analysis of Spatial and Spatio-Temporal Point Patterns, Third Edition, Diggle, pp. 83-86, (2003) [doi:10.1080/13658816.2014.937718](https://doi.org/10.1080/13658816.2014.937718).
