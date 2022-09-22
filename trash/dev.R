req_pkgs <- c(
  'available', 
  'desc', 
  'usethis', 
  'gitlabr',
  'git2r',
  'devtools',
  'roxygen2',
  'roxygen2md',
  'testthat',
  'covr'
)

lapply(req_pkgs, function(pkg) {
  if (system.file(package = pkg) == '') install.packages(pkg)
})

pkgs_imports <- c("RcppParallel")
lapply(pkgs_imports, function(pkg) {
  if (system.file(package = pkg) == '') install.packages(pkg)
})

devtools::check()
devtools::document()


# ************************
# Documentation
# ************************
usethis::use_roxygen_md()
roxygen2md::roxygen2md()
devtools::document()
help("btb_ptsToGrid")
pkgload::dev_help('btb_smooth')

# ************************
# Développement
# ************************


usethis::use_package("sf") # IMporter des fonctions de packages externes
# Au sein d'une fonction
#' @importFrom purrr map 
usethis::use_pipe(export = TRUE) # si besoin

# ************************
# Tests
# ************************
usethis::use_testthat()
usethis::use_test("btb_ptsToGrid")
devtools::test()
# covr::package_coverage()




# Charger toutes les fonctions du package
devtools::load_all()

# Installer le package
devtools::check()
devtools::install()
devtools::build()


# Volet RCPP
Rcpp::Rcpp.package.skeleton("btbis")
Rcpp::compileAttributes()


# 
# 
# iCellSize <- 20L
# iBandwidth <- 41L
# dfObservations <- data.frame(x = c(22, 35), y = c(70, 55), V1 = c(10, 13))
# vQuantiles = NULL
# dfCentroids = NULL
# fUpdateProgress = NULL
# iNeighbor = NULL
# iNbObsMin = 250

btb::btb_smooth(dfObservations = data.frame(x = c(22, 35), y = c(70, 55), V1 = c(10, 13)),
                iCellSize = 20L,iBandwidth = 41L 
                )

pts <- data.frame(x = c(100, 100, 300, 300, 500), y = c(100, 300, 100, 300, 100))
btb::btb_ptsToGrid(pts = pts, sEPSG = "2154", iCellSize = 200)


# vXCentroides <- rep(seq(from = 10, to = 50, by = iCellSize), 4)
# vYCentroides <- rep(seq(from = 30, to = 90, by = iCellSize), each = 3)
# dfCentroids <- data.frame(cbind(x = vXCentroides, y = vYCentroides))

#library(btb)
# data(dfPrix_SP95_2016)
# dfPrix_SP95_2016$nbObs <- 1L
# dfSmoothed <- btb::btb_smooth(dfObservations = dfPrix_SP95_2016, 
#                          sEPSG = "2154", 
#                          iCellSize = 5000L, 
#                          iBandwidth = 30000L)
