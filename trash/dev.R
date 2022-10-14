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

pkgs_imports <- c("RcppParallel","mapsf")
lapply(pkgs_imports, function(pkg) {
  if (system.file(package = pkg) == '') install.packages(pkg)
})

library(Rcpp)
library(devtools)

compileAttributes()
document()

check()

load_all()
install()

# ************************
# Documentation
# ************************
usethis::use_roxygen_md()
roxygen2md::roxygen2md()
devtools::document()
help("constituerGrappes")
pkgload::dev_help('constituerGrappes')

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
usethis::use_test("btb_smooth")
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

# Vignette

usethis::use_vignette("vignette_btb")
build_rmd("vignettes/vignette_btb.Rmd")

# pkgdown website
usethis::use_pkgdown()
pkgdown::build_site()
