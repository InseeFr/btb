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


devtools::check()
devtools::document()


# ************************
# Documentation
# ************************
usethis::use_roxygen_md()
roxygen2md::roxygen2md()
devtools::document()
help("btb_ptsToGrid")

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
# Rcpp::Rcpp.package.skeleton("btbis")
# Rcpp::compileAttributes()


