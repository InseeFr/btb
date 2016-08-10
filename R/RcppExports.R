lissage <- function(x_estim, y_estim, x_pt, y_pt, rayon,var,taille_int) {
    .Call('btb_rcpp_lissage', PACKAGE = 'btb', x_estim, y_estim, x_pt, y_pt, rayon,var,taille_int)
}

