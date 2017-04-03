# library(sp)
####################### Classe Grid #################################################
#
# date          version         auteur                          commentaire
# 2016/08/09      0.0.5      Francois Semecurbe     
# 2016/08/11      0.0.7      Francois Semecurbe     première version déployée sur le CRAN
# 2016/__/__      0.0.8      Arlindo Dos Santos     
# 2016/08/31      0.0.9      Francois Semecurbe     version diffusée à Audric Sophie et Cacheux Lionel
# 2016/10/03      0.1.0      Arlindo Dos Santos     version avec la fonction kernelSmoothingMedian
# 2016/10/04      0.1.1      Arlindo Dos Santos     correction de bugs mineurs
# 2016/10/14      0.1.2      Arlindo Dos Santos     4 modes d'appels pour kernelSmoothing
# 2016/12/01      0.1.4      Arlindo Dos Santos     fonction retour updateProgress pour kernelSmoothing
# 2017/01/01      0.1.6      Arlindo Dos Santos     prise en compte du paramètre neighbor si indiqué (et 0 si userGrid fournie)
# 2017/01/11      0.1.7      Arlindo Dos Santos     fonction retour updateProgress pour smoothingToGrid
# 2017/01/11      0.1.8      Arlindo Dos Santos     ajout de la colonne nbObsMin pour le lissage classique
#
##########################################################################################

#' @useDynLib btb
#' @importFrom Rcpp evalCpp
#' @import methods sp

#Grid est le nom de la classe definie
#les slots sont des variables typées
#' @export
setClass(
  Class = "Grid",
  slots = list(cellSize = "integer", bandwidth = "integer"),
  contains = "data.frame"
)

#' @export
setMethod(
  `[`,
  signature = signature(x = "Grid"),
  function(x, ...){
    # Save the original
    storedtdt <- x
    # Use the fact that x is a subclass to "data.frame"
    Nargs <- nargs()
    hasdrop <- "drop" %in% names(sys.call())
    if (Nargs == 2) {
      tmpDF <- `[.data.frame`(x, i = TRUE, j = i, ..., drop = FALSE)
    } else if ((Nargs == 3 && hasdrop)) {
      tmpDF <- `[.data.frame`(x, i = TRUE, j = i, ..., drop)
    } else if (hasdrop) {
      tmpDF <- `[.data.frame`(x, i, j, ..., drop)
    } else {
      tmpDF <- `[.data.frame`(x, i, j, ...)
    }
    # Reintegrate the results
    if (inherits(x = tmpDF, what = "data.frame"))
    {
      for (sName in names(getSlots("data.frame")))
      {
        slot(storedtdt, sName) <- slot(tmpDF, sName)
      }
      return(storedtdt)
    } else {
      return(tmpDF)
    }
  })

#' @export
setMethod(
  `[<-`,
  signature = signature(x = "Grid"),
  function(x, ..., value){
    # Save the original
    storedtdt <- x
    # Use the fact that x is a subclass to "data.frame"
    Nargs <- nargs()
    if (any(!names(sys.call()) %in% c("", "i", "j", "value"))) {
      stop("extra arguments are not allowed")
    }
    tmpDF <- data.frame(x)
    if (Nargs == 3) {
      if (missing(i)) 
        i <- j
      tmpDF[i] <- value
    } else if (Nargs == 4) {
      tmpDF[i, j] <- value
    }
    # Reintegrate the results
    for (sName in names(getSlots("data.frame")))
    {
      slot(storedtdt, sName) <- slot(tmpDF, sName)
    }
    return(storedtdt)
  })

############################## smoothingToGrid() #######################################################################
# fonction pour transformer un lissage en un fond de carte
#' @export
smoothingToGrid <- function(grid, epsg, fUpdateProgress = NULL)
{

  if (is.null(fUpdateProgress))
  {
    # version optimisee sans compte-rendu d'avancement du traitement
    carreauSpeed <- function(x, y, cellSize) 
    {
      sp::Polygons(list(sp::Polygon(
        cbind(
          c(0, cellSize, cellSize, 0, 0) + x - cellSize / 2,
          c(0, 0, cellSize, cellSize, 0) + y - cellSize / 2
        )
      )
      ), paste(x, y, sep = "_")
      )
    }
    
    grille <- mapply(carreauSpeed, grid[, "x"], grid[, "y"], MoreArgs = list(grid@cellSize))
  }
  else
  {
    carreauSlow <- function(i, grid, n, fUpdateProgress, startTime) {
      x <- grid[i, "x"]
      y <- grid[i, "y"]
      
      if (i == n)
      {# traitement du dernier polygone
        dTempsPasse = Sys.time() - startTime
        message <- paste0("\rGrid progress: ", ceiling(i / n * 100), "% - Elapsed time: ", floor(as.numeric(dTempsPasse, units = "mins") / 60), "m ", floor(as.numeric(dTempsPasse, units = "secs")) %% 60, "s   ")
        cat(message)
        fUpdateProgress(100, message)
      }
      else if (i %% floor( n / 100) == 0)
      {
        iPourcentageEffectue <- ceiling(i / n * 100)
        dTempsPasse = Sys.time() - startTime
        dTempsTotal <- dTempsPasse * 100 / iPourcentageEffectue;
        iTempsRestant <- ceiling(dTempsTotal - dTempsPasse);
        message <- paste0("\rGrid progress: ", iPourcentageEffectue, "% - remaining time: ", floor(as.numeric(iTempsRestant, units = "mins") / 60), "m ", as.numeric(iTempsRestant, units = "secs") %% 60, "s   ")
        cat(message)
        fUpdateProgress(iPourcentageEffectue, message)
      }
      
      sp::Polygons(list(sp::Polygon(
            cbind(
              c(0, grid@cellSize, grid@cellSize, 0, 0) + x - grid@cellSize / 2,
              c(0, 0, grid@cellSize, grid@cellSize, 0) + y - grid@cellSize / 2
            )
          )
        ), paste(x, y, sep = "_")
      )
    }
    
    cat("\n")
    dateDebutTraitement <- Sys.time()
    grille <- lapply(1:nrow(grid), carreauSlow, grid, length(grid@.Data[[1]]), fUpdateProgress, dateDebutTraitement)
    cat("\n")
  }

  # on crée un spatial polygon avec un code epsg de projection defini
  grille_spat <- sp::SpatialPolygons((grille), proj4string = sp::CRS(paste0("+init=epsg:", epsg)))
  df <- data.frame(grid@.Data)
  names(df) <- names(grid)
  data <- data.frame(ID = paste(df[, "x"], df[, "y"], sep = "_"), df)

  # un SpatialPolygonsDataFrame est un SpatialPolygon auquel on attache une table d'attributs
  return(sp::SpatialPolygonsDataFrame(grille_spat, data, match.ID = "ID"))
}

############################## kernelSmoothing() #######################################################################
# constructeur public
# 
# arguments
# dfObservations  : data.frame comportant les coordonnées géographiques (x,y), ainsi que les variables que l'on souhaite lisser
# cellSize        : Taille des carreaux
# bandwidth       : Rayon de lissage
# neighbor        : Paramètre technique pour calculer l'étendue des points d'estimation
# vQuantiles      : vecteur de quantiles à utiliser pour le lissage median
# dfCentroids     : data.frame avec deux colonnes, nommées x et y avec les coordonnees des centroides à utiliser
# fUpdateProgress : fonction permettant d'offrir à l'appelant une estimation de l'avancement du traitement 
# neighbor        : Paramètre technique pour calculer l'étendue des points d'estimations, à ne pas remplir
# 
# retourne
# un objet Grid dont le slot @data contient la valeur des variables lissees
# 
#' @export
kernelSmoothing <-
  function(dfObservations
          , cellSize
          , bandwidth
          , vQuantiles = NULL
          , dfCentroids = NULL
          , fUpdateProgress = NULL
          , neighbor = NULL
  )
  {
    cellSize <- as.integer(cellSize)
    bandwidth <- as.integer(bandwidth)
    dRayonMinimum <- cellSize * sqrt(2) / 2
    
    if(is.null(neighbor))
    {
      if (is.null(dfCentroids)) 
        neighbor <- max(0, ceiling(bandwidth / cellSize / 2L) - 1L) 
      else 
        neighbor <- 0
    }
    
    if (bandwidth < dRayonMinimum)
    {
      stop("bandwidth must be greater than cellSize * sqrt(2) / 2")
    }
    
    if (cellSize <= 0)
    {
      stop("cellSize must be greater than 0")
    }
    
    # les coordonnées des observations ne doivent pas être NA
    if (length(dfObservations[is.na(dfObservations$x), ]$x) + length(dfObservations[is.na(dfObservations$y), ]$y) > 0)
    {
      stop("NA coordinates are not allowed")
    } 
    
    # vérifier la regularite des centroides fournis par l'utilisateur
    if (!is.null(dfCentroids))
    {
      xOffset <- (dfCentroids$x + cellSize / 2) %% cellSize
      yOffset <- (dfCentroids$y + cellSize / 2) %% cellSize
      if (!all(xOffset == xOffset[1]) | !all(yOffset == yOffset[1]) )
      {
        stop("Centroids are not regular")
      }
    }else
    {
      xOffset <- 0
      yOffset <- 0
    }

    if (anyNA(dfObservations) )
    {
      warning("Be careful! NA values detected in your observations")
    }
    
    if (is.null(dfCentroids))
    { 
      # calcul de l'indice des observations - on prend le rectangle englobant et on positionne le debut de la numérotation sur la première observation
      dfObservations$i <- as.integer(floor((dfObservations$x - xOffset[1]) / cellSize) - floor(min(dfObservations$x / cellSize)) + 1)
      dfObservations$j <- as.integer(floor((dfObservations$y - yOffset[1]) / cellSize) - floor(min(dfObservations$y / cellSize)) + 1)
      
      # calcul des centroides
      dfCentroids <- data.frame( x = as.integer(floor(dfObservations$x / cellSize) * cellSize + (cellSize / 2)),
                                 y = as.integer(floor(dfObservations$y / cellSize) * cellSize + (cellSize / 2))
      )

      # les observations sont positionnées sur une matrice. mIndices[i, j] == 1 indique qu'il y a au moins 1 observation pour le carreau (i,j)
      mIndices <- matrix(0L, max(dfObservations$i), max(dfObservations$j))
      mIndices[cbind(dfObservations$i, dfObservations$j)] <- 1L
      
      # construction d'une matrice des indices des centroides étendue au voisinage
      mIndicesEtendus <- matrix(0L, nrow(mIndices) + 2 * neighbor, ncol(mIndices) + 2 * neighbor)
      
      # décalage de la matrice d'index sur la matrice étendue, ce qui permet de compter combien de fois un carreau est nécessaire
      for (voisin_x in -neighbor:neighbor)
      {
        for (voisin_y in -neighbor:neighbor)
        {
          mIndicesEtendus[neighbor + 1:nrow(mIndices) + voisin_x, neighbor + 1:ncol(mIndices) + voisin_y] <- 
            mIndicesEtendus[neighbor + 1:nrow(mIndices) + voisin_x, neighbor + 1:ncol(mIndices) + voisin_y] + mIndices
        }
      }
      
      # la matrice d'indices étendue est transformée en vecteurs de coordonnées
      vIndicesEtendus <- which(mIndicesEtendus > 0, arr.ind = TRUE)
      rm(list = c("mIndicesEtendus", "mIndices"))
      
      # retour aux coordonnées
      vX <- as.integer(round(min(dfCentroids$x) + (vIndicesEtendus[, 1] - 1 - neighbor) * cellSize))
      vY <- as.integer(round(min(dfCentroids$y) + (vIndicesEtendus[, 2] - 1 - neighbor) * cellSize))
      dfCentroidesUniques <- data.frame(x = vX, y = vY, i = vIndicesEtendus[, 1], j = vIndicesEtendus[, 2])
      
      rm(list = c("vX", "vY", "vIndicesEtendus", "dfCentroids"))
    }else
    {
      # Remarque: il n'est pas nécessaire de rapatrier les observations dans les carreaux de la grille fournie.
      # lors du lissage, les observations enverront leur contribution uniquement vers les carreaux fournis
      
      # calcul de l'indice des observations - on commence la numérotation pour la coordonnée minimale (qu'elle soit détenue par une observation ou un centroide)
      obsEtCentroides <- data.frame(x = c(dfObservations$x, dfCentroids$x), y = c(dfObservations$y, dfCentroids$y))
      indiceMinX <- floor(min(obsEtCentroides$x / cellSize))
      indiceMinY <- floor(min(obsEtCentroides$y / cellSize))
        
      dfObservations$i <- as.integer(floor((dfObservations$x - xOffset[1]) / cellSize) - indiceMinX + 1)
      dfObservations$j <- as.integer(floor((dfObservations$y - yOffset[1]) / cellSize) - indiceMinY + 1)
      
      # calcul de l'indice des centroides
      dfCentroids$i <- as.integer(floor(dfCentroids$x / cellSize) - indiceMinX + 1)
      dfCentroids$j <- as.integer(floor(dfCentroids$y / cellSize) - indiceMinY + 1)
      dfCentroidesUniques <- dfCentroids
    }
      
    nomColonnes <- colnames(dfObservations)
    listVar <- nomColonnes[nomColonnes != "x" & nomColonnes != "y" & nomColonnes != "i" & nomColonnes != "j"]
    
    if (is.null(vQuantiles))
    {
      # numérotation des centroides - décalage de -1 pour faire commencer la numerotation des lignes à 0 pour le traitement c++
      dfCentroidesUniques$index <- (1:nrow(dfCentroidesUniques)) - 1L
      
      # transformation en matrice
      mXcentroides = mIcentroides = mYcentroides = matrix(NA, max(dfCentroidesUniques$j), max(dfCentroidesUniques$i))
      mXcentroides[cbind(dfCentroidesUniques$j, dfCentroidesUniques$i)] <- dfCentroidesUniques$x
      mYcentroides[cbind(dfCentroidesUniques$j, dfCentroidesUniques$i)] <- dfCentroidesUniques$y
      mIcentroides[cbind(dfCentroidesUniques$j, dfCentroidesUniques$i)] <- dfCentroidesUniques$index
      
      iNbCentroidesUniques <- nrow(dfCentroidesUniques)
      dfResultat <- data.frame(dfCentroidesUniques[, c("x", "y")])
      
      mVariablesLissees <- rcppLissage(
          dfObservations$x
        , dfObservations$y
        , dfObservations$i
        , dfObservations$j
        , cellSize
        , bandwidth
        , neighbor
        , as.matrix(dfObservations[, listVar])
        , mXcentroides
        , mYcentroides
        , mIcentroides
        , iNbCentroidesUniques
        , fUpdateProgress
      )
      
      rm(list = c("dfObservations", "mXcentroides", "mYcentroides", "mIcentroides"))
      
      dfResultat <- cbind(dfResultat, mVariablesLissees)
      names(dfResultat) <- c("x","y", listVar)
      # names(dfResultat) <- c("x","y", listVar, "nbObsPondere") # version pour calcul de la colonne nbObsPondere
      
      rm(mVariablesLissees)
      
      return(new(Class = "Grid", dfResultat, cellSize = cellSize, bandwidth = bandwidth))    
    }else
    {
      mMedianes <- rcppLissageMedianSort(
          dfObservations$x
        , dfObservations$y
        , bandwidth
        , as.matrix(dfObservations[, listVar])
        , dfCentroidesUniques$x
        , dfCentroidesUniques$y
        , vQuantiles
        , fUpdateProgress
      )
    
      rm(dfObservations)
      
      dfResultat <- data.frame(cbind(dfCentroidesUniques$x, dfCentroidesUniques$y))
      dfResultat <- cbind(dfResultat, mMedianes)
      
      vNomQuantiles <- gsub("\\.", "", as.character(vQuantiles))
      
      dfListeVar <- expand.grid(vNomQuantiles, listVar)
      lNewVarNames <- paste0(dfListeVar[, 2], "_", dfListeVar[, 1])
      lNewVarNames <- c("x", "y", "nbObs", lNewVarNames)    
      colnames(dfResultat) <- lNewVarNames
      
      return(new(Class = "Grid", dfResultat, cellSize = cellSize, bandwidth = bandwidth))
    }
  }
