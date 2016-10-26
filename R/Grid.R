# library(sp)
####################### Classe Grid #################################################
#
# date          version         auteur                          commentaire
# 2016/08/09      0.0.5      Francois Semecurbe     
# 2016/08/11      0.0.7      Francois Semecurbe     premiere version deployee sur le CRAN
# 2016/__/__      0.0.8      Arlindo Dos Santos     
# 2016/08/31      0.0.9      Francois Semecurbe     version diffusee a Audric Sophie et Cacheux Lionel
# 2016/10/03      0.1.0      Arlindo Dos Santos     version avec la fonction kernelSmoothingMedian
# 2016/10/04      0.1.1      Arlindo Dos Santos     correction de bugs mineurs
# 2016/10/14      0.1.2      Arlindo Dos Santos     4 modes d'appels pour kernelSmoothing
#
##########################################################################################

#' @useDynLib btb
#' @importFrom Rcpp evalCpp
#' @import methods sp

#Grid est le nom de la classe definie
#les slots sont des variables typees
#' @export
setClass(
  Class = "Grid",
  slots = list(cellSize = "integer", bandwidth = "integer"),
  contains = "data.frame"
)

#' @export
setMethod(
  `[`,
  signature=signature(x="Grid"),
  function(x, ...){
    # Save the original
    storedtdt <- x
    # Use the fact that x is a subclass to "data.frame"
    Nargs <- nargs()
    hasdrop <- "drop" %in% names(sys.call())
    if(Nargs==2) {
      tmpDF <- `[.data.frame`(x, i=TRUE, j=i, ..., drop=FALSE)
    } else if((Nargs==3 && hasdrop)) {
      tmpDF <- `[.data.frame`(x, i=TRUE, j=i, ..., drop)
    } else if(hasdrop) {
      tmpDF <- `[.data.frame`(x, i, j, ..., drop)
    } else {
      tmpDF <- `[.data.frame`(x, i, j, ...)
    }
    # Reintegrate the results
    if (inherits(x=tmpDF, what="data.frame")){
      for(sName in names(getSlots("data.frame"))){
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
  signature=signature(x="Grid"),
  function(x, ..., value){
    # Save the original
    storedtdt <- x
    # Use the fact that x is a subclass to "data.frame"
    Nargs <- nargs()
    if (any(!names(sys.call()) %in% c("", "i", "j", "value"))) {
      stop("extra arguments are not allowed")
    }
    tmpDF <- data.frame(x)
    if(Nargs==3) {
      if (missing(i)) i <- j
      tmpDF[i] <- value
    } else if(Nargs==4) {
      tmpDF[i, j] <- value
    }
    # Reintegrate the results
    for(sName in names(getSlots("data.frame"))){
      slot(storedtdt, sName) <- slot(tmpDF, sName)
    }
    return(storedtdt)
  })

#creation de fond de carte dont la projection a le code epsg defini
setGeneric(name = "creationFondDeCarte", def = function(object, epsg) {standardGeneric("creationFondDeCarte")})

#on definit plus specifiquement la fonction creationFondDeCarte pour un objet Grid
setMethod(f = "creationFondDeCarte", signature = "Grid", definition = function(object, epsg)
  {
    #fonction de creation de carreaux de cote cellSize
    carreau <-
      function(x, y, cellSize) 
      {
        (sp::Polygons(list(sp::Polygon(
          cbind(
            c(0, cellSize, cellSize, 0, 0) + (x - (cellSize) / 2),
            c(0, 0, cellSize, cellSize, 0) + (y - (cellSize) / 2)
          )
        )), paste(x, y, sep = "_")))
      }
    
    #on applique la fonction carreau
    #MoreArgs contient les parametres supplementaires de la fonction carreau
    grille <- mapply(carreau, object[, "x"], object[, "y"], MoreArgs = list(cellSize = object@cellSize))
    
    #on cree un spatial polygon avec un code epsg de projection defini
    grille_spat <- sp::SpatialPolygons((grille), proj4string = sp::CRS(paste("+init=epsg:", epsg, sep = "")))
    df <- data.frame(object@.Data)
    names(df) <- names(object)
    data <- data.frame(ID = paste(df[, "x"], df[, "y"], sep = "_"), df)
    
    #un SpatialPolygonsDataFrame est un SpatialPolygon auquel on attache une table d'attributs
    return(sp::SpatialPolygonsDataFrame(grille_spat, data, match.ID = "ID"))
  }
)

# fonction pour transformer un lissage en un fond de carte
#' @export
smoothingToGrid <- function(grid, epsg, cellSize = NULL)
{
  if (is.null(cellSize)) {
    cellSize <- grid@cellSize
  }
  bandwidth <- NULL
  xcol <- "x"
  ycol <- "y"
  if (is.null(bandwidth)) {
    bandwidth <- 0L
  }
  car <- grid[, c(xcol, ycol, (liste_var = names(grid)[!(names(grid) %in% c(xcol, ycol))]))]
  names(car) <- c("x", "y", liste_var)
  return(creationFondDeCarte(
    new(
      Class = "Grid",
      car,
      cellSize = cellSize,
      bandwidth = bandwidth
    ),
    epsg
  ))
}

############################## kernelSmoothing() #######################################################################################
# constructeur public
# 
# arguments
# dfObservations  : data.frame comportant les coordonnees geographiques (x,y), ainsi que les variables que l'on souhaite lisser
# cellSize        : Taille des carreaux
# bandwidth       : Rayon de lissage
# neighbor        : Parametre technique pour calculer l'etendue des points d'estimation
# vQuantiles      : vecteur de quantiles a utiliser pour le lissage median
# dfCentroids     : data.frame avec deux colonnes, nommees x et y avec les coordonnees des centroides a utiliser
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
           , neighbor = max(0, ceiling(bandwidth / cellSize / 2L) - 1L)
  )
  {
    cellSize <- as.integer(cellSize)
    bandwidth <- as.integer(bandwidth)
    dRayonMinimum <- cellSize * sqrt(2) / 2
    
    if (bandwidth < dRayonMinimum)
    {
      stop("bandwidth must be greater than cellSize * sqrt(2) / 2")
    }
    
    if (cellSize <= 0)
    {
      stop("cellSize must be greater than 0")
    }
    
    # les coordonnees des observations ne doivent pas etre NA
    if(length(dfObservations[is.na(dfObservations$x), ]$x) + length(dfObservations[is.na(dfObservations$y), ]$y) > 0)
    {
      stop("NA coordinates are not allowed")
    } 
    
    # verifier la regularite des centroides fournis par l'utilisateur
    if(!is.null(dfCentroids))
    {
      xOffset <- (dfCentroids$x + cellSize / 2) %% cellSize
      yOffset <- (dfCentroids$y + cellSize / 2) %% cellSize
      if (!all(xOffset == xOffset[1]) | !all(yOffset == yOffset[1]) )
      {
        stop("Centroids are not regular")
      }
    }
    else
    {
      xOffset <- 0
      yOffset <- 0
    }

    if (anyNA(dfObservations) )
    {
      warning("Be careful! NA values detected in your observations")
    }
    
    # calcul de l'indice des observations - on prend le rectangle englobant et on positionne le debut de la numerotation sur la premiere observation
    dfObservations$i <- as.integer(floor((dfObservations$x - xOffset[1]) / cellSize) - floor(min(dfObservations$x / cellSize)) + 1)
    dfObservations$j <- as.integer(floor((dfObservations$y - yOffset[1]) / cellSize) - floor(min(dfObservations$y / cellSize)) + 1)
    
    if (is.null(dfCentroids))
    { # calcul des centroides
      dfCentroids <- data.frame( x = as.integer(floor(dfObservations$x / cellSize) * cellSize + (cellSize / 2)),
                                 y = as.integer(floor(dfObservations$y / cellSize) * cellSize + (cellSize / 2))
      )

      # les observations sont positionnees sur une matrice. mIndices[i, j] == 1 indique qu'il y a au moins 1 observation pour le carreau (i,j)
      mIndices <- matrix(0L, max(dfObservations$i), max(dfObservations$j))
      mIndices[cbind(dfObservations$i, dfObservations$j)] <- 1L
      
      # construction d'une matrice des indices des centroides etendue au voisinage
      mIndicesEtendus <- matrix(0L, nrow(mIndices) + 2 * neighbor, ncol(mIndices) + 2 * neighbor)
      
      # decalage de la matrice d'index sur la matrice etendue, ce qui permet de compter combien de fois un carreau est necessaire
      for (voisin_x in -neighbor:neighbor)
      {
        for (voisin_y in -neighbor:neighbor)
        {
          mIndicesEtendus[neighbor + (1:nrow(mIndices)) + voisin_x, neighbor + (1:ncol(mIndices)) + voisin_y] <- 
            mIndicesEtendus[neighbor + (1:nrow(mIndices)) + voisin_x, neighbor + (1:ncol(mIndices)) + voisin_y] + mIndices
        }
      }
      
      # la matrice d'indices etendue est transformee en vecteurs de coordonnees
      vIndicesEtendus <- which(mIndicesEtendus > 0, arr.ind = TRUE)
      rm(list = c("mIndicesEtendus", "mIndices"))
      
      # retour aux coordonnees
      vX <- as.integer(round(min(dfCentroids$x) + (vIndicesEtendus[, 1] - 1 - neighbor) * cellSize))
      vY <- as.integer(round(min(dfCentroids$y) + (vIndicesEtendus[, 2] - 1 - neighbor) * cellSize))
      dfCentroidesUniques <- data.frame(x = vX, y = vY, i = vIndicesEtendus[, 1], j = vIndicesEtendus[, 2])
      rm(list = c("vX", "vY", "vIndicesEtendus", "dfCentroids"))
    }
    else
    {
      # calcul de l'indice des centroides
      dfCentroids$i <- as.integer(floor(dfCentroids$x / cellSize) - floor(min(dfCentroids$x / cellSize)) + 1)
      dfCentroids$j <- as.integer(floor(dfCentroids$y / cellSize) - floor(min(dfCentroids$y / cellSize)) + 1)
      dfCentroidesUniques <- dfCentroids
    }
      
    nomColonnes <- colnames(dfObservations)
    listVar <- nomColonnes[nomColonnes !="x" & nomColonnes !="y" & nomColonnes !="i" & nomColonnes !="j"]
    
    if(is.null(vQuantiles))
    {
      # numerotation des centroides - decalage de -1 pour faire commencer la numerotation des lignes a 0 pour le traitement c++
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
        ,dfObservations$y
        ,dfObservations$i
        ,dfObservations$j
        ,cellSize
        ,bandwidth
        ,as.matrix(dfObservations[, listVar])
        ,mXcentroides
        ,mYcentroides
        ,mIcentroides
        ,iNbCentroidesUniques
      )
      
      rm(list = c("dfObservations", "mXcentroides", "mYcentroides", "mIcentroides"))
      
      dfResultat <- cbind(dfResultat, mVariablesLissees)
      names(dfResultat) <- c("x","y", listVar)
      
      rm(mVariablesLissees)
      
      return(new(Class = "Grid", dfResultat, cellSize = cellSize, bandwidth = bandwidth))    
    }
    else
    {
      mMedianes <- rcppLissageMedianSort(
        dfObservations$x
        , dfObservations$y
        , bandwidth
        , as.matrix(dfObservations[, listVar])
        , dfCentroidesUniques$x
        , dfCentroidesUniques$y
        , vQuantiles
      )
    
      rm(dfObservations)
      
      dfResultat <- data.frame(cbind(dfCentroidesUniques$x, dfCentroidesUniques$y))
      dfResultat <- cbind(dfResultat, mMedianes)
      
      vNomQuantiles <- gsub("\\.", "", as.character(vQuantiles))
      
      dfListeVar <- expand.grid(vNomQuantiles, listVar)
      lNewVarNames <- paste(dfListeVar[, 2], "_", dfListeVar[, 1], sep = "")
      lNewVarNames <- c("x", "y", "nbObs", lNewVarNames)    
      colnames(dfResultat) <- lNewVarNames
      
      return(new(Class = "Grid", dfResultat, cellSize = cellSize, bandwidth = bandwidth))
    }
  }
