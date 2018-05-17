library(sf)
########################################################################################################################
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
# 2017/01/01      0.1.6      Arlindo Dos Santos     prise en compte du paramètre iNeighbor si indiqué (et 0 si userGrid fournie)
# 2017/01/11      0.1.7      Arlindo Dos Santos     fonction retour updateProgress pour smoothingToGrid
# 2017/01/11      0.1.8      Arlindo Dos Santos     ajout de la colonne iNbObsMin pour le lissage classique
# 2018/01/04      0.1.8      Arlindo Dos Santos     fusion des fonctions kernelSmoothing et smoothingToGrid
#
########################################################################################################################

#' @useDynLib btb
#' @importFrom Rcpp evalCpp
#' @import methods sf

############################## carreauxSF() #######################################################################
# fonction pour transformer une data.frame lissée en un objet sf
# carreauxSF <- function(df, iCellSize, sEPSG)
# {
#   r <- iCellSize / 2
#   df$geom <- sprintf("POLYGON ((%i %i, %i %i, %i %i, %i %i, %i %i))", df$x-r, df$y+r, df$x+r, df$y+r, df$x+r, df$y-r, df$x-r, df$y-r, df$x-r, df$y+r)
#   sfdf <- st_as_sf(df, wkt = "geom", crs = as.integer(sEPSG))
#   return(sfdf)
# }

############################## dfToGrid() #######################################################################
# arguments
# dfObservations  : data.frame comportant les coordonnées géographiques (x,y), ainsi que les variables que l'on souhaite lisser
# 
# retourne
# un objet Grid dont le slot @data contient la valeur des variables lissees
# 
#' @export
dfToGrid <- function(df, sEPSG, iCellSize = NULL)
{
  if(!is.null(iCellSize))
  {
    r = iCellSize / 2
    df$geometry <- sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))", df$x-r, df$y+r, df$x+r, df$y+r, df$x+r, df$y-r, df$x-r, df$y-r, df$x-r, df$y+r)
  }
  else
    df$geometry <- sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))", df$x-df$iCellSize / 2, df$y+df$iCellSize / 2, df$x+df$iCellSize / 2, df$y+df$iCellSize / 2, df$x+df$iCellSize / 2, df$y-df$iCellSize / 2, df$x-df$iCellSize / 2, df$y-df$iCellSize / 2, df$x-df$iCellSize / 2, df$y+df$iCellSize / 2)
  
  sfdf <- st_as_sf(df, wkt = "geometry", crs = as.integer(sEPSG))
  return(sfdf)
}

############################## kernelSmoothing() #######################################################################
# arguments
# dfObservations  : data.frame comportant les coordonnées géographiques (x,y), ainsi que les variables que l'on souhaite lisser
# iCellSize        : Taille des carreaux
# iBandwidth       : Rayon de lissage
# iNeighbor        : Paramètre technique pour calculer l'étendue des points d'estimation
# vQuantiles      : vecteur de quantiles à utiliser pour le lissage median
# dfCentroids     : data.frame avec deux colonnes, nommées x et y avec les coordonnees des centroides à utiliser
# fUpdateProgress : fonction permettant d'offrir à l'appelant une estimation de l'avancement du traitement 
# iNeighbor        : Paramètre technique pour calculer l'étendue des points d'estimations, à ne pas remplir
# 
# retourne
# un objet Grid dont le slot @data contient la valeur des variables lissees
#
# Remarque : les coordonnées (x;y) des observations sont implicitements arrondies à l'unité. (lors de l'appel aux fct c++)
# Ceci est justifié par le fait que :
# - la géolocalisation est souvent imprécise en dessous de 1 mètre 
# - cela permet d'économiser de la mémoire dans les fct développées en C++
# 
#' @export
kernelSmoothing <-
  function(dfObservations
           , sEPSG
           , iCellSize
           , iBandwidth
           , vQuantiles = NULL
           , dfCentroids = NULL
           , fUpdateProgress = NULL
           , iNeighbor = NULL
           , iNbObsMin = 250
  )
  {
    iCellSize <- as.integer(iCellSize)
    iBandwidth <- as.integer(iBandwidth)
    dRayonMinimum <- iCellSize * sqrt(2) / 2
    
    if(is.null(iNeighbor))
    {
      if (is.null(dfCentroids)) 
        iNeighbor <- max(0, ceiling(iBandwidth / iCellSize / 2L) - 1L) 
      else 
        iNeighbor <- 0
    }
    
    if (iBandwidth < dRayonMinimum)
      stop("iBandwidth must be greater than iCellSize * sqrt(2) / 2")
    
    if (iCellSize <= 0)
      stop("iCellSize must be greater than 0")

    for(colname in colnames(dfObservations))
    {
      iNbNA <- sum(is.na(dfObservations[, colname]))
      if(iNbNA > 0)
      {
        if(colname == "x" || colname == "y")
        {
          # les coordonnées des observations ne doivent pas être NA
          stop("NA coordinates are not allowed")
        }
        else
        { 
          warning("Be careful! NA values detected in your observations\n")
          break;
        }
      }
    }
    
    # vérifier la regularite des centroides fournis par l'utilisateur
    if (!is.null(dfCentroids))
    {
      xOffset <- (dfCentroids$x + iCellSize / 2) %% iCellSize
      yOffset <- (dfCentroids$y + iCellSize / 2) %% iCellSize
      if (!all(xOffset == xOffset[1]) | !all(yOffset == yOffset[1]) )
        stop("Centroids are not regular")
    }else
    {
      xOffset <- 0
      yOffset <- 0
    }
    
    # vérifier la validité des quantiles fournis
    if (!is.null(vQuantiles))
    {
      if (min(vQuantiles) <= 0 || max(vQuantiles) >= 1)
        stop("Invalid quantiles values")
    }    
    
    if (is.null(dfCentroids))
    { 
      # calcul de l'indice des observations - on prend le rectangle englobant et on positionne le debut de la numérotation sur la première observation
      dfObservations$col <- as.integer(floor((dfObservations$x - xOffset[1]) / iCellSize) - floor(min(dfObservations$x / iCellSize)) + 1)
      dfObservations$row <- as.integer(floor((dfObservations$y - yOffset[1]) / iCellSize) - floor(min(dfObservations$y / iCellSize)) + 1)
      
      # calcul des centroides
      dfCentroids <- data.frame( x = as.integer(floor(dfObservations$x / iCellSize) * iCellSize + (iCellSize / 2)),
                                 y = as.integer(floor(dfObservations$y / iCellSize) * iCellSize + (iCellSize / 2))
      )
      
      # les observations sont positionnées sur une matrice. mIndices[col, row] == 1 indique qu'il y a au moins 1 observation pour le carreau (col, row)
      mIndices <- matrix(0L, max(dfObservations$col), max(dfObservations$row))
      mIndices[cbind(dfObservations$col, dfObservations$row)] <- 1L
      
      # construction d'une matrice des indices des centroides étendue au voisinage
      mIndicesEtendus <- matrix(0L, nrow(mIndices) + 2 * iNeighbor, ncol(mIndices) + 2 * iNeighbor)
      
      # décalage de la matrice d'index sur la matrice étendue, ce qui permet de compter combien de fois un carreau est nécessaire
      for (voisin_x in -iNeighbor:iNeighbor)
      {
        for (voisin_y in -iNeighbor:iNeighbor)
        {
          mIndicesEtendus[iNeighbor + 1:nrow(mIndices) + voisin_x, iNeighbor + 1:ncol(mIndices) + voisin_y] <- 
            mIndicesEtendus[iNeighbor + 1:nrow(mIndices) + voisin_x, iNeighbor + 1:ncol(mIndices) + voisin_y] + mIndices
        }
      }
      
      # la matrice d'indices étendue est transformée en vecteurs de coordonnées
      vIndicesEtendus <- which(mIndicesEtendus > 0, arr.ind = TRUE)
      rm(list = c("mIndicesEtendus", "mIndices"))
      
      # retour aux coordonnées
      vX <- as.integer(round(min(dfCentroids$x) + (vIndicesEtendus[, 1] - 1 - iNeighbor) * iCellSize))
      vY <- as.integer(round(min(dfCentroids$y) + (vIndicesEtendus[, 2] - 1 - iNeighbor) * iCellSize))
      dtCentroidesUniques <- data.frame(x = vX, y = vY, col = vIndicesEtendus[, 1], row = vIndicesEtendus[, 2])
      
      rm(list = c("vX", "vY", "vIndicesEtendus", "dfCentroids"))
      
    }else
    {
      # Remarque: il n'est pas nécessaire de rapatrier les observations dans les carreaux de la grille fournie.
      # lors du lissage, les observations enverront leur contribution uniquement vers les carreaux fournis
      
      # calcul de l'indice des observations - on commence la numérotation pour la coordonnée minimale (qu'elle soit détenue par une observation ou un centroide)
      obsEtCentroides <- data.frame(x = c(dfObservations$x, dfCentroids$x), y = c(dfObservations$y, dfCentroids$y))
      indiceMinX <- floor(min(obsEtCentroides$x / iCellSize))
      indiceMinY <- floor(min(obsEtCentroides$y / iCellSize))
      
      dfObservations$col <- as.integer(floor((dfObservations$x - xOffset[1]) / iCellSize) - indiceMinX + 1)
      dfObservations$row <- as.integer(floor((dfObservations$y - yOffset[1]) / iCellSize) - indiceMinY + 1)
      
      # calcul de l'indice des centroides
      dfCentroids$col <- as.integer(floor(dfCentroids$x / iCellSize) - indiceMinX + 1)
      dfCentroids$row <- as.integer(floor(dfCentroids$y / iCellSize) - indiceMinY + 1)
      dtCentroidesUniques <- dfCentroids
      rm(dfCentroids)
      rm(obsEtCentroides)
    }
    
    nomColonnes <- colnames(dfObservations)
    listVar <- nomColonnes[nomColonnes != "x" & nomColonnes != "y" & nomColonnes != "col" & nomColonnes != "row"]
    
    if (is.null(vQuantiles))
    {
      # numérotation des centroides - décalage de -1 pour faire commencer la numerotation des lignes à 0 pour le traitement c++
      dtCentroidesUniques$index <- (1:nrow(dtCentroidesUniques)) - 1L
      
      # transformation en matrice
      mIcentroides = matrix(-1, max(dtCentroidesUniques$row), max(dtCentroidesUniques$col))
      mIcentroides[cbind(dtCentroidesUniques$row, dtCentroidesUniques$col)] <- dtCentroidesUniques$index
      
      iNbCentroidesUniques <- nrow(dtCentroidesUniques)
      dfResultat <- data.frame(dtCentroidesUniques[, c("x", "y")])
      
      mVariablesLissees <- rcppLissage(
          dfObservations$x
        , dfObservations$y
        , dfObservations$row + iNeighbor
        , dfObservations$col + iNeighbor
        , iCellSize
        , iBandwidth
        , iNeighbor
        , as.matrix(dfObservations[, listVar])
        , max(dtCentroidesUniques$col)
        , max(dtCentroidesUniques$row)
        , min(dtCentroidesUniques$x)
        , min(dtCentroidesUniques$y)
        , mIcentroides
        , iNbCentroidesUniques
        , fUpdateProgress
      )
      
      rm(list = c("dfObservations", "dtCentroidesUniques", "mIcentroides"))
      
      dfResultat <- cbind(dfResultat, mVariablesLissees)
      names(dfResultat) <- c("x", "y", listVar)
      # names(dfResultat) <- c("x","y", listVar, "nbObsPondere") # version pour calcul de la colonne nbObsPondere
      
      rm(mVariablesLissees)
     
      return(dfToGrid(dfResultat, sEPSG = sEPSG, iCellSize = iCellSize))
    }else
    {
      dfObservations$col <- as.integer(floor((dfObservations$x - xOffset[1]) / iCellSize) - floor(min(dfObservations$x / iCellSize)) + 1)
      dfObservations$row <- as.integer(floor((dfObservations$y - yOffset[1]) / iCellSize) - floor(min(dfObservations$y / iCellSize)) + 1)
      
      dtCentroidesUniques$col <- as.integer(floor((dtCentroidesUniques$x - xOffset[1]) / iCellSize) - floor(min(dtCentroidesUniques$x / iCellSize)) + 1)
      dtCentroidesUniques$row <- as.integer(floor((dtCentroidesUniques$y - yOffset[1]) / iCellSize) - floor(min(dtCentroidesUniques$y / iCellSize)) + 1)
      
      mMedianes <- rcppLissageMedianGrappe(
          iNbObsMin
        , dfObservations$x
        , dfObservations$y
        , dfObservations$row + iNeighbor
        , dfObservations$col + iNeighbor
        , iCellSize
        , iBandwidth
        , as.matrix(dfObservations[, listVar])
        , dtCentroidesUniques$x
        , dtCentroidesUniques$y
        , dtCentroidesUniques$row - 1L
        , dtCentroidesUniques$col - 1L
        , vQuantiles
      )
      
      rm(dfObservations)
      rm(dtCentroidesUniques)
      
      vNomQuantiles <- gsub("\\.", "", as.character(vQuantiles))
      
      dfListeVar <- expand.grid(vNomQuantiles, listVar)
      lNewVarNames <- paste0(dfListeVar[, 2], "_", dfListeVar[, 1])
      lNewVarNames <- c("nbObs", lNewVarNames, "x", "y")    
      colnames(mMedianes) <- lNewVarNames
      
      return(dfToGrid(data.frame(mMedianes), sEPSG = sEPSG, iCellSize = iCellSize))
    }
  }

# .onUnload <- function (libpath) {
#   library.dynam.unload("btb", libpath)
# }
