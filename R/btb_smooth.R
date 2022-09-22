
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
#' @title btb_smooth 
#' 
#' @description 
#' Smoothing function with a bisquare kernel or median.
#' 
#' (Fonction de lissage à partir d'un noyau bisquare ou de la médiane.)
#' 
#' @usage 
#' 
#' # Call mode 1: bisquare kernel smoothing - automatic grid
#'btb_smooth(  dfObservations
#'                  , sEPSG
#'                  , iCellSize
#'                  , iBandwidth
#'                  , vQuantiles = NULL
#'                  , dfCentroids = NULL
#'                  , fUpdateProgress = NULL
#'                  , iNeighbor = NULL
#'                  , iNbObsMin = 250
#')
#'
#'# Call mode 2: median smoothing - automatic grid
#'btb_smooth(  dfObservations
#'                  , sEPSG
#'                  , iCellSize
#'                  , iBandwidth
#'                  , vQuantiles
#'                  , dfCentroids = NULL
#'                  , fUpdateProgress = NULL
#'                  , iNeighbor = NULL
#'                  , iNbObsMin = 250
#')
#'
#'# Call mode 3: bisquare kernel smoothing - user grid
#'btb_smooth(  dfObservations
#'                  , sEPSG
#'                  , iCellSize
#'                  , iBandwidth
#'                  , vQuantiles = NULL
#'                  , dfCentroids
#'                  , fUpdateProgress = NULL
#'                  , iNeighbor = NULL
#'                  , iNbObsMin = 250
#')
#'
#'# Call mode 4: median smoothing - user grid
#'btb_smooth(  dfObservations
#'                  , sEPSG
#'                  , iCellSize
#'                  , iBandwidth
#'                  , vQuantiles
#'                  , dfCentroids
#'                  , fUpdateProgress = NULL
#'                  , iNeighbor = NULL
#'                  , iNbObsMin = 250
#')
#' 
#' @param dfObservations 
#' A `data.frame` with cartesian geographical coordinates and variables to smooth. (x, y, var1, var2, ...)
#'  
#'  (Un `data.frame` comportant les coordonnées géographiques cartésiennes (x,y), ainsi que les variables que l'on souhaite lisser. (x, y, var1, var2, ...)
#'  
#' @param sEPSG 
#' EPSG code of projection (`character`). For example, the RGF93 / Lambert-93 projection has "2154" code.
#' 
#' (code EPSG de la projection (`character`). Par exemple, la projection RGF93 / Lambert-93 a pour code "2154".)`
#'  
#' @param iCellSize
#' Cell size of the grid (`integer`). The unit of measurement is free. It must be the same as the unit of 
#' `iBandwidth` variable. 
#' 
#' (Taille des carreaux (`integer`). Le choix de l'unité de mesure est laissé libre à l'utilisateur. Elle doit seulement être la même que celle de la variable `iBandwidth`.)
#'  
#' @param iBandwidth
#' Radius of the Kernel Density Estimator (`integer`). This bandwidth acts as a smoothing parameter, controlling the balance between bias and variance. A large bandwidth leads to a very smooth (i.e. high-bias) density distribution. A small bandwidth leads to an unsmooth (i.e. high-variance) density distribution. The unit of measurement is free. It must be the same as the unit of `iCellSize` variable. 
#' 
#' (Rayon de lissage de l'estimation d'intensité par noyau (`integer`). Cette bande-passante se comporte comme un paramètre de lissage, controlant l'équilibre entre biais et variance. Un rayon élevé conduit à une densité tres lissée, avec un biais élevé. Un petit rayon génère une densité peu lissée avec une forte variance. Le choix de l'unité de mesure est laissé libre à l'utilisateur. Elle doit seulement être la même que celle de la variable `iCellSize`.
#' 
#' @param vQuantiles 
#' Percentile vector to calculate. For example c(0.1, 0.25, 0.5) will calculate the first decile, the first quartile and the median. 
#'  
#' (Vecteur des quantiles à calculer. Par exemple c(0.1, 0.25, 0.5) retournera le premier décile, le premier quartile et la mediane.)`
#' 
#' @param dfCentroids
#' A `data.frame` with two columns (x, y) containing coordinates of the user's centroids. The coordinates must be in the same projection than (`dfObservations`).
#'  
#' (Un `data.frame` avec deux colonnes (x, y) contenant les coordonnées des centroides de l'utilisateur. Les coordonnées doivent être dans le même système de coordonnées que (`dfObservations`).)
#'
#' @param fUpdateProgress
#' A `function` to see compute progress. 
#'  
#' (Une `function` pour voir la progression du calcul.)
#'
#' @param iNeighbor
#' Technical parameter, leave empty. (`integer`) 
#'  
#' (Paramètre technique pour calculer l'étendue des points d'estimations, à ne pas remplir. (`integer`))
#' 
#' @param iNbObsMin
#' Minimum size of constituted grappes for median smoothing.  (`integer`) 
#'  
#' (Taille minimale des grappes constituées pour le lissage "médian" (géographiquement pondéré). (`integer`))
#' 
#' @details 
#' Returns an object inheriting from the `data.frame` class. (Retourne un objet qui se comporte comme un `data.frame`, par heritage.)
#' 
#' - Smoothing covers a set of methods to extract pertinent and structuring information from noisy data.  In the field of spatial analysis, and most widely in quantitative geography, smoothing is used to modelise density variations of a population distribution in geographical space. Kernel smoothing methods are widely used.
#' In this method, for each location x, we count the number of events of a process within a distance h of x, and weighted by the square reciprocal of the radius h. We apply a edge-correction to deal with edge-effects. So the method is conservative..
#' 
#' - Le lissage recouvre un ensemble de méthodes pour extraire d'une source de données bruitées une information pertinente et structurante. Dans le champ de l'analyse spatiale et plus largement de la géographie quantitative, le lissage est principalement utilisé pour modéliser les variations de densités d'une distribution de population dans l'espace géographique. On utilise principalement des méthodes de lissage par noyau.
#' Il s'agit ici, pour chaque point x, de comptabliser le nombre d' "évènements" d'un processus à une distance h de ce point, tout en ponderant ce nombre par l'inverse de la distance h au carré. On applique une correction à la ponderation afin de traiter les effets de bord. Cette méthode est conservative.
#' 
#' @references 
#' 
#' - "Geographically weighted summary statistics : a framework for localised exploratory data analysis", C.Brunsdon & al., in Computers, Environment and Urban Systems 2002
#' - Statistical Analysis of Spatial and Spatio-Temporal Point Patterns, Third Edition, Diggle, 2003, pp. 83-86
#' 
#' @author 
#' - Psar Analyse Urbaine Insee 
#' - Thierry Cornely
#' - Laure Genebes
#' - Arlindo Dos Santos
#' - Cynthia Faivre
#' - Auriane Renaud 
#' - Francois Semecurbe
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' ######### example 1 #########
#' data(dfPrix_SP95_2016)
#' dfPrix_SP95_2016$nbObs <- 1L
#' dfSmoothed <- btb_smooth(dfObservations = dfPrix_SP95_2016
#'                               , sEPSG = "2154"
#'                               , iCellSize = 5000L
#'                               , iBandwidth = 30000L)
#' dfSmoothed$prix95 <- dfSmoothed$SP95 / dfSmoothed$nbObs * 100
#' library(cartography)
#' choroLayer(dfSmoothed
#'            , var = "prix95"
#'            , nclass = 5
#'            , method = "fisher-jenks"
#'            , border = NA
#'            , legend.title.txt = "prix du SP95 en centimes")
#' ######### example 2 #########
#' library(sp)
#' library(cartography)
#' data(reunion)
#' # Smoothing all variables for Reunion (Lissage de toutes les variables pour la Reunion)
#' # Call mode 1: classic smoothing - automatic grid
#' reunionSmoothed <- btb_smooth( dfObservations = reunion
#'                                     , sEPSG = "32740"
#'                                     , iCellSize = 200L
#'                                     , iBandwidth = 400L)
#' # preview (Apercu)
#' choroLayer(reunionSmoothed, var = "houhold", nclass = 5, method = "fisher-jenks", border = NA)
#' # Call mode 2: median smoothing - automatic grid
#' reunionSmoothed <- btb_smooth( dfObservations = reunion
#'                                     , sEPSG = "32740"
#'                                     , iCellSize = 200L
#'                                     , iBandwidth = 400L
#'                                     , vQuantiles = c(0.1, 0.5, 0.9))
#' # preview (Apercu)
#' choroLayer(reunionSmoothed, var = "houhold_05", nclass = 5, method = "fisher-jenks", border = NA)
#' # Call mode 3: classic smoothing - user grid
#' dfCentroidsUser <- merge( x = seq(from =  314400L, to =  378800L, by = 200L)
#'                           , y = seq(from = 7634000L, to = 7691200L, by = 200L))
#' reunionSmoothed <- btb_smooth( dfObservations = reunion
#'                                     , sEPSG = "32740"
#'                                     , iCellSize = 200L
#'                                     , iBandwidth = 400L
#'                                     , dfCentroids = dfCentroidsUser)
#' # preview (Apercu)
#' reunionSmoothed <- reunionSmoothed[reunionSmoothed$houhold > 0, ]
#' choroLayer(reunionSmoothed, var = "houhold", nclass = 5, method = "fisher-jenks", border = NA)
#' # Call mode 4: median smoothing - user grid
#' reunionSmoothed <- btb_smooth( dfObservations = reunion
#'                                     , sEPSG = "32740"
#'                                     , iCellSize = 200L
#'                                     , iBandwidth = 400L
#'                                     , vQuantiles = c(0.1, 0.5, 0.9)
#'                                     , dfCentroids = dfCentroidsUser)
#' # preview (Apercu)
#' reunionSmoothed <- reunionSmoothed[reunionSmoothed$nbObs > 0, ]
#' choroLayer(reunionSmoothed, var = "houhold_05", nclass = 5, method = "fisher-jenks", border = NA)
#' }




btb_smooth <-
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
    # Numeric format
    iCellSize <- as.integer(iCellSize)
    iBandwidth <- as.integer(iBandwidth)
    
    
    # *************************************************************
    #                           Checks
    # *************************************************************
    
    # Number of neighbors pixels 
    if(is.null(iNeighbor))
    {
      if (is.null(dfCentroids)) 
        iNeighbor <- max(0, ceiling(iBandwidth / iCellSize / 2L) - 1L) 
      else 
        iNeighbor <- 0
    }
    #/!\ what if !is.null(dfCentroids) & !is.null(iNeighbor) ?
    
    # Controling size of iBandwidth
    if (iBandwidth < iCellSize * sqrt(2) / 2)
      stop("iBandwidth must be greater than iCellSize * sqrt(2) / 2")
    
    if (iCellSize <= 0)
      stop("iCellSize must be greater than 0")
    
    # No NA values accepted for coordinates
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
    
    
    # *********************************************************
    #  Launching smoothing with rcpp btb functions
    # *********************************************************
    # - Creates inputs for rcppLissage or rcppLissageMedianGrappe
    # - lauches rcppLissage / rcppLissageMedianGrappe
    # - Returns a grid with results (using btb_ptsToGrid)
    
    
    if (is.null(dfCentroids)) # Building automatic grid
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
    
    # List of numeric variables to smooth
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
     
      return(btb::btb_ptsToGrid(dfResultat, sEPSG = sEPSG, iCellSize = iCellSize))
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
      
      return(btb::btb_ptsToGrid(data.frame(mMedianes), sEPSG = sEPSG, iCellSize = iCellSize))
    }
  }

# .onUnload <- function (libpath) {
#   library.dynam.unload("btb", libpath)
# }



############################## carreauxSF() #######################################################################
# fonction pour transformer une data.frame lissée en un objet sf
# carreauxSF <- function(df, iCellSize, sEPSG)
# {
#   r <- iCellSize / 2
#   df$geom <- sprintf("POLYGON ((%i %i, %i %i, %i %i, %i %i, %i %i))", df$x-r, df$y+r, df$x+r, df$y+r, df$x+r, df$y-r, df$x-r, df$y-r, df$x-r, df$y+r)
#   sfdf <- st_as_sf(df, wkt = "geom", crs = as.integer(sEPSG))
#   return(sfdf)
# }



# iCellSize <- 20L
# iBandwidth <- 41L
# dfObservations <- data.frame(x = c(22, 35), y = c(70, 55), V1 = c(10, 13))
# vQuantiles = NULL
# dfCentroids = NULL
# fUpdateProgress = NULL
# iNeighbor = NULL
# iNbObsMin = 250
# 
# vXCentroides <- rep(seq(from = 10, to = 50, by = iCellSize), 4)
# vYCentroides <- rep(seq(from = 30, to = 90, by = iCellSize), each = 3)
# dfCentroids <- data.frame(cbind(x = vXCentroides, y = vYCentroides))

