#' @title Smoothing with a bisquare kernel or median
#' @name btb_smooth 
#' @description 
#' Smoothing function with a bisquare kernel or median.
#' 
#' (Fonction de lissage à partir d'un noyau bisquare ou de la médiane.)
#' 
#'  
#' @param pts 
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
#' A `data.frame` with two columns (x, y) containing coordinates of the user's centroids. The coordinates must be in the same projection than (`pts`).
#'  
#' (Un `data.frame` avec deux colonnes (x, y) contenant les coordonnées des centroides de l'utilisateur. Les coordonnées doivent être dans le même système de coordonnées que (`pts`).)
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
#' @param inspire (boolean) : if TRUE, returns a column for Inspire grid names.
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
#' 
#' @export
#' 
#' 
#' @examples 
#' \dontrun{
#' # Example 1
#' data(dfPrix_SP95_2016)
#' dfPrix_SP95_2016$nbObs <- 1L
#' dfSmoothed <- btb::btb_smooth(pts = dfPrix_SP95_2016,
#'                               sEPSG = "2154",
#'                               iCellSize = 5000L,
#'                               iBandwidth = 30000L, 
#'                               inspire = TRUE)
#' dfSmoothed$prix95 <- dfSmoothed$SP95 / dfSmoothed$nbObs * 100
#' library(mapsf)
#' mf_map(dfSmoothed,
#'       type = "choro",
#'       var = "prix95",
#'       breaks = "fisher",
#'       nbreaks = 5,
#'       border = NA,
#'       leg_title = "prix du SP95 en centimes")
#' # Example 2
#' data(reunion)
#' # Call mode 1: classic smoothing - automatic grid
#' reunionSmoothed <- btb_smooth( pts = reunion,
#'                                     sEPSG = "32740",
#'                                     iCellSize = 200L,
#'                                     iBandwidth = 400L)
#' library(mapsf)
#' mf_map(reunionSmoothed,
#'       type = "choro",
#'       var = "houhold",
#'       breaks = "fisher",
#'       nbreaks = 5,
#'       border = NA)
#' # Call mode 2: median smoothing - automatic grid
#' reunionSmoothed <- btb_smooth( pts = reunion,
#'                                      sEPSG = "32740",
#'                                      iCellSize = 200L,
#'                                      iBandwidth = 400L,
#'                                      vQuantiles = c(0.1, 0.5, 0.9))
#' mf_map(reunionSmoothed,
#'       type = "choro",
#'       var = "houhold_05",
#'       breaks = "fisher",
#'       nbreaks = 5,
#'       border = NA)
#' # Call mode 3: classic smoothing - user grid
#' dfCentroidsUser <- merge( x = seq(from =  314400L, to =  378800L, by = 200L),
#'                           y = seq(from = 7634000L, to = 7691200L, by = 200L))
#' reunionSmoothed <- btb_smooth( pts = reunion,
#'                                     sEPSG = "32740",
#'                                     iCellSize = 200L,
#'                                     iBandwidth = 400L,
#'                                     dfCentroids = dfCentroidsUser)
#' reunionSmoothed <- reunionSmoothed[reunionSmoothed$houhold > 0, ]
#' mf_map(reunionSmoothed,
#'       type = "choro",
#'       var = "houhold",
#'       breaks = "fisher",
#'       nbreaks = 5,
#'       border = NA)
#' # Call mode 4: median smoothing - user grid
#' reunionSmoothed <- btb_smooth( pts = reunion,
#'                                     sEPSG = "32740",
#'                                     iCellSize = 200L,
#'                                     iBandwidth = 400L,
#'                                     vQuantiles = c(0.1, 0.5, 0.9),
#'                                     dfCentroids = dfCentroidsUser)
#' reunionSmoothed <- reunionSmoothed[reunionSmoothed$nbObs > 0, ]
#' mf_map(reunionSmoothed,
#'       type = "choro",
#'       var = "houhold_05",
#'       breaks = "fisher",
#'       nbreaks = 5,
#'       border = NA)
#' }




btb_smooth <-
  function(pts ,
           sEPSG = NA,
           iCellSize = NA ,
           iBandwidth ,
           vQuantiles = NULL ,
           dfCentroids = NULL ,
           iNeighbor = NULL,
           inspire = F,
           iNbObsMin = 250
  )
  {
    # Numeric format
    iCellSize <- as.integer(iCellSize)
    iBandwidth <- as.integer(iBandwidth)
    
    # If pts is a sf objet
    if("sf" %in% class(pts)){
      stopifnot("sf pts : geometry type must be POINT"= identical(as.vector(sf::st_geometry_type(pts,by_geometry=F)),"POINT") )
      
      # Retrieve epsg code from sf points (or NA)
      if(is.na(sEPSG)) sEPSG <- (sf::st_crs(pts))$epsg
      
      # As df with coordinates
      pts$x <- sf::st_coordinates(pts)[,1]
      pts$y <- sf::st_coordinates(pts)[,2]
      pts <- sf::st_drop_geometry(pts) 
      
      stopifnot("sEPSG not found" = !is.na(sEPSG))
      stopifnot("sEPSG not valid "= identical(nchar(as.character(sEPSG)),4L))
    }
    
    
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
    for(colname in colnames(pts))
    {
      iNbNA <- sum(is.na(pts[, colname]))
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
    
    # Absence de variables characters
    vartest <- colnames(pts)[!colnames(pts) %in% c("x","y")]
    if(length(vartest)>0){
      nb_charvar <- pts[,vartest] %>% lapply(is.character) %>% unlist() %>% sum()
      stopifnot("No character variables permitted"=identical(nb_charvar,0L))
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
      pts$col <- as.integer(floor((pts$x - xOffset[1]) / iCellSize) - floor(min(pts$x / iCellSize)) + 1)
      pts$row <- as.integer(floor((pts$y - yOffset[1]) / iCellSize) - floor(min(pts$y / iCellSize)) + 1)
      
      # calcul des centroides
      # dfCentroids <- data.frame( x = as.integer(floor(pts$x / iCellSize) * iCellSize + (iCellSize / 2)),
      #                            y = as.integer(floor(pts$y / iCellSize) * iCellSize + (iCellSize / 2))
      # )
      dfCentroids <- btb::btb_add_centroids(pts,iCellSize,names_centro = c("x","y"),add=F)
      
      # les observations sont positionnées sur une matrice. mIndices[col, row] == 1 indique qu'il y a au moins 1 observation pour le carreau (col, row)
      mIndices <- matrix(0L, max(pts$col), max(pts$row))
      mIndices[cbind(pts$col, pts$row)] <- 1L
      
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
      obsEtCentroides <- data.frame(x = c(pts$x, dfCentroids$x), y = c(pts$y, dfCentroids$y))
      indiceMinX <- floor(min(obsEtCentroides$x / iCellSize))
      indiceMinY <- floor(min(obsEtCentroides$y / iCellSize))
      
      pts$col <- as.integer(floor((pts$x - xOffset[1]) / iCellSize) - indiceMinX + 1)
      pts$row <- as.integer(floor((pts$y - yOffset[1]) / iCellSize) - indiceMinY + 1)
      
      # calcul de l'indice des centroides
      dfCentroids$col <- as.integer(floor(dfCentroids$x / iCellSize) - indiceMinX + 1)
      dfCentroids$row <- as.integer(floor(dfCentroids$y / iCellSize) - indiceMinY + 1)
      dtCentroidesUniques <- dfCentroids
      rm(dfCentroids)
      rm(obsEtCentroides)
    }
    
    # List of numeric variables to smooth
    nomColonnes <- colnames(pts)
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
          pts$x
        , pts$y
        , pts$row + iNeighbor
        , pts$col + iNeighbor
        , iCellSize
        , iBandwidth
        , iNeighbor
        , as.matrix(pts[, listVar])
        , max(dtCentroidesUniques$col)
        , max(dtCentroidesUniques$row)
        , min(dtCentroidesUniques$x)
        , min(dtCentroidesUniques$y)
        , mIcentroides
        , iNbCentroidesUniques
        , NA #old message update param in rcppLissage
      )
      
      rm(list = c("pts", "dtCentroidesUniques", "mIcentroides"))
      
      dfResultat <- cbind(dfResultat, mVariablesLissees)
      names(dfResultat) <- c("x", "y", listVar)
      # names(dfResultat) <- c("x","y", listVar, "nbObsPondere") # version pour calcul de la colonne nbObsPondere
      
      rm(mVariablesLissees)
     
      return(btb::btb_ptsToGrid(dfResultat, sEPSG = sEPSG, iCellSize = iCellSize,names_centro = c("x","y"), inspire = inspire))
    }else
    {
      pts$col <- as.integer(floor((pts$x - xOffset[1]) / iCellSize) - floor(min(pts$x / iCellSize)) + 1)
      pts$row <- as.integer(floor((pts$y - yOffset[1]) / iCellSize) - floor(min(pts$y / iCellSize)) + 1)
      
      dtCentroidesUniques$col <- as.integer(floor((dtCentroidesUniques$x - xOffset[1]) / iCellSize) - floor(min(dtCentroidesUniques$x / iCellSize)) + 1)
      dtCentroidesUniques$row <- as.integer(floor((dtCentroidesUniques$y - yOffset[1]) / iCellSize) - floor(min(dtCentroidesUniques$y / iCellSize)) + 1)
      
      mMedianes <- rcppLissageMedianGrappe(
          iNbObsMin
        , pts$x
        , pts$y
        , pts$row + iNeighbor
        , pts$col + iNeighbor
        , iCellSize
        , iBandwidth
        , as.matrix(pts[, listVar])
        , dtCentroidesUniques$x
        , dtCentroidesUniques$y
        , dtCentroidesUniques$row - 1L
        , dtCentroidesUniques$col - 1L
        , vQuantiles
      )
      
      rm(pts)
      rm(dtCentroidesUniques)
      
      vNomQuantiles <- gsub("\\.", "", as.character(vQuantiles))
      
      dfListeVar <- expand.grid(vNomQuantiles, listVar)
      lNewVarNames <- paste0(dfListeVar[, 2], "_", dfListeVar[, 1])
      lNewVarNames <- c("nbObs", lNewVarNames, "x", "y")    
      colnames(mMedianes) <- lNewVarNames
      
      return(btb::btb_ptsToGrid(data.frame(mMedianes), sEPSG = sEPSG, iCellSize = iCellSize,names_centro = c("x","y"), inspire = inspire))
    }
  }


