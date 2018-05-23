\name{constituerGrappes}
\alias{constituerGrappes}
\docType{package}

\title{Quadtree variant function (Variante de la fonction Quadtree)}

\description{
Quadtree variant function. 

(Variante de la fonction Quadtree.)
}

\usage{constituerGrappes(iNbObsMin, mEffectifs)}

\arguments{
  \item{iNbObsMin}{An \code{integer} representing the minimum number of elements in each cluster
  
  (Un \code{data.frame} représentant le nombre minimum d'éléments dans chaque grappe.)}
  
  \item{mEffectifs}{staffing matrix for each cell (\code{matrix}). 
  
  (matrice des effectifs pour chaque case (\code{matrix}). .)}

}

\value{
Returns a \code{matrix} with the cluster number for each cell . 

(Retourne une \code{matrix} avec le numéro de grappe pour chaque cellule.)
} 
 
\author{Psar Analyse Urbaine Insee - Arlindo Dos Santos and Francois Semecurbe}

\examples{
dfObservations <- data.frame(x = c(15, 35, 15, 25, 35, 55, 45, 45, 55, 65, 70, 75, 85, 90,
                                    65, 75, 85, 65, 70, 75, 85, 90, 65, 70, 75)
                          ,  y = c(10, 10, 30, 30, 35, 35, 45, 55, 55, 65, 65, 65, 65, 65,
                                    70, 70, 70, 75, 75, 75, 75, 75, 85, 85, 85)
)

cellSize <- 20L

# calcul de l'indice des observations 
# on prend le rectangle englobant 
# et on positionne le debut de la numérotation sur la première observation
dfObservations$col <- as.integer(floor((dfObservations$x) / cellSize)
                               - floor(min(dfObservations$x / cellSize)) + 1)
dfObservations$row <- as.integer(floor((dfObservations$y) / cellSize) 
                               - floor(min(dfObservations$y / cellSize)) + 1)

mEffectifs <- constituerMatriceEffectifs(dfObservations$row - 1, dfObservations$col - 1)

#### matrice des grappes
mGrappes <- constituerGrappes(1, mEffectifs)
}