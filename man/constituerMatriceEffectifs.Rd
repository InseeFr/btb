\name{constituerMatriceEffectifs}
\alias{constituerMatriceEffectifs}
\docType{package}

\title{Function constituting a membership matrix (Fonction constituant une matrice des effectifs)}

\description{
Function constituting a membership matrix. 

(Fonction constituant une matrice des effectifs.)
}

\usage{constituerMatriceEffectifs(vLigneObservation, vColonneObservation)}

\arguments{
  \item{vLigneObservation}{A \code{vector} containing the line number of each observation
  
  (Un \code{vector} contenant le numéro de ligne de chaque observation.)
  }
  
  \item{vColonneObservation}{A \code{vector} containing the column number of each observation
  
  (Un \code{vector} contenant le numéro de colonne de chaque observation.)
  }
  
}

\value{
Returns a \code{matrix} with the number of observations for each cell.

(Retourne une \code{matrix} avec le nombre d'observations pour chaque cellule.)
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

}