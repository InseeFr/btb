\docType{class}
\name{Grid-class}
\alias{Grid-class}
\title{Grid class (Classe Grid)}
\alias{[,Grid-method}
\alias{[<-,Grid-method}

\description{
A S4 class to define a grid (Une classe S4 pour definir un carroyage)
}

\value{
an object of class "Grid"
}

\section{Slots}{
\describe{
\item{\code{cellSize}}{Cell size of the grid (taille des carreaux)}
\item{\code{bandwidth}}{Radius of the Kernel Density Estimator (rayon de lissage)}
}}

\section{Methods}{
\describe{
\item{\code{[}}{Gets the value of an object (Recupere la valeur d un objet)}
\item{\code{[<-}}{Sets the value of an object (Modifie la valeur d un objet)}
}
}