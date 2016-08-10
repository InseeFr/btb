\name{Carroyage}
\alias{Carroyage}
\alias{Carroyage-class}
\alias{[,Carroyage-method}
\alias{[<-,Carroyage-method}

\docType{package}

\title{
Carroyage class (Classe Carroyage)
}

\description{
A S4 class to define a grid (Une classe S4 pour definir un carroyage)
}

\value{
Carroyage returns an object of class "Carroyage". 
}

\section{Slots}{
\describe{
\item{\code{taille_carreau}}{Cell size of the grid (Taille des carreaux)}
\item{\code{rayon_de_lissage}}{Radius of the Kernel Density Estimator (rayon de lissage)}
}
}

\section{Methods}{
\describe{
\item{\code{[}}{Gets the value of an object (Recupere la valeur d un objet)}
\item{\code{[<-}}{Sets the value of an object (Modifie la valeur d un objet)}
}
}
