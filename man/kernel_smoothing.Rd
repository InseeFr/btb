\name{kernel_smoothing}
\alias{kernel_smoothing}
\docType{package}

\title{
Smoothing function (Fonction de lissage)
}

\description{
Smoothing function with a bisquare kernel. (Fonction de lissage a partir d un noyau bisquare.)
}

\usage{
  kernel_smoothing(df,cell_size,bandwith,list_var,neighbor)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{A \code{data.frame} with geographical coordinates (x,y) and variables to smooth. 
  
  (Un \code{data.frame} comportant les coordonnees geographiques (x,y), ainsi que les variables que l on souhaite lisser.)}
  
  \item{cell_size}{Cell size of the grid. The unit of measurement is free. It must be the same as the unit of 
  \code{bandwith} variable. 
  
  (Taille des carreaux. Le choix de l unite de mesure est laisse libre a l utilisateur. Elle doit seulement etre la meme que celle de la variable \code{bandwith}.)}
  
  \item{bandwith}{Radius of the Kernel Density Estimator. This bandwidth acts as a smoothing parameter, controlling the balance between bias and variance. A large bandwidth leads to a very smooth (i.e. high-bias) density distribution. A small bandwidth leads to an unsmooth (i.e. high-variance) density distribution. The unit of measurement is free. It must be the same as the unit of \code{cell_size} variable. 
  
  (Rayon de lissage de l estimation d intensite par noyau. Cette bande-passante se comporte comme un parametre de lissage, controlant l equilibre entre biais et variance. Un rayon eleve conduit a une densite tres lissee, avec un biais eleve. Un petit rayon genere une densite peu lisse avec une forte variance. Le choix de l unite de mesure est laisse libre a l utilisateur. Elle doit seulement etre la meme que celle de la variable \code{cell_size}.}
  
  \item{list_var}{List of variables (Liste des variables)}
  
  \item{neighbor}{Technical parameter, leave empty. (Parametre technique pour calculer l etendue des points d estimations, a ne pas remplir.)}
 }

\details{
Returns an object inheriting from the \code{data.frame} class. (Retourne un objet qui se comporte comme un \code{data.frame}, par heritage.)

\itemize{
\item{Smoothing covers a set of methods to extract pertinent and structuring information from noisy data.  In the field of spatial analysis, and most widely in quantiative geography, smoothing is used to modelise density variations of a population distribution in geographical space. Kernel smootthing methods are widely used.
In this method, for each location x, we count the number of events of a process within a distance h of x, and weighted by the square reciprocal of the radius h. We apply a edge-correction to deal with edge-effects. So the method is conservative.}.

\item{Le lissage recouvre un ensemble de methodes pour extraire d une source de donnees bruitees une information pertinente et structurante. Dans le champ de l analyse spatiale et plus largement de la geographie quantitative, le lissage est principalement utilise pour modeliser les variations de densites d une distribution de population dans l espace geographique. On utilise principalement des methodes de lissage par noyau.
Il s agit ici, pour chaque point x, de comptabliser le nombre d "evenements" d un processus a une distance h de ce point, tout en ponderant ce nombre par l inverse de la distance h au carre. On applique une correction a la ponderation afin de traiter les effets de bord. Cette methode est conservative.}
} 
} 
 
\references{
Statistical Analysis of Spatial and Spatio-Temporal Point Patterns, Third Edition, Diggle, 2003, pp. 83-86
}

\author{Psar Analyse Urbaine Insee - Thierry Cornely, Laure Genebes, Arlindo Dos Santos, Cynthia Faivre, Auriane Renaud and Francois Semecurbe
}

\examples{
\dontrun{
data(reunion)
# Smoothing of the variable houhold (households) and phouhold (poor households) for Reunion
# Lissage de la variable houhold(menages) et phouhold(menages pauvres) pour La Reunion
reunion_smooth=kernel_smoothing(reunion,200,800,c("houhold","phouhold"))

# Calculating the poverty rate (Calcul du taux de pauvrete)
reunion_smooth$ratio=reunion_smooth$phouhold/reunion_smooth$houhold*100

# Building of the associated basemap (Creation du fond de carte associe)
grid=smoothing_to_grid(reunion_smooth,"32740")

#library(rgdal)
# Export of the basemap in shapefile format (Export du fond de carte au format shapefile)
#writeOGR(grid,"reunion.shp","reunion",driver="ESRI Shapefile")
}}