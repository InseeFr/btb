#Carroyage est le nom de la classe definie
#les slots sont des variables typees
setClass("Carroyage",
         slots = list(taille_carreau = "numeric",rayon_de_lissage= "numeric"),
         contains = "data.frame")


#fonction generique = utilisable pour plusieurs classes differentes 
setGeneric(
  #creation de fond de carte dont la projection a le code epsg defini
  name="creation_fond_de_carte",
  def=function(object,epsg){standardGeneric("creation_fond_de_carte")}
)

#on definit plus specifiquement la fonction creation_fond_de_carte pour un objet carroyage
setMethod(
  f="creation_fond_de_carte",
  signature="Carroyage",
  definition= function(object,epsg)
  {
    #fonction de creation de carreaux de cote taille_carreau
    carreau<-function(x,y,taille_carreau) {(Polygons(list(Polygon(cbind(c(0,taille_carreau,taille_carreau,0,0)+(x-(taille_carreau)/2),c(0,0,taille_carreau,taille_carreau,0)+(y-(taille_carreau)/2)))),paste(x,y,sep="_")))}
    
    #on applique la fonction carreau
    #MoreArgs contient les parametres supplementaires de la fonction carreau
    grille=mapply(carreau,object[,"x"],object[,"y"],MoreArgs=list(taille_carreau=object@taille_carreau))
    
    #on cree un spatial polygon avec un code epsg de projection defini
    grille_spat = SpatialPolygons((grille),proj4string=CRS(paste("+init=epsg:",epsg,sep="")))
    df=data.frame(object@.Data);names(df)=names(object)
    data=data.frame(ID=paste(df[,"x"],df[,"y"],sep="_"),df)
    
    #un SpatialPolygonsDataFrame est un SpatialPolygon auquel on attache une table d attributs
    return(SpatialPolygonsDataFrame(grille_spat, data, match.ID = "ID"))
  }
)

setMethod(
  `[`,
  signature=signature(x="Carroyage"),
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

setMethod(
  `[<-`,
  signature=signature(x="Carroyage"),
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


#fonction pour transformer un lissage en un carroyage (fond de carte)

smoothing_to_grid<-function(df,epsg,cell_size=NULL)
{
  if (is.null(cell_size)){cell_size=df@taille_carreau}
  rayon_de_lissage=NULL
  xcol="x";ycol="y"
  if (is.null(rayon_de_lissage)){rayon_de_lissage=numeric(0)}
  car=df[,c(xcol,ycol,(liste_var=names(df)[!(names(df)%in% c(xcol,ycol))]))]
  names(car)=c("x","y",liste_var)
  return(creation_fond_de_carte(new(Class = "Carroyage",car,taille_carreau=cell_size,rayon_de_lissage=rayon_de_lissage),epsg))
}



#constructeur grand public
kernel_smoothing<-function(df,cell_size,bandwith,list_var,neighbor=max(0,ceiling(bandwith/cell_size/2)-1))
{

  coord_estimation=data.frame(x=floor(df$x/cell_size)*cell_size+(cell_size/2),y=floor(df$y/cell_size)*cell_size+(cell_size/2))
  coord_estimation=unique(coord_estimation)
  
  list_coord_estimation=list()
  compteur=1
  for (voisin_x in -neighbor:neighbor)
  {
    for (voisin_y in -neighbor:neighbor)
    {
      compteur=compteur+1
      list_coord_estimation[[compteur]]=coord_estimation
      list_coord_estimation[[compteur]]$x=coord_estimation$x+voisin_x*cell_size #on construit une "fenetre" reduite autour des pts
      list_coord_estimation[[compteur]]$y=coord_estimation$y+voisin_y*cell_size #on construit une "fenetre" reduite autour des pts
    }
  }
  
  coord_estimation=unique(do.call(rbind,list_coord_estimation))
  
  lissage=lissage(coord_estimation$x,coord_estimation$y,df$x,df$y,bandwith,as.matrix(df[,list_var]),ceiling(((bandwith/cell_size)*2+3)^2))
  resultat=cbind(coord_estimation,data.frame(lissage))
  names(resultat)=c(names(coord_estimation),list_var)
  
  
  return(new(Class = "Carroyage",resultat,taille_carreau=cell_size,rayon_de_lissage=bandwith))
}


