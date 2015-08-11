###############################################################################
###############################################################################
#map of the genotypic richness
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")

library(spatstat)
library(RColorBrewer)
library(plotrix)
library(maptools)

#change the default value of the resolution of 'spatstat' package
spatstat.options(npixel=c(nx=750, ny=750)) #default values are 100,100

###############################################################################
#Functions to plot the scale and arrow
###############################################################################

#function for a scale, found in "Auxiliary Cartographic Functions in R: 
#North Arrow, Scale Bar, and Label with a Leader Arrow", Tanimura et al 2007, 
#J of Statistical software
#The code has been slightly modified in order to convert the meter in km
scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])/1000
  labels <- append(labels,paste((x[5]-loc[1])/1000,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=.5,cex=division.cex)
}

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
  # checking arguments
  if(missing(loc)) stop("loc is missing")
  if(missing(size)) stop("size is missing")
  # default colors are white and black
  if(missing(cols)) cols <- rep(c("white","black"),8)
  # calculating coordinates of polygons
  radii <- rep(size/c(1,4,2,4),4)
  x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
  y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
  # drawing polygons
  for (i in 1:15) {
    x1 <- c(x[i],x[i+1],loc[1])
    y1 <- c(y[i],y[i+1],loc[2])
    polygon(x1,y1,col=cols[i])
  }
  # drawing the last polygon
  polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
  # drawing letters
  b <- c("E","N","W","S")
  for (i in 0:3) text((size+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                      (size+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                      cex=cex)
}


###############################################################################
#Loading and subsetting the data before analyses
###############################################################################

#split the shapefile between islands of the Aland Archipelago where something 
#is happenning, and islands where no data were collected
Aland_lim<-Aland[Aland$RECNO %in% c(1,4,5,12,15,21,23,25,29,31,32,34,35,36,39,
                                    40,42,43,45,47,49,52),]
Aland_sup<-Aland[!(Aland$RECNO %in% c(1,4,5,12,15,21,23,25,29,31,32,34,35,36,
                                      39,40,42,43,45,47,49,52)),]
#extract the borders of the shapefile so we can plot the border independently
boundary<-as(as(Aland_lim, "SpatialPolygons"), "owin")
plot(boundary)

#select patches for which we got a Genotypic richness (ie with at least 3 
#successfully genotyped samples)
coinf2012t<-coinf2012[!is.na(coinf2012$Gr2012),]
coinf2013t<-coinf2013[!is.na(coinf2013$Gr2013),]


###############################################################################
#computing and mapping the smoothing kernel of the genotypic richness
###############################################################################

#Maps for 2012

#myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)

#preparing the dataset for 'spatstat' package
samppatch2012<-data.frame("x"=coinf2012t$Longitude,"y"=coinf2012t$Latitude,
                          "marks"=as.numeric(coinf2012t$Gr2012))
#removing points with missing things
samppatch2012<-samppatch2012[complete.cases(samppatch2012),]
samppatch2012ppp<-ppp(samppatch2012[,1], samppatch2012[,2],
                      marks=samppatch2012[,3],window=boundary) 
#computing the smoothing kernel of the Genotypic richness
#hmax controls the range of trial values of smoothing bandwidth
smoothGr2012<-(Smooth.ppp(samppatch2012ppp,hmax=8500))
#sigma control the bandwidth of the kernel, here we use the same 
#bandwidth as in the coinfection maps smoothing
smoothGr2012<-(Smooth.ppp(samppatch2012ppp,sigma=1157))

#preparing the dataset for 'spatstat' package
samppatch2013<-data.frame("x"=coinf2013t$Longitude,"y"=coinf2013t$Latitude,
                          "marks"=as.numeric(coinf2013t$Gr2013))
#removing points with missing things
samppatch2013<-samppatch2013[complete.cases(samppatch2013),]
samppatch2013ppp<-ppp(samppatch2013[,1], samppatch2013[,2],
                      marks=samppatch2013[,3],window=boundary) 
#computing the smoothing kernel of the Genotypic richness
#hmax controls the range of trial values of smoothing bandwidth
smoothGr2013<-(Smooth.ppp(samppatch2013ppp,hmax=8500))
#sigma control the bandwidth of the kernel, here we use the same 
#bandwidth as in the coinfection maps smoothing
smoothGr2013<-(Smooth.ppp(samppatch2013ppp,sigma=1251))


###############################################################################
#Map of the two relative risk surface in 2012 and 2013
###############################################################################

myPal<-attr(plot(smoothGr2012,zlim=c(1,3)),"stuff")$outputs

op<-par(mfrow=c(1,2))

#map with the color scale not starting as white (so we can distinguish between 
#no data and no coinfection)
col.labels<-c("1","1.5","2","2.5","3")

plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection in 2012")
plot(smoothGr2012,zlim=c(1,3),col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(coinf2012t[,c("Longitude","Latitude")],pch=19,cex=1.2,col=grey(0.5))
#adding the scalebar
scalebar(c(86000,6667000),20000,"km",division.cex=1.5)
#color.legend(165000,6670000,167000,6715000,col.labels,
#             myPal,gradient="y",align="rb",cex=2) 

plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection in 2013")
plot(smoothGr2013,zlim=c(1,3),col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(coinf2013t[,c("Longitude","Latitude")],pch=19,cex=1.2,col=grey(0.5))
#adding the scalebar
scalebar(c(86000,6667000),20000,"km",division.cex=1.5)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) 

par(op)


###############################################################################
#Cleaning the environment
###############################################################################

rm(Aland_sup,Aland_lim,samppoint,boundary,samppatch2012,samppatch2012ppp,
   samppatch2013,samppatch2013ppp,coinf2012t,coinf2013t,smoothGr2012,
   smoothGr2013,myPal)


###############################################################################
#END
###############################################################################