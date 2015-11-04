###############################################################################
###############################################################################
#map of the kernel density of coinfection
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")

library(spatstat)
library(RColorBrewer)
library(plotrix)
library(maptools)


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

#select a subset of the genotype according to their sampling year and genotype 
#quality: first 2012
geno2012t<-geno_hom[!is.na(geno_hom$leave_ID) & geno_hom$survey=="MS12" & 
                    geno_hom$nb_missing<1 & geno_hom$nb_snp_typed>18 & 
                    geno_hom$duplicate!="DUP" & geno_hom$leave_ID=="a",]
geno2012t[(geno2012t$nb_snp_het>0),"nb_snp_het"]<-"COINF"
geno2012t[(geno2012t$nb_snp_het==0),"nb_snp_het"]<-"PURE"

#select a subset of the genotype according to their sampling year and genotype 
#quality: then 2013
geno2013t<-geno_hom[!is.na(geno_hom$leave_ID) & geno_hom$survey=="MS13" & 
                    geno_hom$nb_missing<1 & geno_hom$nb_snp_typed>18 & 
                    geno_hom$duplicate!="DUP" & geno_hom$leave_ID=="a",]
geno2013t[(geno2013t$nb_snp_het>0),"nb_snp_het"]<-"COINF"
geno2013t[(geno2013t$nb_snp_het==0),"nb_snp_het"]<-"PURE"


#change the default value of the resolution of 'spatstat' package
spatstat.options(npixel=c(nx=750, ny=750)) #default values are 100,100


###############################################################################
#Plotting maps of samples
###############################################################################

op<-par(mfrow=c(1,2))

#mapping coinfected and not coinfected samples in 2012
plot(boundary,
     main="Coinfected (red) and single infected (green) samples in 2012")
#plotting of the samples where a no coinfection was detected
points(geno2012t[geno2012t$nb_snp_het=="PURE",]$long_plant,
       geno2012t[geno2012t$nb_snp_het=="PURE",]$lat_plant,
       cex=2.5,bg=rgb(0,0.7,0.4,0.2),pch=21,col="transparent")
#plotting of the samples where a coinfection was detected
points(geno2012t[geno2012t$nb_snp_het=="COINF",]$long_plant,
       geno2012t[geno2012t$nb_snp_het=="COINF",]$lat_plant,
       cex=2.5,bg=rgb(1,0,0,0.2),pch=21,col="transparent")

#mapping coinfected and not coinfected samples in 2013
plot(boundary,
     main="Coinfected (red) and single infected (green) samples in 2013")
#plotting of the samples where a no coinfection was detected
points(geno2013t[geno2013t$nb_snp_het=="PURE",]$long_plant,
       geno2013t[geno2013t$nb_snp_het=="PURE",]$lat_plant,
       cex=2.5,bg=rgb(0,0.7,0.4,0.2),pch=21,col="transparent")
#plotting of the samples where a coinfection was detected
points(geno2013t[geno2013t$nb_snp_het=="COINF",]$long_plant,
       geno2013t[geno2013t$nb_snp_het=="COINF",]$lat_plant,
       cex=2.5,bg=rgb(1,0,0,0.2),pch=21,col="transparent")

par(op)


###############################################################################
#computing and mapping the relative risk surface of coinfection
###############################################################################

#Maps for 2012

myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)

#map of the estimate of the intensity function of the sampling process
#coordinates of the sampling points
samppoint<-cbind("x"=geno2012t$long_plant,"y"=geno2012t$lat_plant)
#removing points without coordinates
samppoint<-samppoint[complete.cases(samppoint),]

samppoint<-ppp(samppoint[,1], samppoint[,2], window=boundary) 
kdensityppp<-density.ppp(samppoint, 1157,diggle=TRUE)
breakdens<-plot(kdensityppp,col=myPal)
breakdens<-attr(breakdens,"stuff")$breaks[c(1,5,10,15,20,25)]
breakdens<-as.character(round(round(breakdens*100000000,digits=0)/1000,
                              digits=1))
#computing the relative risk surface of coinfection compare to pure infection
purecoinf2012<-data.frame("x"=geno2012t$long_plant,"y"=geno2012t$lat_plant,
                      "marks"=as.factor(geno2012t$nb_snp_het))
#removing points with missing things
purecoinf2012<-purecoinf2012[complete.cases(purecoinf2012),]
purecoinfppp2012<-ppp(purecoinf2012[,1], purecoinf2012[,2],
                      window=boundary,marks=purecoinf2012[,3]) 
prisk2012<-relrisk(purecoinfppp2012,1157,at="pixels",diggle=TRUE,case=1)
breakcoin2012<-plot(prisk2012,col=myPal)
breakcoin2012<-attr(breakcoin2012,"stuff")$breaks[c(1,5,10,15,20,25)]
breakcoin2012<-as.character(round(breakcoin2012,digits=1))

#code for the map
op<-par(mfrow=c(1,2))
#map of the estimate of the intensity function of the sampling process
col.labels<-paste(breakdens,"e-02")
myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)
plot(Aland,lty=0)
title(main="Map of the estimate of the intensity function of the sampling 
      process")
plot(kdensityppp,col=myPal,add=TRUE)
plot(Aland,add=TRUE)
scalebar(c(86000,6667000),20000,"km")
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb")
points(purecoinf2012[,1:2],pch=19,cex=0.3,col=grey(0.5))
#mapping of the relative risk surface of coinfection compare to pure infection
col.labels<-breakcoin2012
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)
plot(Aland,lty=0)
title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk2012,col=myPal,add=TRUE)
plot(Aland,add=TRUE)
points(purecoinf2012[,1:2],pch=19,cex=0.3,col=grey(0.5))
scalebar(c(86000,6667000),20000,"km")
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb")
par(op) #export the map 32 x 10 inches

#map with the color scale not starting as white (so we can distinguish between 
#no data and no coinfection)
col.labels<-breakcoin2012
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)[3:25]
plot(Aland,lty=0)
title(main="Map of the relative risk surface of coinfection vs infection 
      in 2012")
plot(prisk2012,col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(purecoinf2012[,1:2],pch=19,cex=1.2,col=grey(0.2))
#adding the scalebar
scalebar(c(86000,6667000),20000,"km",division.cex=1.5)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) 
#export the map 16 x 10 inches or 1600 x 1000 for jpeg


#Maps for 2013

myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)

#map of the estimate of the intensity function of the sampling process
#coordinates of the sampling points
samppoint<-cbind("x"=geno2013t$long_plant,"y"=geno2013t$lat_plant)
#removing points without coordinates
samppoint<-samppoint[complete.cases(samppoint),]

samppoint<-ppp(samppoint[,1], samppoint[,2], window=boundary) 
kdensityppp<-density.ppp(samppoint, 1251,diggle=TRUE)
breakdens<-plot(kdensityppp,col=myPal)
breakdens<-attr(breakdens,"stuff")$breaks[c(1,5,10,15,20,25)]
breakdens<-as.character(round(round(breakdens*100000000,digits=0)/1000,
                              digits=1))
#computing the relative risk surface of coinfection compare to pure infection
purecoinf2013<-data.frame("x"=geno2013t$long_plant,"y"=geno2013t$lat_plant,
                          "marks"=as.factor(geno2013t$nb_snp_het))
#removing points with missing things
purecoinf2013<-purecoinf2013[complete.cases(purecoinf2013),]
purecoinfppp2013<-ppp(purecoinf2013[,1], purecoinf2013[,2],
                      window=boundary,marks=purecoinf2013[,3]) 
prisk2013<-relrisk(purecoinfppp2013,1251,at="pixels",diggle=TRUE,case=1)
breakcoin2013<-plot(prisk2013,col=myPal)
breakcoin2013<-attr(breakcoin2013,"stuff")$breaks[c(1,5,10,15,20,25)]
breakcoin2013<-as.character(round(breakcoin2013,digits=1))

#code for the map
op<-par(mfrow=c(1,2))
#map of the estimate of the intensity function of the sampling process
col.labels<-paste(breakdens,"e-02")
myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)
plot(Aland,lty=0)
title(main=
      "Map of the estimate of the intensity function of the sampling process")
plot(kdensityppp,col=myPal,add=TRUE)
plot(Aland,add=TRUE)
scalebar(c(86000,6667000),20000,"km")
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb")
points(purecoinf2013[,1:2],pch=19,cex=0.3,col=grey(0.5))
#mapping of the relative risk surface of coinfection compare to pure infection
col.labels<-breakcoin2013
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)
plot(Aland,lty=0)
title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk2013,col=myPal,add=TRUE)
plot(Aland,add=TRUE)
points(purecoinf2013[,1:2],pch=19,cex=0.3,col=grey(0.5))
scalebar(c(86000,6667000),20000,"km")
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb")
par(op)
#export the map 32 x 10 inches

#map with the color scale not starting as white (so we can distinguish between 
#no data and no coinfection)
col.labels<-breakcoin2013
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)[3:25]
plot(Aland,lty=0)
title(main="Map of the relative risk surface of coinfection vs infection 
      in 2013")
plot(prisk2013,col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(purecoinf2013[,1:2],pch=19,cex=1.2,col=grey(0.2))
#adding the scalebar
scalebar(c(86000,6667000),20000,"km",division.cex=1.5)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) 
#export the map 16 x 10 inches or 1600 x 1000 for jpeg


###############################################################################
#Map of the two relative risk surface in 2012 and 2013
###############################################################################

op<-par(mfrow=c(1,2))

#map with the color scale not starting as white (so we can distinguish between 
#no data and no coinfection)
col.labels<-c("0","0.2","0.4","0.6","0.8")
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)[3:25]

plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection 
#      in 2012")
plot(prisk2012,col=myPal,zlim=c(0,0.8),add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(purecoinf2012[,1:2],pch=19,cex=1.2,col=grey(0.2))
#adding the scalebar
scalebar(c(86000,6667000),20000,"km",division.cex=1.5)
#color.legend(165000,6670000,167000,6715000,col.labels,
#             myPal,gradient="y",align="rb",cex=2) 

plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection 
#      in 2013")
plot(prisk2013,col=myPal,zlim=c(0,0.8),add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(purecoinf2013[,1:2],pch=19,cex=1.2,col=grey(0.2))
#adding the scalebar
scalebar(c(86000,6667000),20000,"km",division.cex=1.5)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) 

par(op)


###############################################################################
#Map for the final figure of the manuscript
###############################################################################

#map with the color scale not starting as white (so we can distinguish between 
#no data and no coinfection)
col.labels<-c("0","0.2","0.4","0.6","0.8")
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)[3:25]

plot(Aland,lty=0)
plot(prisk2012,col=myPal,zlim=c(0,0.8),add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(purecoinf2012[,1:2],pch=19,cex=1.2,col=grey(0.2))
#adding the scalebar
scalebar(c(86000,6668000),20000,"km",division.cex=1.5)

#export as a pdf file 12 x 8 inches

plot(Aland,lty=0)
plot(prisk2013,col=myPal,zlim=c(0,0.8),add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#adding the sampling point
points(purecoinf2013[,1:2],pch=19,cex=1.2,col=grey(0.2))
#adding the scalebar
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) 

#export as a pdf file 12 x 8 inches

#import the two files in gimp and convert them to .tiff after removing alpha, 
#flatten and resize them to 800 x 556 for 2012 map and 931 x 556 for 2013 map. 
#then combine them in powerpoint, add the letter for the 2 panels and convert 
#the file to a pdf that you will turn again to a .tiff of 2049 x 700 pixels


###############################################################################
#Cleaning the environment
###############################################################################

rm(Aland_sup,Aland_lim,myPal,samppoint,purecoinfppp2012,purecoinfppp2013,
   kdensityppp,breakdens,breakcoin2012,breakcoin2013,col.labels,boundary,
   purecoinf2012,purecoinf2013,prisk2012,prisk2013)


###############################################################################
#END
###############################################################################