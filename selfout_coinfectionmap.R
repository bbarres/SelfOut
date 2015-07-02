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
spatstat.options(npixel=c(nx=500, ny=500)) #default values are 100,100


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

myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)
#map of the estimate of the intensity function of the sampling process
#coordinates of the sampling points
samppoint<-cbind("x"=geno2012t$long_plant,"y"=geno2012t$lat_plant)
samppoint<-samppoint[complete.cases(samppoint),] #removing points without coordinates
samppointppp <- ppp(samppoint[,1], samppoint[,2], window=boundary) 
kdensity<-density.ppp(samppointppp, 1268,diggle=TRUE)
breakdens<-plot(kdensity,col=myPal)
breakdens<-attr(breakdens,"stuff")$breaks[c(1,5,10,15,20,25)]
breakdens<-as.character(round(round(breakdens*100000000,digits=0)/1000,digits=1))
#computing the relative risk surface of coinfection compare to pure infection
purecoinf<-data.frame("x"=geno2012t$long_plant,"y"=geno2012t$lat_plant,
                      "marks"=as.factor(geno2012t$nb_snp_het))
purecoinf<-purecoinf[complete.cases(purecoinf),] #removing points with missing things
purecoinfppp <- ppp(purecoinf[,1], purecoinf[,2], window=boundary,marks=purecoinf[,3]) 
prisk <- relrisk(purecoinfppp, 1268,at="pixels",diggle=TRUE,case=1)
breakcoin<-plot(prisk,col=myPal)
breakcoin<-attr(breakcoin,"stuff")$breaks[c(1,5,10,15,20,25)]
breakcoin<-as.character(round(breakcoin,digits=1))

#code for the map
op<-par(mfrow=c(1,2))
#map of the estimate of the intensity function of the sampling process
col.labels<-paste(breakdens,"e-02")
myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)
plot(Aland,lty=0)
title(main="Map of the estimate of the intensity function of the sampling process")
plot(kdensity,col=myPal,add=TRUE)
plot(Aland,add=TRUE)
scalebar(c(86000,6667000),20000,"km")
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb")
points(purecoinf[,1:2],pch=19,cex=0.3,col=grey(0.5))
#mapping of the relative risk surface of coinfection compare to pure infection
col.labels<-breakcoin
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)
plot(Aland,lty=0)
title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk,col=myPal,add=TRUE)
plot(Aland,add=TRUE)
points(purecoinf[,1:2],pch=19,cex=0.3,col=grey(0.5))
scalebar(c(86000,6667000),20000,"km")
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb")
par(op) #export the map 32 x 10 inches

#map send on the 150514 for the first version of the draft
col.labels<-breakcoin
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)
plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk,col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#points(purecoinf[,1:2],pch=19,cex=1,col=grey(0.3))
#points(purecoinf[,1:2],pch=19,cex=1.2,col=colours()[81])
points(purecoinf[,1:2],pch=19,cex=1.2,col=grey(0.2))
#scalebar(c(86000,6667000),20000,"km")
scalebar(c(86000,6667000),20000,"km",division.cex=2)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) #export the map 16 x 10 inches or 1600 x 1000 for jpeg

#map send on the 180814 for the revised version of the manuscript
#first changing the scale so that close to 0% aren't white like no-data
col.labels<-breakcoin
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)[4:25]
plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk,col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#points(purecoinf[,1:2],pch=19,cex=1,col=grey(0.3))
#points(purecoinf[,1:2],pch=19,cex=1.2,col=colours()[81])
points(purecoinf[,1:2],pch=19,cex=1.2,col=grey(0.2))
#scalebar(c(86000,6667000),20000,"km")
scalebar(c(86000,6667000),20000,"km",division.cex=2)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) #export the map 16 x 10 inches or 1600 x 1000 for jpeg

#second changing the colour of no-data so that they are grey and then different than 0% coinfection
col.labels<-breakcoin
myPal<-colorRampPalette(brewer.pal(9,"Reds"))(25)
plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk,col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
plot(Aland_sup,add=TRUE,lwd=3, col=grey(0.8))
#points(purecoinf[,1:2],pch=19,cex=1,col=grey(0.3))
#points(purecoinf[,1:2],pch=19,cex=1.2,col=colours()[81])
points(purecoinf[,1:2],pch=19,cex=1.2,col=grey(0.2))
#scalebar(c(86000,6667000),20000,"km")
scalebar(c(86000,6667000),20000,"km",division.cex=2)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) #export the map 16 x 10 inches or 1600 x 1000 for jpeg

#third, changing the colour scale using a divergent scale instead of one gradient color
col.labels<-breakcoin
myPal<- colorRampPalette(brewer.pal(11,"Spectral"))(25)[25:1]
plot(Aland,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk,col=myPal,add=TRUE)
plot(Aland,add=TRUE,lwd=3)
#points(purecoinf[,1:2],pch=19,cex=1,col=grey(0.3))
#points(purecoinf[,1:2],pch=19,cex=1.2,col=colours()[81])
points(purecoinf[,1:2],pch=19,cex=1.2,col=grey(0.2))
#scalebar(c(86000,6667000),20000,"km")
scalebar(c(86000,6667000),20000,"km",division.cex=2)
color.legend(165000,6670000,167000,6715000,col.labels,
             myPal,gradient="y",align="rb",cex=2) #export the map 16 x 10 inches or 1600 x 1000 for jpeg


#Thumbnail image code
Aland_limlim<-Aland[Aland$RECNO %in% c(1,4,5,12,15,21,23,25,29,34,35,39,40,47,52),]
boundary<-as(as(Aland_limlim, "SpatialPolygons"), "owin")
spatstat.options(npixel=c(nx=500, ny=500)) #default values are 100,100
samppoint<-cbind("x"=geno2012t$long_plant,"y"=geno2012t$lat_plant)
samppoint<-samppoint[complete.cases(samppoint),] #removing points without coordinates
samppointppp <- ppp(samppoint[,1], samppoint[,2], window=boundary) 
kdensity<-density.ppp(samppointppp, 1268,diggle=TRUE)
myPal<-colorRampPalette(brewer.pal(9,"Blues"))(25)
breakdens<-plot(kdensity,col=myPal)
breakdens<-attr(breakdens,"stuff")$breaks[c(1,5,10,15,20,25)]
breakdens<-as.character(round(round(breakdens*100000000,digits=0)/1000,digits=1))
#computing the relative risk surface of coinfection compare to pure infection
purecoinf<-data.frame("x"=geno2012t$long_plant,"y"=geno2012t$lat_plant,
                      "marks"=as.factor(geno2012t$nb_snp_het))
purecoinf<-purecoinf[complete.cases(purecoinf),] #removing points with missing things
purecoinfppp <- ppp(purecoinf[,1], purecoinf[,2], window=boundary,marks=purecoinf[,3]) 
prisk <- relrisk(purecoinfppp, 1268,at="pixels",diggle=TRUE,case=1)
breakcoin<-plot(prisk,col=myPal)
breakcoin<-attr(breakcoin,"stuff")$breaks[c(1,5,10,15,20,25)]
breakcoin<-as.character(round(breakcoin,digits=1))
col.labels<-breakcoin
myPal<- colorRampPalette(brewer.pal(11,"Spectral"))(25)[25:1]
plot(Aland_limlim,lty=0)
#title(main="Map of the relative risk surface of coinfection vs infection")
plot(prisk,col=myPal,add=TRUE)
plot(Aland_limlim,add=TRUE,lwd=3)

