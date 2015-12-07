###############################################################################
###############################################################################
#map of the sampling site across years
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first


###############################################################################
#Functions to plot the scale
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
  text(x[c(1,3,5)],y[4],labels=labels,adj=c(0.5,0),cex=division.cex)
}


###############################################################################
#The maps
###############################################################################

op<-par(mfrow=c(2,2),mar=c(1,1,2,1),oma=c(0,0,0,0))

plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Samples 2010")
points(sample_info[sample_info$survey=="MS10",]$long_plant,
       sample_info[sample_info$survey=="MS10",]$lat_plant,
       cex=1,bg=rgb(0,0,1,0.1),pch=21,col="transparent")
mtext(side=1,text=expression(paste("-394 samples\n-89 populations")),
      font=2,cex=0.7,adj=0.17,padj=-3)

plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Samples 2011")
points(sample_info[sample_info$survey=="MS11",]$long_plant,
       sample_info[sample_info$survey=="MS11",]$lat_plant,
       cex=1,bg=rgb(1,0,0,0.1),pch=21,col="transparent")
mtext(side=1,text=expression(paste("-452 samples\n-96 populations")),
      font=2,cex=0.7,adj=0.17,padj=-3)

plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Samples 2012")
points(sample_info[sample_info$survey=="MS12",]$long_plant,
       sample_info[sample_info$survey=="MS12",]$lat_plant,
       cex=1,bg=rgb(0,0.5,0.1,0.1),pch=21,col="transparent")
mtext(side=1,text=expression(paste("-4374 samples\n-641 populations")),
      font=2,cex=0.7,adj=0.17,padj=-3)

plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Samples 2013")
points(sample_info[sample_info$survey=="MS13",]$long_plant,
       sample_info[sample_info$survey=="MS13",]$lat_plant,
       cex=1,bg=rgb(0.5,0,0.5,0.1),pch=21,col="transparent")
mtext(side=1,text=expression(paste("-2848 samples\n-722 populations")),
      font=2,cex=0.7,adj=0.17,padj=-3)
scalebar(c(145000,6665200),20000,"kms",division.cex=0.8,xpd=NA)

par(op)

#export at pdf format, 12 X 8 inches


###############################################################################
#END
###############################################################################