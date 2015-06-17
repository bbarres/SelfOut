###############################################################################
###############################################################################
#map of the Allelic and Genotypic richness across years
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first

#plotting Ar et Gr on the map for 2010
op<-par(mfrow=c(1,2),mar=c(1,1,3,1))
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Ar2010")
points(coinf2010[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2010[,c("Ar2010")])),
       bg=rgb(0,0,1,0.2),pch=21,col="transparent")
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Gr2010")
points(coinf2010[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2010[,c("Gr2010")])),
       bg=rgb(1,0,0,0.2),pch=21,col="transparent")
par(op)

#plotting Ar et Gr on the map for 2011
op<-par(mfrow=c(1,2),mar=c(1,1,3,1))
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Ar2011")
points(coinf2011[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2011[,c("Ar2011")])),
       bg=rgb(0,0,1,0.2),pch=21,col="transparent")
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Gr2011")
points(coinf2011[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2011[,c("Gr2011")])),
       bg=rgb(1,0,0,0.2),pch=21,col="transparent")
par(op)

#plotting Ar et Gr on the map for 2012
op<-par(mfrow=c(1,2),mar=c(1,1,3,1))
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Ar2012")
points(coinf2012[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2012[,c("Ar2012")])),
       bg=rgb(0,0,1,0.2),pch=21,col="transparent")
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Gr2012")
points(coinf2012[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2012[,c("Gr2012")])),
       bg=rgb(1,0,0,0.2),pch=21,col="transparent")
par(op)

#plotting Ar et Gr on the map for 2013
op<-par(mfrow=c(1,2),mar=c(1,1,3,1))
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Ar2013")
points(coinf2013[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2013[,c("Ar2013")])),
       bg=rgb(0,0,1,0.2),pch=21,col="transparent")
plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Gr2013")
points(coinf2013[,c("Longitude","Latitude")],
       cex=as.numeric(as.character(coinf2013[,c("Gr2013")])),
       bg=rgb(1,0,0,0.2),pch=21,col="transparent")
par(op)

#export each map to a 21 X 7 inches pdf file

###############################################################################
#END
###############################################################################