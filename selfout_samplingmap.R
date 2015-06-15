###############################################################################
###############################################################################
#map of the sampling site across years
###############################################################################
###############################################################################



op<-par(mfrow=c(2,2),mar=c(2,1,2,1))

plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Samples 2010")
points(sample_info[sample_info$survey=="MS10",]$long_plant,
       sample_info[sample_info$survey=="MS10",]$lat_plant,
       cex=1,bg=rgb(0,0,1,0.1),pch=21,col="transparent")
mtext(side=1,text=expression(paste("-385 samples\n-82 populations")),
      font=2,cex=0.7,adj=0.17,padj=-3)

plot(Aland,col=grey(0.7),axes=FALSE,lty=0,main="Samples 2011")
points(sample_info[sample_info$survey=="MS11",]$long_plant,
       sample_info[sample_info$survey=="MS11",]$lat_plant,
       cex=1,bg=rgb(1,0,0,0.1),pch=21,col="transparent")
mtext(side=1,text=expression(paste("-427 samples\n-91 populations")),
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
mtext(side=1,text=expression(paste("-2847 samples\n-722 populations")),
      font=2,cex=0.7,adj=0.17,padj=-3)

par(op)

#export at pdf format, 12 X 8 inches
