###############################################################################
###############################################################################
#barplot of the effect of coinfection on different parameters (Figure 2)
###############################################################################
###############################################################################


###############################################################################
#code for the barplot
###############################################################################

op<-par(mfrow=c(1,3),mar=c(5.1,16,3.1,2.1), oma=c(4,1,1,0))

coolcol<-c("grey70","grey30")

#first plot for the colonization success
datmat<-matrix(c(13,14,12,17),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,font=2,las=1,
              cex.axis=3.5,cex.names=1,lwd=5,ylim=c(0,22.5),xpd=FALSE)
axis(1,at=c(0,6),lwd=5,labels=FALSE,lwd.ticks=0)
axis(1,at=c(2,5),labels=list("2012","2013"),lwd=5,font=2,cex.axis=3,padj=1)
mtext(side=2,text="% of colonized populations",line=7,font=2,cex=2.5)
text(temp[1],datmat[1]+1,"n=1263",font=3,cex=2.5)
text(temp[2],datmat[2]+1,"n=2044",font=3,cex=2.5)
text((temp[1]+temp[2])/2,max(datmat[1:2])+2.5,"*",font=2,cex=5)
text(temp[3],datmat[3]+1,"n=962",font=3,cex=2.5)
text(temp[4],datmat[4]+1,"n=1184",font=3,cex=2.5)
text((temp[3]+temp[4])/2,max(datmat[3:4])+2.5,"*",font=2,cex=5)
mtext(side=2,text="A)",font=2,cex=2.5,adj=5,padj=-11,las=1)

#second plot for the overwintering success
datmat<-matrix(c(51,74,48,65),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,font=2,las=1,
              cex.axis=3.5,cex.names=1,lwd=5,ylim=c(0,91),xpd=FALSE)
axis(1,at=c(0,6),lwd=5,labels=FALSE,lwd.ticks=0)
axis(1,at=c(2,5),labels=list("Winter\n2012/2013","Winter\n2013/2014"),lwd=5,
     font=2,cex.axis=3,padj=1)
mtext(side=2,text="% of overwintering success",line=7,font=2,cex=2.5)
text(temp[1],datmat[1]+3,"n=299",font=3,cex=2.5)
text(temp[2],datmat[2]+3,"n=338",font=3,cex=2.5)
text((temp[1]+temp[2])/2,max(datmat[1:2])+8,"***",font=2,cex=5)
text(temp[3],datmat[3]+3,"n=372",font=3,cex=2.5)
text(temp[4],datmat[4]+3,"n=346",font=3,cex=2.5)
text((temp[3]+temp[4])/2,max(datmat[3:4])+8,"*",font=2,cex=5)
mtext(side=2,text="B)",font=2,cex=2.5,adj=5,padj=-11,las=1)

#third plot for the production of overwintering structures
datmat<-matrix(c(98,98,95,97),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,font=2,las=1,
              cex.axis=3.5,cex.names=1,lwd=5,ylim=c(0,111),xpd=FALSE)
axis(1,at=c(0,6),lwd=5,labels=FALSE,lwd.ticks=0)
axis(1,at=c(2,5),labels=list("2012","2013"),lwd=5,font=2,cex.axis=3,padj=1)
mtext(side=2,text="% of populations with\nresting structures",
      line=7,font=2,cex=2.5)
text(temp[1],datmat[1]+4,"n=274",font=3,cex=2.5)
text(temp[2],datmat[2]+4,"n=315",font=3,cex=2.5)
text((temp[1]+temp[2])/2,max(datmat[1:2])+12,"ns",font=3,cex=2.5,xpd=TRUE)
text(temp[3],datmat[3]+4,"n=374",font=3,cex=2.5)
text(temp[4],datmat[4]+4,"n=343",font=3,cex=2.5)
text((temp[3]+temp[4])/2,max(datmat[3:4])+12,"ns",font=3,cex=2.5)
mtext(side=2,text="C)",font=2,cex=2.5,adj=5.5,padj=-11,las=1)

par(op)

#export to pdf 24.5 x 8 inches


###############################################################################
#END
###############################################################################