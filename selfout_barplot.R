###############################################################################
###############################################################################
#barplot of the effect of coinfection on different parameters
###############################################################################
###############################################################################


###############################################################################
#first try
###############################################################################

op<-par(mfrow=c(1,3),mar=c(5.1,16,6.1,2.1), oma=c(0,0,0,0))

#first plot for the overwintering success
coolcol<-c("darkblue","deeppink3")
datmat<-matrix(c(51,74,48,65),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("Winter 2012/2013","Winter 2013/2014"),
              main=coolcol,font=2,las=1,cex.lab=2,lwd=5,ylim=c(0,110))
mtext(side=2,text="Percentage of overwintering success",line=4,font=2,cex=1.2)
text(temp[1],datmat[1]+3,"n=299",font=2)
text(temp[2],datmat[2]+3,"n=338",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+10,"***",font=2,cex=4)
text(temp[3],datmat[3]+3,"n=372",font=2)
text(temp[4],datmat[4]+3,"n=346",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+10,"*",font=2,cex=4)
mtext(side=2,text="A)",font=2,cex=1.2,adj=3,padj=-18,las=1)

#second plot for the colonization success
coolcol<-c("dodgerblue3","deeppink")
datmat<-matrix(c(13,14,12,17),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("2012","2013"),
              main=coolcol,font=2,las=1,cex.lab=2,lwd=5,ylim=c(0,110))
mtext(side=2,text="Percentage of colonyzed populations",line=4,font=2,cex=1.2)
text(temp[1],datmat[1]+3,"n=1263",font=2)
text(temp[2],datmat[2]+3,"n=2044",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+10,"*",font=2,cex=4)
text(temp[3],datmat[3]+3,"n=962",font=2)
text(temp[4],datmat[4]+3,"n=1184",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+10,"*",font=2,cex=4)
mtext(side=2,text="B)",font=2,cex=1.2,adj=3,padj=-18,las=1)

#third plot for the production of overwintering structures
coolcol<-c("royalblue","hotpink")
datmat<-matrix(c(98,98,95,97),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("2012","2013"),
              main=coolcol,font=2,las=1,cex.lab=2,lwd=5,ylim=c(0,110))
mtext(side=2,text="Percentage of populations with resting structures",
      line=4,font=2,cex=1.2)
text(temp[1],datmat[1]+3,"n=274",font=2)
text(temp[2],datmat[2]+3,"n=315",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+10,"ns",font=3,cex=2)
text(temp[3],datmat[3]+3,"n=374",font=2)
text(temp[4],datmat[4]+3,"n=343",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+10,"ns",font=3,cex=2)
mtext(side=2,text="C)",font=2,cex=1.2,adj=3,padj=-18,las=1)

par(op)

#export to pdf 16 X 9 inches


###############################################################################
#second try
###############################################################################

op<-par(mfrow=c(1,3),mar=c(5.1,16,6.1,2.1), oma=c(0,0,0,0))

#first plot for the overwintering success
coolcol<-c("blue","red")
datmat<-matrix(c(51,74,48,65),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("Winter 2012/2013","Winter 2013/2014"),
              main=coolcol,font=2,las=1,cex.lab=2,lwd=5,ylim=c(0,90))
mtext(side=2,text="% of overwintering success",line=4,font=2,cex=1.2)
text(temp[1],datmat[1]+2,"n=299",font=2)
text(temp[2],datmat[2]+2,"n=338",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+6,"***",font=2,cex=4)
text(temp[3],datmat[3]+2,"n=372",font=2)
text(temp[4],datmat[4]+2,"n=346",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+6,"*",font=2,cex=4)
mtext(side=2,text="A)",font=2,cex=1.2,adj=3,padj=-18,las=1)

#second plot for the colonization success
coolcol<-c("midnightblue","magenta")
datmat<-matrix(c(13,14,12,17),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("2012","2013"),
              main=coolcol,font=2,las=1,cex.lab=2,lwd=5,ylim=c(0,25))
mtext(side=2,text="% of colonyzed populations",line=4,font=2,cex=1.2)
text(temp[1],datmat[1]+1,"n=1263",font=2)
text(temp[2],datmat[2]+1,"n=2044",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+3,"*",font=2,cex=4)
text(temp[3],datmat[3]+1,"n=962",font=2)
text(temp[4],datmat[4]+1,"n=1184",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+3,"*",font=2,cex=4)
mtext(side=2,text="B)",font=2,cex=1.2,adj=3,padj=-18,las=1)

#third plot for the production of overwintering structures
coolcol<-c("steelblue2","pink3")
datmat<-matrix(c(98,98,95,97),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("2012","2013"),
              main=coolcol,font=2,las=1,cex.lab=2,lwd=5,ylim=c(0,110))
mtext(side=2,text="% of populations with resting structures",
      line=4,font=2,cex=1.2)
text(temp[1],datmat[1]+3,"n=274",font=2)
text(temp[2],datmat[2]+3,"n=315",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+10,"ns",font=3,cex=2)
text(temp[3],datmat[3]+3,"n=374",font=2)
text(temp[4],datmat[4]+3,"n=343",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+10,"ns",font=3,cex=2)
mtext(side=2,text="C)",font=2,cex=1.2,adj=3,padj=-18,las=1)

par(op)

#export to pdf 16 X 9 inches

which(colors()=="hotpink")


###############################################################################
#third try, this time not in color, but with 2 shades of greys
###############################################################################

op<-par(mfrow=c(1,3),mar=c(5.1,16,6.1,2.1), oma=c(0,0,0,0))

coolcol<-c("grey70","grey30")

#first plot for the colonization success
datmat<-matrix(c(13,14,12,17),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("2012","2013"),
              font=2,las=1,cex.axis=2,cex.names=1,lwd=5,ylim=c(0,22.5))
mtext(side=2,text="% of colonyzed populations",line=6,font=2,cex=1.2)
text(temp[1],datmat[1]+1,"n=1263",font=2)
text(temp[2],datmat[2]+1,"n=2044",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+2.5,"*",font=2,cex=4)
text(temp[3],datmat[3]+1,"n=962",font=2)
text(temp[4],datmat[4]+1,"n=1184",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+2.5,"*",font=2,cex=4)
mtext(side=2,text="A)",font=2,cex=1.2,adj=5,padj=-20,las=1)


#second plot for the overwintering success
datmat<-matrix(c(51,74,48,65),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("Winter 2012/2013","Winter 2013/2014"),
              font=2,las=1,cex.axis=2,cex.names=1,lwd=5,ylim=c(0,91))
mtext(side=2,text="% of overwintering success",line=6,font=2,cex=1.2)
text(temp[1],datmat[1]+2,"n=299",font=2)
text(temp[2],datmat[2]+2,"n=338",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+6,"***",font=2,cex=4)
text(temp[3],datmat[3]+2,"n=372",font=2)
text(temp[4],datmat[4]+2,"n=346",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+6,"*",font=2,cex=4)
mtext(side=2,text="B)",font=2,cex=1.2,adj=5,padj=-20,las=1)


#third plot for the production of overwintering structures
datmat<-matrix(c(98,98,95,97),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("2012","2013"),
              font=2,las=1,cex.axis=2,cex.names=1,lwd=5,ylim=c(0,111))
mtext(side=2,text="% of populations with resting structures",
      line=6,font=2,cex=1.2)
text(temp[1],datmat[1]+3,"n=274",font=2)
text(temp[2],datmat[2]+3,"n=315",font=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+10,"ns",font=3,cex=2)
text(temp[3],datmat[3]+3,"n=374",font=2)
text(temp[4],datmat[4]+3,"n=343",font=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+10,"ns",font=3,cex=2)
mtext(side=2,text="C)",font=2,cex=1.2,adj=5,padj=-20,las=1)

par(op)

#export to pdf 16 X 9 inches


###############################################################################
#fourth try, still not in color, with bigger text font
###############################################################################

op<-par(mfrow=c(1,3),mar=c(5.1,16,3.1,2.1), oma=c(4,0,1,0))

coolcol<-c("grey70","grey30")

#first plot for the colonization success
datmat<-matrix(c(13,14,12,17),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,font=2,las=1,
              cex.axis=3,cex.names=1,lwd=5,ylim=c(0,22.5),xpd=FALSE)
axis(1,at=c(0,6),lwd=5,labels=FALSE,lwd.ticks=0)
axis(1,at=c(2,5),labels=list("2012","2013"),lwd=5,font=2,cex.axis=2.5,padj=1)
mtext(side=2,text="% of colonyzed populations",line=7,font=2,cex=2)
text(temp[1],datmat[1]+1,"n=1263",font=3,cex=2)
text(temp[2],datmat[2]+1,"n=2044",font=3,cex=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+2.5,"*",font=2,cex=4)
text(temp[3],datmat[3]+1,"n=962",font=3,cex=2)
text(temp[4],datmat[4]+1,"n=1184",font=3,cex=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+2.5,"*",font=2,cex=4)
mtext(side=2,text="A)",font=2,cex=2,adj=5,padj=-10,las=1)


#second plot for the overwintering success
datmat<-matrix(c(51,74,48,65),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,font=2,las=1,
              cex.axis=3,cex.names=1,lwd=5,ylim=c(0,91),xpd=FALSE)
axis(1,at=c(0,6),lwd=5,labels=FALSE,lwd.ticks=0)
axis(1,at=c(2,5),labels=list("Winter\n2012/2013","Winter\n2013/2014"),lwd=5,
     font=2,cex.axis=2.5,padj=1)
mtext(side=2,text="% of overwintering success",line=7,font=2,cex=2)
text(temp[1],datmat[1]+3,"n=299",font=3,cex=2)
text(temp[2],datmat[2]+3,"n=338",font=3,cex=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+8,"***",font=2,cex=4)
text(temp[3],datmat[3]+3,"n=372",font=3,cex=2)
text(temp[4],datmat[4]+3,"n=346",font=3,cex=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+8,"*",font=2,cex=4)
mtext(side=2,text="B)",font=2,cex=2,adj=5,padj=-10,las=1)


#third plot for the production of overwintering structures
datmat<-matrix(c(98,98,95,97),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,font=2,las=1,
              cex.axis=3,cex.names=1,lwd=5,ylim=c(0,111),xpd=FALSE)
axis(1,at=c(0,6),lwd=5,labels=FALSE,lwd.ticks=0)
axis(1,at=c(2,5),labels=list("2012","2013"),lwd=5,font=2,cex.axis=2.5,padj=1)
mtext(side=2,text="% of populations with\nresting structures",
      line=7,font=2,cex=2)
text(temp[1],datmat[1]+4,"n=274",font=3,cex=2)
text(temp[2],datmat[2]+4,"n=315",font=3,cex=2)
text((temp[1]+temp[2])/2,max(datmat[1:2])+12,"ns",font=3,cex=2,xpd=TRUE)
text(temp[3],datmat[3]+4,"n=374",font=3,cex=2)
text(temp[4],datmat[4]+4,"n=343",font=3,cex=2)
text((temp[3]+temp[4])/2,max(datmat[3:4])+12,"ns",font=3,cex=2)
mtext(side=2,text="C)",font=2,cex=2,adj=5,padj=-10,las=1)

par(op)

#export to pdf 24.5 x 7 inches
#export to TIFF 2049 x 600 (but problem of resolution for publication
#72 ppi instead of 300-600 ppi required)



###############################################################################
#barplot for the mono-MLG vs several-MLG MLG survival rate
###############################################################################

op<-par(mfrow=c(1,1),mar=c(5.1,6,6.1,2.1), oma=c(0,0,0,0))

coolcol<-c("dodgerblue3","deeppink")
datmat<-matrix(c(94,63,84,39),nrow=2)
temp<-barplot(datmat,beside=TRUE,border=NA,col=coolcol,
              names.arg=list("winter 2011/2012","winter 2012/2013"),
              main=coolcol,font=2,las=1,cex.lab=2,lwd=5,ylim=c(0,100))
mtext(side=2,text="% of MLG overwintering successfully",line=4,font=2,cex=1.2)
text(temp[1],datmat[1]+3,"n=18",font=2)
text(temp[2],datmat[2]+3,"n=29",font=2)
#text((temp[1]+temp[2])/2,max(datmat[1:2])+3,"*",font=2,cex=4)
text(temp[3],datmat[3]+3,"n=80",font=2)
text(temp[4],datmat[4]+3,"n=236",font=2)
#text((temp[3]+temp[4])/2,max(datmat[3:4])+3,"*",font=2,cex=4)
mtext(side=2,text="D)",font=2,cex=1.2,adj=3,padj=-18,las=1)
mtext(side=1,text="n stands for the number of patches that overwintered 
successfully and that were included in this analysis. In blue the patches 
      with 1 MLG before the winter and in red the patches with more than 
      1 MLG before winter",
      font=2,cex=0.5,adj=NA,padj=1.5,las=0)

par(op)

#export to pdf 6 X 9 inches


###############################################################################
#END
###############################################################################