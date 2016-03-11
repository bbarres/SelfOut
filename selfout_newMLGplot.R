###############################################################################
###############################################################################
#plot of the effect of % of coinfection on the number of new MLG following year
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' and 
#'selfout_glm.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")

library(visreg)

###############################################################################
#New MLG as a function of the % of coinfection the previous year (Figure 3)
###############################################################################

op<-par(mar=c(6,8,4,2))
visreg(modnewgeno_1,"percoinf12",scale="response",lwd=10,font.lab=2,font=2,
       cex.lab=2,cex.axis=2,ann=FALSE,nn=10000,
       font.axis=2,xlab="% of co-infection in 2012",line.par=list(col="black"),
       ylab="Number of new MLG in 2013",bty="l",rug=FALSE,xlim=c(0,100),
       axes=FALSE,xpd=FALSE)
box(bty="l",lwd=4)
axis(1,lwd=3,font=2,cex.axis=3,padj=0.4)
axis(2,lwd=3,font=2,las=1,cex.axis=3)
mtext(side=1,text="% of co-infection in 2012",
      line=1,font=2,cex=3,padj=2)
mtext(side=2,text="Number of new MLG in 2013",
      line=1,font=2,cex=3,padj=-2.5)
rug(jitter(datanewgeno$percoinf12,amount=1),col="violetred",pos=0.065)
par(op)
#transparent colors don't seem to work with 'rug'. In case this is fixed at 
#some point, here is an example of transparent color: 
#rgb(208,32,144,alpha=0.3,maxColorValue=255)

#export to a pdf file of 14 x 10 inches


###############################################################################
#An alternative, less fancy plot to display the increase of new MLG with coinf
###############################################################################

breaks<-c(-1,10,20,30,40,50)
freq.cut<-cut(datanewgeno$percoinf12,breaks)
boxplot(datanewgeno$number_pure_new.y~freq.cut)


###############################################################################
#END
###############################################################################