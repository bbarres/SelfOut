###############################################################################
###############################################################################
#barplot of the effect of coinfection on different parameters
###############################################################################
###############################################################################


op<-par(mfrow=c(1,5),mar=c(5.1,5.1,4.1,2.1), oma=c(0,0,0,0))

coolcol<-c("dodgerblue3","deeppink")
temp<-barplot(matrix(c(14,16,23,56),nrow=2),beside=TRUE,border=NA,col=coolcol,
              lwd=5,ylim=c(0,70),names.arg=list("2012/2013","2013/2014"),
              main=coolcol,font=2,las=1,cex.lab=1.5,
              ylab="Percentage of overwintering success")

coolcol<-c("darkblue","deeppink2")
temp<-barplot(matrix(c(14,16,23,56),nrow=2),beside=TRUE,border=NA,col=coolcol,
              lwd=5,ylim=c(0,70),names.arg=list("2012/2013","2013/2014"),
              main=coolcol,font=2,las=1,cex.lab=1.5,
              ylab="Percentage of overwintering success")

coolcol<-c("dodgerblue3","deeppink2")
temp<-barplot(matrix(c(14,16,23,56),nrow=2),beside=TRUE,border=NA,col=coolcol,
              lwd=5,ylim=c(0,70),names.arg=list("2012/2013","2013/2014"),
              main=coolcol,font=2,las=1,cex.lab=1.5,
              ylab="Percentage of overwintering success")

coolcol<-c("darkblue","pink2")
temp<-barplot(matrix(c(14,16,23,56),nrow=2),beside=TRUE,border=NA,col=coolcol,
              lwd=5,ylim=c(0,70),names.arg=list("2012/2013","2013/2014"),
              main=coolcol,font=2,las=1,cex.lab=1.5,
              ylab="Percentage of overwintering success")

coolcol<-c("royalblue","hotpink")
temp<-barplot(matrix(c(14,16,23,56),nrow=2),beside=TRUE,border=NA,col=coolcol,
              lwd=5,ylim=c(0,70),names.arg=list("2012/2013","2013/2014"),
              main=coolcol,font=2,las=1,cex.lab=1.5,
              ylab="Percentage of overwintering success")

text(temp[1],matrix(c(14,16,23,56),nrow=2)[1]+3,"n=214",font=2)
text(temp[2],matrix(c(14,16,23,56),nrow=2)[2]+3,"n=312",font=2)

text((temp[1]+temp[2])/2,max(matrix(c(14,16,23,56),nrow=2)[1:2])+10,"*",font=2,cex=4)

text(temp[3],matrix(c(14,16,23,56)+3,nrow=2)[3],"n=682",font=2)
text(temp[4],matrix(c(14,16,23,56)+3,nrow=2)[4],"n=454",font=2)

text((temp[3]+temp[4])/2,max(matrix(c(14,16,23,56),nrow=2)[3:4])+10,"***",font=2,cex=4)

par(op)


which(colors()=="hotpink")







