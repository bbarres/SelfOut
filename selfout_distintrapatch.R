###############################################################################
###############################################################################
#Mean distance between samples inside patches
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")

#the aim here is to evaluate the mean distance between genetic samples at the 
#intra-patch level

library(gdata)

###############################################################################
#formating the data
###############################################################################

intradist10<-geno_hom[!is.na(geno_hom$leave_ID) & geno_hom$survey=="MS10" & 
                      geno_hom$nb_missing<1 & geno_hom$nb_snp_typed>18 & 
                      geno_hom$duplicate!="DUP" & geno_hom$leave_ID=="a",]

intradist11<-geno_hom[!is.na(geno_hom$leave_ID) & geno_hom$survey=="MS11" & 
                      geno_hom$nb_missing<1 & geno_hom$nb_snp_typed>18 & 
                      geno_hom$duplicate!="DUP" & geno_hom$leave_ID=="a",]

intradist12<-geno_hom[!is.na(geno_hom$leave_ID) & geno_hom$survey=="MS12" & 
                      geno_hom$nb_missing<1 & geno_hom$nb_snp_typed>18 & 
                      geno_hom$duplicate!="DUP" & geno_hom$leave_ID=="a",]

intradist13<-geno_hom[!is.na(geno_hom$leave_ID) & geno_hom$survey=="MS13" & 
                      geno_hom$nb_missing<1 & geno_hom$nb_snp_typed>18 & 
                      geno_hom$duplicate!="DUP" & geno_hom$leave_ID=="a",]


###############################################################################
#define a function that returns the mean and standard deviation of 
###############################################################################

CompIntraDist<-function(intradist){
  #we list the patches ID in the dataset that have at least 2 samples
  samppatch<-summary(drop.levels(intradist$patche_ID),maxsum=10000)
  samppatch<-as.numeric(attr(samppatch[samppatch>1],"names"))
  #we build two lists, the first is the list of all pairwise distances between 
  #samples within patch for every patch, the second is the list of the max 
  #distance between two samples within a patch for every patch
  listintradist<-c()
  listmaxintradist<-c()
  for (i in 1:length(samppatch)){
    listintradist<-c(listintradist,
                     dist(intradist[intradist$patche_ID==samppatch[i],
                                    c("long_plant","lat_plant")]))
    listmaxintradist<-rbind(listmaxintradist,c(samppatch[i],
                       max(dist(intradist[intradist$patche_ID==samppatch[i],
                                          c("long_plant","lat_plant")]),
                           na.rm=TRUE)))
    
  }
  #we compute the mean and the standard deviation of the intra patch distances
  MeInPaDi<-mean(listintradist,na.rm=TRUE)
  SDInPaDi<-sqrt(var(listintradist,na.rm=TRUE))
  #we group the different results in a list
  rez<-list("Mean intra patch distance (m)"=MeInPaDi,
            "Standard deviation of intra patch distance"=SDInPaDi,
            "List of intra patch distances"=listintradist,
            "List of maximal intra patch distances"=listmaxintradist)
  return(rez)
}


###############################################################################
#computation for each year
###############################################################################

op<-par(mfrow=c(4,2),mar=c(5,5,3,2))

temp<-CompIntraDist(intradist10)
plot(temp$`List of maximal intra patch distances`[,2], 
     main="Distance between samples within patch in 2010",
     ylab="Distance (m)",cex.main=0.7)
plot(density(temp$`List of intra patch distances`,na.rm=TRUE),
     xlim=c(0,100),main=paste(sep="","Mean distance in 2010 = ",
     round(temp$`Mean intra patch distance (m)`,2), "m (+/-",
     round(temp$`Standard deviation of intra patch distance`,2),")"),
     xlab="Distance (m)",cex.main=0.7)
abline(v=temp$`Mean intra patch distance (m)`,col="red")

temp<-CompIntraDist(intradist11)
plot(temp$`List of maximal intra patch distances`[,2], 
     main="Distance between samples within patch in 2011",
     ylab="Distance (m)",cex.main=0.7)
plot(density(temp$`List of intra patch distances`,na.rm=TRUE),
     xlim=c(0,100),main=paste(sep="","Mean distance in 2011 = ",
     round(temp$`Mean intra patch distance (m)`,2), "m (+/-",
     round(temp$`Standard deviation of intra patch distance`,2),")"),
     xlab="Distance (m)",cex.main=0.7)
abline(v=temp$`Mean intra patch distance (m)`,col="red")

temp<-CompIntraDist(intradist12)
plot(temp$`List of maximal intra patch distances`[,2], 
     main="Distance between samples within patch in 2012",
     ylab="Distance (m)",cex.main=0.7)
plot(density(temp$`List of intra patch distances`,na.rm=TRUE),
     xlim=c(0,100),main=paste(sep="","Mean distance in 2012 = ",
     round(temp$`Mean intra patch distance (m)`,2), "m (+/-",
     round(temp$`Standard deviation of intra patch distance`,2),")"),
     xlab="Distance (m)",cex.main=0.7)
abline(v=temp$`Mean intra patch distance (m)`,col="red")

temp<-CompIntraDist(intradist13)
plot(temp$`List of maximal intra patch distances`[,2], 
     main="Distance between samples within patch in 2013",
     ylab="Distance (m)",cex.main=0.7)
plot(density(temp$`List of intra patch distances`,na.rm=TRUE),
     xlim=c(0,100),main=paste(sep="","Mean distance in 2013 = ",
     round(temp$`Mean intra patch distance (m)`,2), "m (+/-",
     round(temp$`Standard deviation of intra patch distance`,2),")"),
     xlab="Distance (m)",cex.main=0.7)
abline(v=temp$`Mean intra patch distance (m)`,col="red")

par(op)

#cleanning the environment
rm(intradist10,intradist11,intradist12,intradist13,temp,CompIntraDist)

#export to pdf file 8 X 16 inches


###############################################################################
#END
###############################################################################