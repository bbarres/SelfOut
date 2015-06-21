###############################################################################
###############################################################################
#Intra seasonal colonisation of the patches in the vicinity of focal patches
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")

library(gdata)
library(lme4)
library(visreg)

###############################################################################
#Compute the distances between every patches
###############################################################################

#to compute the distances between every patches saves some computation time, 
#because you only compute the distance matrix once. After that, when a distance 
#matrix between a subset of patches is needed, you just need to select a subset 
#of the complete distance matrix

#load the coordinates of the patches
coord <- patche_info[,c(5:6)]
row.names(coord) <- patche_info[,1]
#compute distances (km) between patches
vecdistan<-dist(coord[1:2],
                method = "euclidean",diag = FALSE,upper = FALSE)/1000
#turn the distance matrix into a dataframe
vecdistan<-data.frame(t(combn(rownames(coord),2)),as.numeric(vecdistan))


###############################################################################
#Identify the closest focal patch to other patches
###############################################################################

#we build a small function that seek for the closest focal patch to other 
#patches

closerList<-function(datatab,spring_PA,fall_PA,all_patch){
  #extract the list of the focal patches (infected in June and Sept)
  foc_patches<-datatab[spring_PA==1 & fall_PA==1,"patche_ID"]
  foc_patches<-levels(drop.levels(foc_patches))
  #extract the list of the colonized patches (not infected 
  #in June but infected in Sept)
  colo_patches<-datatab[spring_PA!=1 | is.na(spring_PA),"patche_ID"]
  colo_patches<-levels(drop.levels(colo_patches))
  #extract distances between focal patches and other patches
  vecdistanlim<-vecdistan[(vecdistan$X1 %in% foc_patches | 
                             vecdistan$X2 %in% foc_patches),]
  #list of patches excluding focal patches
  other_patches<-patche_info[!is.na(all_patch),"ID"]
  other_patches<-setdiff(other_patches,foc_patches)
  #include patches without information
  #other_patches<-setdiff(patche_info[,1],foc_patches) 
  
  closer<-data.frame("patche1_ID"=character(),"patche2_ID"=character(),
                     "dist"=character())
  for (i in 1:length(other_patches)) {
    test<-vecdistanlim[(vecdistanlim$X1==other_patches[i] | 
                          vecdistanlim$X2==other_patches[i]),]
    test2<-test[test[,3]==min(test[,3]),]
    colnames(test2)<-c("patche1_ID","patche2_ID","dist")
    closer<-rbind(closer,test2)
  }
  return(closer)
}

#for 2013####
closer2013<-closerList(coinf2013,coinf2013$PA_S2013,coinf2013$PA_2013,
                       patche_info$PA_2013)

#for 2012####
closer2012<-closerList(coinf2012,coinf2012$PA_S2012,coinf2012$PA_2012,
                       patche_info$PA_2012)



eval(parse(text=paste(temp,temp2,sep="")))

#extract the list of the focal patches (infected in June and Sept)
foc_patches<-levels(drop.levels(coinf2013$patche_ID[coinf2013$PA_S2013==1 & 
                                                      coinf2013$PA_2013==1]))
#extract the list of the colonized patches (not infected 
#in June but infected in Sept)
colo_patches<-coinf2013[coinf2013$PA_S2013!=1 | is.na(coinf2013$PA_S2013),]
colo_patches<-levels(drop.levels(colo_patches$patche_ID))
#extract distances between focal patches and other patches
vecdistanlim<-vecdistan[(vecdistan$X1 %in% foc_patches | 
                           vecdistan$X2 %in% foc_patches),]
#list of patches excluding focal patches
other_patches<-patche_info[!is.na(patche_info$PA_2013),1]
other_patches<-setdiff(other_patches,foc_patches)
#include patches without information
#other_patches<-setdiff(patche_info[,1],foc_patches) 

closer<-data.frame("patche1_ID"=character(),"patche2_ID"=character(),
                   "dist"=character())
for (i in 1:length(other_patches)) {
  test<-vecdistanlim[(vecdistanlim$X1==other_patches[i] | 
                        vecdistanlim$X2==other_patches[i]),]
  test2<-test[test[,3]==min(test[,3]),]
  colnames(test2)<-c("patche1_ID","patche2_ID","dist")
  closer<-rbind(closer,test2)
}
closer2013<-closer


###############################################################################
#Modelisation of the colonization success: distance, coinfection...
###############################################################################

####don´t forget to set the correct foc_patches!!!!#####

#2013####
foc_patches<-levels(drop.levels(coinf2013$patche_ID[coinf2013$PA_S2013==1 & 
                                                      coinf2013$PA_2013==1]))
colo_patches<-coinf2013[coinf2013$PA_S2013!=1 | is.na(coinf2013$PA_S2013),]
colo_patches<-levels(drop.levels(colo_patches$patche_ID))
temp<-cbind(closer2013,
            "foc_patche"=ifelse(as.character(closer2013$patche1_ID) %in% foc_patches,as.character(closer2013$patche1_ID),
                                           as.character(closer2013$patche2_ID)))
temp<-cbind(temp,"prox_patche"=ifelse(as.character(closer2013$patche1_ID) %in% foc_patches,as.character(closer2013$patche2_ID),
                                      as.character(closer2013$patche1_ID)))
temp<-merge(temp,coinf2013,by.x="foc_patche",by.y="patche_ID")
temp<-merge(temp,patche_info,by.x="prox_patche",by.y="ID")
temp<-merge(temp,coinf2013[,1:19],by.x="prox_patche",by.y="patche_ID",all.x=TRUE)
temp<-data.frame(temp,"coinfYN"=temp$number_coinf.x)
temp$coinfYN[(temp$coinfYN)>0]<-1
temp$coinfYN<-as.factor(temp$coinfYN)
temp$road_PA.y<-as.factor(temp$road_PA.y)
temp<-data.frame(temp,"colonized"=ifelse(as.character(temp$prox_patche) %in% colo_patches,1,0))
temp<-data.frame(temp,"Gdiv"=temp$number_MLG.x/temp$number_genotyped.x)
temp<-temp[!is.na(temp$PLM2_Sept2013.y),]
temp<-temp[!is.na(temp$road_PA.y),]
#here there is a particular problem with patches in Kökar. This place probably wasn't surveyed in July or 
#no infected patch was detected, therefore the closest focal patch is more than 20 km away. 
plot(temp$dist,col=((as.numeric(temp$dist)>20)+1))
#here the probability that the "closest" focal patch is the patch of origin is very small, so we remove these 
#observation from the dataset
temp<-temp[temp$dist<20,]
prevalcolo<-glm(colonized~cumulative_sum.x+Gdiv+connec2013.x+coinfYN+AA_S2013.x+log(dist),
                family=binomial,data=temp)
summary(prevalcolo)
prevalcolo<-glm(colonized~PLM2_Sept2013.y+coinfYN+AA_S2013.x+log(dist),
                family=binomial,data=temp)
summary(prevalcolo)
prevalcolo<-glm(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+AA_S2013.x+log(dist),
                family=binomial,data=temp)
summary(prevalcolo)
op<-par(mfrow=c(1,3))
visreg(prevalcolo,"coinfYN",scale="response",overlay=TRUE,
       xlab="No coinf (0)/ coinf (1)",ylab="P(colonize)")
visreg(prevalcolo,"AA_S2013.x",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Abundance Spring",ylab="P(colonize)")
visreg(prevalcolo,"PLM2_Sept2013.y",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Net trapping area",ylab="P(colonize)")
par(op)
mixprevalcolo<-glmer(colonized~cumulative_sum.x+Gdiv+connec2013.x+coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp)
summary(mixprevalcolo)
mixprevalcolo<-glmer(colonized~connec2013.x+coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp)
summary(mixprevalcolo)
mixprevalcolo<-glmer(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp)
summary(mixprevalcolo)
#remove focal patches where nothing happens (not a single colonization event)
active<-levels(drop.levels(temp[temp$colonized==1,]$foc_patch))
temp2<-temp[temp$foc_patch %in% active,]
prevalcolo<-glm(colonized~PLM2_Sept2013.y+coinfYN+AA_S2013.x+log(dist),
                family=binomial,data=temp2)
summary(prevalcolo)
prevalcolo<-glm(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+AA_S2013.x+log(dist),
                family=binomial,data=temp2)
summary(prevalcolo)
op<-par(mfrow=c(1,3))
visreg(prevalcolo,"coinfYN",scale="response",overlay=TRUE,
       xlab="No coinf (0)/ coinf (1)",ylab="P(colonize)")
visreg(prevalcolo,"dist",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Abundance Spring",ylab="P(colonize)")
visreg(prevalcolo,"PLM2_Sept2013.y",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Net trapping area",ylab="P(colonize)")
par(op)

prevalcolo_1<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2013.y+coinfYN+AA_S2013.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_1)
prevalcolo_2<-glm(colonized~Sh.y+road_PA.y+PLM2_Sept2013.y+coinfYN+AA_S2013.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_2)
prevalcolo_3<-glm(colonized~shadow.y+road_PA.y+PLM2_Sept2013.y+coinfYN+AA_S2013.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_3)
prevalcolo_4<-glm(colonized~shadow.y+Sh.y+PLM2_Sept2013.y+coinfYN+AA_S2013.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_4)
prevalcolo_5<-glm(colonized~shadow.y+Sh.y+road_PA.y+coinfYN+AA_S2013.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_5)
prevalcolo_6<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2013.y+AA_S2013.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_6)
prevalcolo_7<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2013.y+coinfYN+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_7)
prevalcolo_8<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2013.y+coinfYN+AA_S2013.x,
                  family=binomial,data=temp2)
summary(prevalcolo_8)
anova(prevalcolo_2,prevalcolo_1,test="Chisq")

mixprevalcolo<-glmer(colonized~cumulative_sum.x+Gdiv+connec2013.x+coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp2)
summary(mixprevalcolo)
mixprevalcolo<-glmer(colonized~PLM2_Sept2013.y+coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp2)
summary(mixprevalcolo)
mixprevalcolo_1<-glmer(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                       family=binomial,data=temp2)
summary(mixprevalcolo_1)
mixprevalcolo_2<-glmer(colonized~coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                       family=binomial,data=temp2)
summary(mixprevalcolo_2)
mixprevalcolo_3<-glmer(colonized~I(sqrt(PLM2_Sept2013.y))+AA_S2013.x+log(dist)+(1|foc_patche),
                       family=binomial,data=temp2)
summary(mixprevalcolo_3)
mixprevalcolo_4<-glmer(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+log(dist)+(1|foc_patche),
                       family=binomial,data=temp2)
summary(mixprevalcolo_4)
mixprevalcolo_5<-glmer(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+AA_S2013.x+(1|foc_patche),
                       family=binomial,data=temp2)
summary(mixprevalcolo_5)
anova(mixprevalcolo_5,mixprevalcolo_1)

plot(temp[temp$colonized==1,]$dist,col=temp[temp$colonized==1,]$coinfYN)
#remove focal patches belonging to SIN where less than 3 colonization events happen
#first, find the SINs with more than 3 colonization events
SINgood<-attr(table(temp$colonized,temp$SIN_86.x)[2,],"names")[table(temp$colonized,temp$SIN_86.x)[2,]>2]
temp3<-temp[temp$SIN_86.x %in% SINgood,]
prevalcolo<-glm(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+AA_S2013.x+log(dist),
                family=binomial,data=temp3)
summary(prevalcolo)
mixprevalcolo<-glmer(colonized~I(sqrt(PLM2_Sept2013.y))+coinfYN+AA_S2013.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp3)
summary(mixprevalcolo)


###############################################################################
#Figure to examplify the main result of the analysis
###############################################################################

#plotting example of patch colonization
#first load data from the year you are interested in
foc_patches<-levels(drop.levels(coinf2013$patche_ID[coinf2013$PA_S2013==1 & coinf2013$PA_2013==1]))
inf_patches<-levels(drop.levels(coinf2013$patche_ID[coinf2013$PA_2013==1]))
ext_patches<-levels(drop.levels(as.factor(patche_info$ID[patche_info$PA_S2013==1 & patche_info$PA_2013==0])))
plot(Aland,col=grey(0.85),lty=0)
plot(patchshape,col="white",lty=1,lwd=0.1,add=TRUE)
plot(patchshape,col="white",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% inf_patches,1],col="lightblue",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% foc_patches,1],col="darkblue",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% ext_patches,1],col="orange",lty=0,add=TRUE)
#then we can emphasize the focal patch with coinfection
foc_patches_coin<-levels(drop.levels(coinf2013$patche_ID[coinf2013$PA_S2013==1 & coinf2013$PA_2013==1 & coinf2013$number_coinf!=0]))
temp<-cbind(closer2013,"foc_patche"=ifelse(as.character(closer2013$patche1_ID) %in% foc_patches,as.character(closer2013$patche1_ID),
                                           as.character(closer2013$patche2_ID)))
temp<-cbind(temp,"prox_patche"=ifelse(as.character(closer2013$patche1_ID) %in% foc_patches,as.character(closer2013$patche2_ID),
                                      as.character(closer2013$patche1_ID)))
temp<-merge(temp,coinf2013,by.x="foc_patche",by.y="patche_ID")
temp<-merge(temp,patche_info,by.x="prox_patche",by.y="ID")
temp<-merge(temp,coinf2013[,1:19],by.x="prox_patche",by.y="patche_ID",all.x=TRUE)
prox_patches_coin<-levels(drop.levels(temp$prox_patche[temp$PA_2013.y==1 & temp$foc_patche %in% foc_patches_coin]))
plot(patchshape[patchshape[[3]] %in% foc_patches_coin,1],col="darkred",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% prox_patches_coin,1],col="red",lty=0,add=TRUE)




foc_patches<-levels(drop.levels(coinf2012$patche_ID[coinf2012$PA_S2012==1 & coinf2012$PA_2012==1]))
inf_patches<-levels(drop.levels(coinf2012$patche_ID[coinf2012$PA_2012==1]))
ext_patches<-levels(drop.levels(as.factor(patche_info$ID[patche_info$PA_S2012==1 & patche_info$PA_2012==0])))
plot(Aland,col=grey(0.85),lty=0)
plot(patchshape,col="white",lty=1,lwd=0.1,add=TRUE)
plot(patchshape[patchshape[[3]] %in% inf_patches,1],col="blue",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% foc_patches,1],col="green",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% ext_patches,1],col="orange",lty=0,add=TRUE)



###############################################################################
#2012
###############################################################################




#extract the list of the focal patches (infected in June and Sept)
foc_patches<-levels(drop.levels(coinf2012$patche_ID[coinf2012$PA_S2012==1 & coinf2012$PA_2012==1]))
#extract the list of the colonized patches (not infected in June but infected in Sept)
colo_patches<-coinf2012[coinf2012$PA_S2012!=1 | is.na(coinf2012$PA_S2012),]
colo_patches<-levels(drop.levels(colo_patches$patche_ID))
#extract distances between focal patches and other patches
vecdistanlim<-vecdistan[(vecdistan$X1 %in% foc_patches | vecdistan$X2 %in% foc_patches),]
#list of patches excluding focal patches
other_patches<-patche_info[!is.na(patche_info$PA_2012),1]
other_patches<-setdiff(other_patches,foc_patches)
#other_patches<-setdiff(patche_info[,1],foc_patches) #include patches without information

closer<-data.frame("patche1_ID"=character(),"patche2_ID"=character(),"dist"=character())
for (i in 1:length(other_patches)) {#
  test<-vecdistanlim[(vecdistanlim$X1==other_patches[i] | vecdistanlim$X2==other_patches[i]),]
  test2<-test[test[,3]==min(test[,3]),]
  colnames(test2)<-c("patche1_ID","patche2_ID","dist")
  closer<-rbind(closer,test2)
}
closer2012<-closer

#2012####
foc_patches<-levels(drop.levels(coinf2012$patche_ID[coinf2012$PA_S2012==1 & coinf2012$PA_2012==1]))
colo<-coinf2012[coinf2012$PA_S2012!=1 | is.na(coinf2012$PA_S2012),]
colo_patches<-levels(drop.levels(colo$patche_ID))
temp<-cbind(closer2012,"foc_patche"=ifelse(as.character(closer2012$patche1_ID) %in% foc_patches,as.character(closer2012$patche1_ID),
                                           as.character(closer2012$patche2_ID)))
temp<-cbind(temp,"prox_patche"=ifelse(as.character(closer2012$patche1_ID) %in% foc_patches,as.character(closer2012$patche2_ID),
                                      as.character(closer2012$patche1_ID)))
temp<-merge(temp,coinf2012,by.x="foc_patche",by.y="patche_ID")
temp<-merge(temp,patche_info,by.x="prox_patche",by.y="ID")
temp<-merge(temp,coinf2012[,1:19],by.x="prox_patche",by.y="patche_ID",all.x=TRUE)
temp<-data.frame(temp,"coinfYN"=temp$number_coinf.x)
temp$coinfYN[(temp$coinfYN)>0]<-1
temp$coinfYN<-as.factor(temp$coinfYN)
temp$road_PA.y<-as.factor(temp$road_PA.y)
temp<-data.frame(temp,"colonized"=ifelse(as.character(temp$prox_patche) %in% colo_patches,1,0))
temp<-data.frame(temp,"Gdiv"=temp$number_MLG.x/temp$number_genotyped.x)
temp<-temp[!is.na(temp$PLM2_Sept2012.y),]
temp<-temp[!is.na(temp$road_PA.y),]
#if you want to limit the radius of the interaction you are interested in
plot(temp$dist,col=((as.numeric(temp$dist)>20)+1))
temp<-temp[temp$dist<20,] #maybe a good idea to limit below 5 kms
prevalcolo<-glm(colonized~cumulative_sum.x+Gdiv+connec2012.x+coinfYN+AA_S2012.x+log(dist),
                family=binomial,data=temp)
summary(prevalcolo)
prevalcolo<-glm(colonized~PLM2_Sept2012.y+coinfYN+AA_S2012.x+log(dist),
                family=binomial,data=temp)
summary(prevalcolo)
prevalcolo<-glm(colonized~I(sqrt(PLM2_Sept2012.y))+coinfYN+AA_S2012.x+log(dist),
                family=binomial,data=temp)
summary(prevalcolo)
op<-par(mfrow=c(1,3))
visreg(prevalcolo,"coinfYN",scale="response",overlay=TRUE,
       xlab="No coinf (0)/ coinf (1)",ylab="P(colonize)")
visreg(prevalcolo,"dist",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Distance to focal",ylab="P(colonize)")
visreg(prevalcolo,"PLM2_Sept2012.y",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Net trapping area",ylab="P(colonize)")
par(op)
mixprevalcolo<-glmer(colonized~cumulative_sum.x+Gdiv+connec2012.x+coinfYN+AA_S2012.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp)
summary(mixprevalcolo)
mixprevalcolo<-glmer(colonized~I(sqrt(PLM2_Sept2012.y))+coinfYN+AA_S2012.x+log(dist)+(1|foc_patche),
                     family=binomial,data=temp)
summary(mixprevalcolo)
#remove focal patches where nothing happens
active<-levels(drop.levels(temp[temp$colonized==1,]$foc_patch))
temp2<-temp[temp$foc_patch %in% active,]
prevalcolo<-glm(colonized~cumulative_sum.x+Gdiv+connec2012.x+coinfYN+AA_S2012.x+log(dist),
                family=binomial,data=temp2)
summary(prevalcolo)
prevalcolo<-glm(colonized~PLM2_Sept2012.y+coinfYN+AA_S2012.x+log(dist),
                family=binomial,data=temp2)
summary(prevalcolo)
prevalcolo<-glm(colonized~I(sqrt(PLM2_Sept2012.y))+coinfYN+AA_S2012.x+log(dist),
                family=binomial,data=temp2)
summary(prevalcolo)
op<-par(mfrow=c(1,3))
visreg(prevalcolo,"coinfYN",scale="response",overlay=TRUE,
       xlab="No coinf (0)/ coinf (1)",ylab="P(colonize)")
visreg(prevalcolo,"dist",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Distance to focal",ylab="P(colonize)")
visreg(prevalcolo,"PLM2_Sept2012.y",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Net trapping area",ylab="P(colonize)")
par(op)

prevalcolo_1<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2012.y+coinfYN+AA_S2012.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_1)
prevalcolo_2<-glm(colonized~Sh.y+road_PA.y+PLM2_Sept2012.y+coinfYN+AA_S2012.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_2)
prevalcolo_3<-glm(colonized~shadow.y+road_PA.y+PLM2_Sept2012.y+coinfYN+AA_S2012.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_3)
prevalcolo_4<-glm(colonized~shadow.y+Sh.y+PLM2_Sept2012.y+coinfYN+AA_S2012.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_4)
prevalcolo_5<-glm(colonized~shadow.y+Sh.y+road_PA.y+coinfYN+AA_S2012.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_5)
prevalcolo_6<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2012.y+AA_S2012.x+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_6)
prevalcolo_7<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2012.y+coinfYN+log(dist),
                  family=binomial,data=temp2)
summary(prevalcolo_7)
prevalcolo_8<-glm(colonized~shadow.y+Sh.y+road_PA.y+PLM2_Sept2012.y+coinfYN+AA_S2012.x,
                  family=binomial,data=temp2)
summary(prevalcolo_8)
anova(prevalcolo_2,prevalcolo_1,test="Chisq")









############
rm(coord,vecdistan, )





###############################################################################
#END
###############################################################################