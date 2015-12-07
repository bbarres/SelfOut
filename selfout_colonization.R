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

#Compute the distances between every patches saves some computation time, 
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


#for 2013#### takes some time to run: ~30 min
closer2013<-closerList(coinf2013,coinf2013$PA_S2013,coinf2013$PA_2013,
                       patche_info$PA_2013)

#for 2012#### takes some time to run: ~30 min
closer2012<-closerList(coinf2012,coinf2012$PA_S2012,coinf2012$PA_2012,
                       patche_info$PA_2012)



eval(parse(text=paste(temp,temp2,sep="")))


###############################################################################
#Modelisation of the colonization success
###############################################################################

#preparing the table before glm analysis for 2013####

foc_patches<-coinf2013[coinf2013$PA_S2013 == 1 & coinf2013$PA_2013 == 1,]
foc_patches<-levels(drop.levels(foc_patches$patche_ID))
colo_patches<-coinf2013[coinf2013$PA_S2013 != 1 | is.na(coinf2013$PA_S2013),]
colo_patches<-levels(drop.levels(colo_patches$patche_ID))
temp <- cbind(
  closer2013,
  "foc_patche" = ifelse(
    as.character(closer2013$patche1_ID) %in% foc_patches,
    as.character(closer2013$patche1_ID),
    as.character(closer2013$patche2_ID)
  )
)
temp <-
  cbind(temp,"prox_patche" = ifelse(
    as.character(closer2013$patche1_ID) %in% foc_patches,
    as.character(closer2013$patche2_ID),
    as.character(closer2013$patche1_ID)
  ))
temp <- merge(temp,coinf2013,by.x = "foc_patche",by.y = "patche_ID")
temp <- merge(temp,patche_info,by.x = "prox_patche",by.y = "ID")
temp <- merge(temp,coinf2013[,1:19],by.x = "prox_patche",by.y = "patche_ID",
              all.x = TRUE)
temp <- data.frame(temp,"coinfYN" = temp$number_coinf.x)
temp$coinfYN[(temp$coinfYN) > 0] <- 1
temp$coinfYN <- as.factor(temp$coinfYN)
temp$road_PA.y <- as.factor(temp$road_PA.y)
temp <- data.frame(temp,
                   "colonized" = ifelse(as.character(temp$prox_patche) %in% 
                                          colo_patches,1,0))
temp <- data.frame(temp,"Gdiv" = temp$number_MLG.x / temp$number_genotyped.x)
#here there is a particular problem with patches in Kökar. This place probably 
#wasn't surveyed in July or no infected patch was detected, therefore the 
#closest focal patch is more than 20 km away. 
plot(temp$dist,col=((as.numeric(temp$dist)>20)+1))
#here the probability that the "closest" focal patch is the patch of origin 
#is very small, so we remove these observation from the dataset
temp<-temp[temp$dist<20,]
#remove focal patches where nothing happens (no colonization event at all). 
#the reason why we remove these patches is because the complete lack 
#of colonization event might have been the consequence of an absence of reel 
#infection in June in the first place (false positive during the spring survey)
active<-levels(drop.levels(temp[temp$colonized==1,]$foc_patch))
temp<-temp[temp$foc_patch %in% active,]

#here is the generalyzed linear model for 2013
modcolo_13 <- glm(
  colonized ~ shadow.y + Sh.y + road_PA.y + PLM2_Sept2013.y + coinfYN + 
              AA_S2013.x + log(dist),
  family = binomial,data = temp
)
summary(modcolo_13)
drop1(modcolo_13,test="Chisq")
#cleaning the environment
rm(temp,active,foc_patches,colo_patches)


#preparing the table before glm analysis for 2012####

foc_patches<-coinf2012[coinf2012$PA_S2012 == 1 & coinf2012$PA_2012 == 1,]
foc_patches<-levels(drop.levels(foc_patches$patche_ID))
colo_patches<-coinf2012[coinf2012$PA_S2012 != 1 | is.na(coinf2012$PA_S2012),]
colo_patches<-levels(drop.levels(colo_patches$patche_ID))
temp <-
  cbind(
    closer2012,"foc_patche" = ifelse(
      as.character(closer2012$patche1_ID) %in% foc_patches,
      as.character(closer2012$patche1_ID),
      as.character(closer2012$patche2_ID)
    )
  )
temp <-
  cbind(temp,"prox_patche" = ifelse(
    as.character(closer2012$patche1_ID) %in% foc_patches,
    as.character(closer2012$patche2_ID),
    as.character(closer2012$patche1_ID)
  ))
#
temp <- merge(temp,coinf2012,by.x = "foc_patche",by.y = "patche_ID")
temp <- merge(temp,patche_info,by.x = "prox_patche",by.y = "ID")
temp <- merge(temp,coinf2012[,1:19],by.x = "prox_patche",by.y = "patche_ID",
              all.x = TRUE)
temp <- data.frame(temp,"coinfYN" = temp$number_coinf.x)
temp$coinfYN[(temp$coinfYN) > 0] <- 1
temp$coinfYN <- as.factor(temp$coinfYN)
temp$road_PA.y <- as.factor(temp$road_PA.y)
temp <- data.frame(temp,
                   "colonized" = ifelse(as.character(temp$prox_patche) %in% 
                                          colo_patches,1,0))
temp <- data.frame(temp,"Gdiv" = temp$number_MLG.x / temp$number_genotyped.x)
#checking that we don't have the same problem as in 2013 with isolated patches
plot(temp$dist,col = ((as.numeric(temp$dist) > 20) + 1))
#remove focal patches where nothing happens (no colonization event at all). 
#the reason why we remove these patches is because the complete lack 
#of colonization event might have been the consequence of an absence of reel 
#infection in June in the first place (false positive during the spring survey)
active <- levels(drop.levels(temp[temp$colonized == 1,]$foc_patch))
temp <- temp[temp$foc_patch %in% active,]

#here is the generalyzed linear model for 2012
modcolo_12 <- glm(
    colonized ~ shadow.y + Sh.y + road_PA.y + PLM2_Sept2012.y + coinfYN + 
      AA_S2012.x + log(dist),
    family = binomial,data = temp
  )
summary(modcolo_12)
drop1(modcolo_12,test="Chisq")
#cleaning the environment
rm(temp,active,foc_patches,colo_patches)

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
  text(x[c(1,3,5)],y[4],labels=labels,adj=c(0.5,0.3),cex=division.cex)
}


###############################################################################
#Figure to examplify the main result of the analysis
###############################################################################

#2013 plotting example of patch colonization####
#we build a table with the focal and proximal patches included in the model
foc_patches<-coinf2013[coinf2013$PA_S2013 == 1 & coinf2013$PA_2013 == 1,]
foc_patches<-levels(drop.levels(foc_patches$patche_ID))
colo_patches<-coinf2013[coinf2013$PA_S2013 != 1 | is.na(coinf2013$PA_S2013),]
colo_patches<-levels(drop.levels(colo_patches$patche_ID))
temp <- cbind(
  closer2013,
  "foc_patche" = ifelse(
    as.character(closer2013$patche1_ID) %in% foc_patches,
    as.character(closer2013$patche1_ID),
    as.character(closer2013$patche2_ID)
  )
)
temp <-
  cbind(temp,"prox_patche" = ifelse(
    as.character(closer2013$patche1_ID) %in% foc_patches,
    as.character(closer2013$patche2_ID),
    as.character(closer2013$patche1_ID)
  ))
temp <- merge(temp,coinf2013,by.x = "foc_patche",by.y = "patche_ID")
temp <- merge(temp,patche_info,by.x = "prox_patche",by.y = "ID")
temp <- merge(temp,coinf2013[,1:19],by.x = "prox_patche",by.y = "patche_ID",
              all.x = TRUE)
temp <- data.frame(temp,"coinfYN" = temp$number_coinf.x)
temp$coinfYN[(temp$coinfYN) > 0] <- 1
temp$coinfYN <- as.factor(temp$coinfYN)
temp$road_PA.y <- as.factor(temp$road_PA.y)
temp <- data.frame(temp,
                   "colonized" = ifelse(as.character(temp$prox_patche) %in% 
                                          colo_patches,1,0))
temp <- data.frame(temp,"Gdiv" = temp$number_MLG.x / temp$number_genotyped.x)
#here there is a particular problem with patches in Kökar. This place probably 
#wasn't surveyed in July or no infected patch was detected, therefore the 
#closest focal patch is more than 20 km away. 
plot(temp$dist,col=((as.numeric(temp$dist)>20)+1))
#here the probability that the "closest" focal patch is the patch of origin 
#is very small, so we remove these observation from the dataset
temp<-temp[temp$dist<20,]

#update the list of focal and proximal patches
foc_patche<-levels(drop.levels(temp$foc_patche))
prox_patche<-levels(drop.levels(temp$prox_patche))
#the combined list of patches
all_patche<-c(foc_patche,prox_patche)
#then we have the patches that are infected at the end of the season
inf_patche<-levels(drop.levels(temp$prox_patche[temp$PA_2013.y==1]))
#then we can emphasize the focal patch with coinfection
#we list the focal patches with coinfection
fopa_coin<-levels(drop.levels(temp$foc_patche[temp$PA_S2013.x==1 & 
                                                temp$PA_2013.x==1 & 
                                                temp$number_coinf.x!=0]))
#and the colonized patches as from a coinfected patch
prpa_coin<-levels(drop.levels(temp$prox_patche[temp$PA_2013.y==1 & 
                                                 temp$foc_patche %in% 
                                                 fopa_coin]))

op<-par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(0,1,0,1))
#we plot the general map
plot(Aland,col="white",lty=1)
points(patche_info[patche_info$ID %in% all_patche,]$Longitude,
       patche_info[patche_info$ID %in% all_patche,]$Latitude,
       cex=0.5,bg="transparent",pch=21,col="black")
points(patche_info[patche_info$ID %in% inf_patche,]$Longitude,
       patche_info[patche_info$ID %in% inf_patche,]$Latitude,
       cex=0.5,bg="lightblue",pch=21,col="black")
points(patche_info[patche_info$ID %in% foc_patche,]$Longitude,
       patche_info[patche_info$ID %in% foc_patche,]$Latitude,
       cex=0.5,bg="darkblue",pch=21,col="black")
points(patche_info[patche_info$ID %in% prpa_coin,]$Longitude,
       patche_info[patche_info$ID %in% prpa_coin,]$Latitude,
       cex=0.5,bg="red",pch=21,col="black")
points(patche_info[patche_info$ID %in% fopa_coin,]$Longitude,
       patche_info[patche_info$ID %in% fopa_coin,]$Latitude,
       cex=0.5,bg="darkred",pch=21,col="black")
polygon(x=c(108700,110800,110800,108700),
        y=c(6673200,6673200,6675100,6675100),
        border="darkorange",lwd=1.5)
mtext(side=2,text="A)",font=2,cex=1,adj=0,padj=-8,las=1,outer=TRUE,line=-1)

#we plot the background of the map
plot(Aland,col="white",lty=0,xlim=c(108700,110800),ylim=c(6673200,6675100))
#we add the entire set of patches
plot(patchshape[patchshape[[3]] %in% all_patche,1]
     ,col="white",lty=1,lwd=1,add=TRUE)
plot(patchshape[patchshape[[3]] %in% all_patche,1]
     ,col="white",lty=0,add=TRUE)
#we add the patches that are infected by the end of the survey
plot(patchshape[patchshape[[3]] %in% inf_patche,1],col="lightblue",
     lty=0,add=TRUE)
#we superimpose the focal patches to the infected patches in the fall
plot(patchshape[patchshape[[3]] %in% foc_patche,1],col="darkblue",
     lty=0,add=TRUE)
#Here we superimposed the coinfection information on the previous map
plot(patchshape[patchshape[[3]] %in% fopa_coin,1],
     col="darkred",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% prpa_coin,1],
     col="red",lty=0,add=TRUE)

#finaly we can represent the relationship between the focal patches and 
#their related candidate for colonization neighbors by a dash line
for (i in 1:length(foc_patches)) {
  temp2<-temp[as.character(temp$foc_patche)==foc_patches[i],]
  for (j in 1:dim(temp2)[1]){
    lines(c(temp2[j,"Longitude.x"],
            temp2[j,"Longitude.y"]),
          c(temp2[j,"Latitude.x"],
            temp2[j,"Latitude.y"]),lwd=1,lty=2,
          col=temp2[j,"PA_2013.y"]+1)
  }
}
polygon(x=c(108700,110800,110800,108700),
        y=c(6673200,6673200,6675100,6675100),
        border="darkorange",lwd=2)
scalebar(c(109700,6673240),1000,"km",division.cex=0.5,xpd=NA)
mtext(side=2,text="B)",font=2,cex=1,adj=0,padj=-8,las=1,outer=TRUE,line=-20)

par(op)

rm(prpa_coin,fopa_coin,temp,foc_patches,colo_patches,
   inf_patche,foc_patche,prox_patche,all_patche)

#export the map to a pdf file 7.93 x 3.64 inches


#2012 plotting example of patch colonization####
#we build a table with the focal and proximal patches included in the model
foc_patches<-coinf2012[coinf2012$PA_S2012 == 1 & coinf2012$PA_2012 == 1,]
foc_patches<-levels(drop.levels(foc_patches$patche_ID))
colo_patches<-coinf2012[coinf2012$PA_S2012 != 1 | is.na(coinf2012$PA_S2012),]
colo_patches<-levels(drop.levels(colo_patches$patche_ID))
temp <- cbind(
  closer2012,
  "foc_patche" = ifelse(
    as.character(closer2012$patche1_ID) %in% foc_patches,
    as.character(closer2012$patche1_ID),
    as.character(closer2012$patche2_ID)
  )
)
temp <-
  cbind(temp,"prox_patche" = ifelse(
    as.character(closer2012$patche1_ID) %in% foc_patches,
    as.character(closer2012$patche2_ID),
    as.character(closer2012$patche1_ID)
  ))
temp <- merge(temp,coinf2012,by.x = "foc_patche",by.y = "patche_ID")
temp <- merge(temp,patche_info,by.x = "prox_patche",by.y = "ID")
temp <- merge(temp,coinf2012[,1:19],by.x = "prox_patche",by.y = "patche_ID",
              all.x = TRUE)
temp <- data.frame(temp,"coinfYN" = temp$number_coinf.x)
temp$coinfYN[(temp$coinfYN) > 0] <- 1
temp$coinfYN <- as.factor(temp$coinfYN)
temp$road_PA.y <- as.factor(temp$road_PA.y)
temp <- data.frame(temp,
                   "colonized" = ifelse(as.character(temp$prox_patche) %in% 
                                          colo_patches,1,0))
temp <- data.frame(temp,"Gdiv" = temp$number_MLG.x / temp$number_genotyped.x)
#here there is a particular problem with patches in Kökar. This place probably 
#wasn't surveyed in July or no infected patch was detected, therefore the 
#closest focal patch is more than 20 km away. 
plot(temp$dist,col=((as.numeric(temp$dist)>20)+1))
#here the probability that the "closest" focal patch is the patch of origin 
#is very small, so we remove these observation from the dataset
temp<-temp[temp$dist<20,]

#update the list of focal and proximal patches
foc_patche<-levels(drop.levels(temp$foc_patche))
prox_patche<-levels(drop.levels(temp$prox_patche))
#the combined list of patches
all_patche<-c(foc_patche,prox_patche)
#then we have the patches that are infected at the end of the season
inf_patche<-levels(drop.levels(temp$prox_patche[temp$PA_2012.y==1]))

#we plot the background of the map
plot(Aland,col="white",lty=0)
#we add the entire set of patches
plot(patchshape[patchshape[[3]] %in% all_patche,1]
     ,col="white",lty=1,lwd=0.1,add=TRUE)
plot(patchshape[patchshape[[3]] %in% all_patche,1]
     ,col="white",lty=0,add=TRUE)
#we add the patches that are infected by the end of the survey
plot(patchshape[patchshape[[3]] %in% inf_patche,1],col="lightblue",
     lty=0,add=TRUE)
#we superimpose the focal patches to the infected patches in the fall
plot(patchshape[patchshape[[3]] %in% foc_patche,1],col="darkblue",
     lty=0,add=TRUE)

#then we can emphasize the focal patch with coinfection
#we list the focal patches with coinfection
fopa_coin<-levels(drop.levels(temp$foc_patche[temp$PA_S2012.x==1 & 
                                                temp$PA_2012.x==1 & 
                                                temp$number_coinf.x!=0]))
#and the colonized patches as from a coinfected patch
prpa_coin<-levels(drop.levels(temp$prox_patche[temp$PA_2012.y==1 & 
                                                 temp$foc_patche %in% 
                                                 fopa_coin]))
#Here we superimposed the coinfection information on the previous map
plot(patchshape[patchshape[[3]] %in% fopa_coin,1],
     col="darkred",lty=0,add=TRUE)
plot(patchshape[patchshape[[3]] %in% prpa_coin,1],
     col="red",lty=0,add=TRUE)

#finaly we can represent the relationship between the focal patche and 
#their related candidate for colonization neighbors by a dash line
for (i in 1:length(foc_patches)) {
  temp2<-temp[as.character(temp$foc_patche)==foc_patches[i],]
  for (j in 1:dim(temp2)[1]){
    lines(c(temp2[j,"Longitude.x"],
            temp2[j,"Longitude.y"]),
          c(temp2[j,"Latitude.x"],
            temp2[j,"Latitude.y"]),lwd=0.05,lty=2,
          col=temp2[j,"PA_2012.y"]+1)
  }
}

rm(prpa_coin,fopa_coin,temp,foc_patches,colo_patches,
   inf_patche,foc_patche,prox_patche,all_patche)

#export the map to a pdf file 20 x 15 inches


############
rm(coord,vecdistan)


###############################################################################
#END
###############################################################################