###############################################################################
###############################################################################
#Importing and preparing the data for analysis
###############################################################################
###############################################################################

setwd("~/work/Rfichiers/Githuber/SelfOut_data")

library(raster)
library(maptools)
library(vegan)
library(mapplots)
library(RColorBrewer)
library(rgdal)
library(gdata)


###############################################################################
#loading the different datasets
###############################################################################

#the code to deal with shapefile instead of RData file for GIS purpose
Aland<-readShapePoly("RECNO.shp",proj4string=CRS("+init=epsg:2393"))
Aland<-spTransform(Aland,CRS("+init=epsg:3067"))
plot(Aland,col=grey(0.7),axes=FALSE,lty=0)
patchshape<-readShapePoly("All patches.shp",proj4string=CRS("+init=epsg:3067"))

plot(Aland,col=grey(0.85),lty=0)
plot(patchshape,col="red",lty=0,add=TRUE)

#loading patches informations
patche_info<-read.table("coord_all_patch12.txt",header=TRUE, sep="\t",dec=".")
#we reorganize the coord data file
patche_info<-patche_info[,c(1:12,14:dim(patche_info)[2])]
points(patche_info$Longitude,patche_info$Latitude,cex=1,bg="white",pch=21)

#loading the sample information data
sample_info<-read.table("sample_info4.txt", header=TRUE, sep="\t", dec=".")

#loading genotypes data
geno_hom<-read.table("geno_hom10_13.txt",header=TRUE,sep="\t",
                     stringsAsFactors=FALSE)
geno_hom<-merge(geno_hom,sample_info,by.x="UNIC_ID",by.y="FIMM_ID")
geno_hom<-merge(geno_hom,patche_info,by.x="patche_ID",by.y="ID",all.x=TRUE)
geno_hom<-replace(geno_hom,geno_hom=="CA","AC")
geno_hom<-replace(geno_hom,geno_hom=="GA","AG")
geno_hom<-replace(geno_hom,geno_hom=="TA","AT")
geno_hom<-replace(geno_hom,geno_hom=="GC","CG")
geno_hom<-replace(geno_hom,geno_hom=="TC","CT")
geno_hom<-replace(geno_hom,geno_hom=="TG","GT")
geno_hom<-droplevels(geno_hom)


###############################################################################
#adding the MLG column
###############################################################################

#we add a column with the multilocus genotype (MLG)
MLG<-geno_hom
MLG<-replace(MLG,MLG=="AA",1)
MLG<-replace(MLG,MLG=="CC",2)
MLG<-replace(MLG,MLG=="GG",3)
MLG<-replace(MLG,MLG=="TT",4)
MLGmat<-vector()
nb_SNP<- 19 #the number of markers used
name_SNP<-dimnames(geno_hom)[[2]][3:(nb_SNP+3-1)]
for (i in 3:(nb_SNP+3-1)) MLGmat<-paste(MLGmat,MLG[,i],sep="/")
geno_hom<-data.frame(geno_hom,"MLG"=MLGmat,stringsAsFactors=FALSE)


###############################################################################
#computing connectivity indices of the patches
###############################################################################

#let's compute the connectivity index in 2010
coor2010<-patche_info[!is.na(patche_info$PA_2010),]
coor2010<-coor2010[coor2010$PA_2010==1,]
coor2010<-coor2010[order(coor2010$ID),]
#xy<-as.matrix(data2010pop@other$xy)
xy<-as.matrix(coor2010[,c("Longitude","Latitude")])
dimnames(xy)[[1]]<-coor2010$ID
Weidist<-vegdist(xy,method="euclidean",diag=TRUE, upper=TRUE)
#'alpha' parameter of Laine and Hanski 2006, here we put 1km (take care that 
#the same unit is used in the Weidist matrix)
alpha<-1
Weidist<-as.matrix(Weidist/1000) #conversion from m to km
#conversion of area in square meter to the square root of the area in square km
racArea<-sqrt(coor2010$fallPLM2_2001_2008/1000000)
#if some area informations are lacking, we assume that their area is very small
racArea[is.na(racArea)]<-0 
Weidist<-exp(-alpha*Weidist)
connec2010<-matrix(nrow=(dim(Weidist)[1]),ncol=(dim(Weidist)[2]))
dimnames(connec2010)<-dimnames(Weidist)
for (i in 1:length(racArea)) {
  connec2010[,i]<-Weidist[,i]*racArea
}
for (i in 1:length(racArea)) {
  connec2010[i,i]<-0
}
connec2010<-apply(connec2010,2,sum)
connec2010<-cbind("patche_ID"=attr(connec2010,"names"),"connec2010"=connec2010)

#let's compute the connectivity index in 2011
coor2011<-patche_info[!is.na(patche_info$PA_2011),]
coor2011<-coor2011[coor2011$PA_2011==1,]
coor2011<-coor2011[order(coor2011$ID),]
#xy<-as.matrix(data2011pop@other$xy)
xy<-as.matrix(coor2011[,c("Longitude","Latitude")])
dimnames(xy)[[1]]<-coor2011$ID
Weidist<-vegdist(xy,method="euclidean",diag=TRUE, upper=TRUE)
#'alpha' parameter of Laine and Hanski 2006, here we put 1km (take care that 
#the same unit is used in the Weidist matrix)
alpha<-1
Weidist<-as.matrix(Weidist/1000) #conversion from m to km
#conversion of area in square meter to the square root of the area in square km
racArea<-sqrt(coor2011$fallPLM2_2001_2008/1000000)
#if some area informations are lacking, we assume that their area is very small
racArea[is.na(racArea)]<-0 
Weidist<-exp(-alpha*Weidist)
connec2011<-matrix(nrow=(dim(Weidist)[1]),ncol=(dim(Weidist)[2]))
dimnames(connec2011)<-dimnames(Weidist)
for (i in 1:length(racArea)) {
  connec2011[,i]<-Weidist[,i]*racArea
}
for (i in 1:length(racArea)) {
  connec2011[i,i]<-0
}
connec2011<-apply(connec2011,2,sum)
connec2011<-cbind("patche_ID"=attr(connec2011,"names"),"connec2011"=connec2011)

#let's compute the connectivity index in fall 2012
coor2012<-patche_info[!is.na(patche_info$PA_2012),]
coor2012<-coor2012[coor2012$PA_2012==1,]
coor2012<-coor2012[order(coor2012$ID),]
#xy<-as.matrix(data2012pop@other$xy)
xy<-as.matrix(coor2012[,c("Longitude","Latitude")])
dimnames(xy)[[1]]<-coor2012$ID
Weidist<-vegdist(xy,method="euclidean",diag=TRUE, upper=TRUE)
#'alpha' parameter of Laine and Hanski 2006, here we put 1km (take care that 
#the same unit is used in the Weidist matrix)
alpha<-1
Weidist<-as.matrix(Weidist/1000) #conversion from m to km
#conversion of area in square meter to the square root of the area in square km
racArea<-sqrt(coor2012$fallPLM2_2001_2008/1000000)
#if some area informations are lacking, we assume that their area is very small
racArea[is.na(racArea)]<-0 
Weidist<-exp(-alpha*Weidist)
connec2012<-matrix(nrow=(dim(Weidist)[1]),ncol=(dim(Weidist)[2]))
dimnames(connec2012)<-dimnames(Weidist)
for (i in 1:length(racArea)) {
  connec2012[,i]<-Weidist[,i]*racArea
}
for (i in 1:length(racArea)) {
  connec2012[i,i]<-0
}
connec2012<-apply(connec2012,2,sum)
connec2012<-cbind("patche_ID"=attr(connec2012,"names"),"connec2012"=connec2012)

#let's compute the connectivity index in spring 2012
coor2012s<-patche_info[!is.na(patche_info$PA_S2012),]
coor2012s<-coor2012s[coor2012s$PA_S2012==1,]
coor2012s<-coor2012s[order(coor2012s$ID),]
#xy<-as.matrix(data2012pop@other$xy)
xy<-as.matrix(coor2012s[,c("Longitude","Latitude")])
dimnames(xy)[[1]]<-coor2012s$ID
Weidist<-vegdist(xy,method="euclidean",diag=TRUE, upper=TRUE)
#'alpha' parameter of Laine and Hanski 2006, here we put 1km (take care that 
#the same unit is used in the Weidist matrix)
alpha<-1
Weidist<-as.matrix(Weidist/1000) #conversion from m to km
#conversion of area in square meter to the square root of the area in square km
racArea<-sqrt(coor2012s$fallPLM2_2001_2008/1000000)
#if some area informations are lacking, we assume that their area is very small
racArea[is.na(racArea)]<-0 
Weidist<-exp(-alpha*Weidist)
connec2012s<-matrix(nrow=(dim(Weidist)[1]),ncol=(dim(Weidist)[2]))
dimnames(connec2012s)<-dimnames(Weidist)
for (i in 1:length(racArea)) {
  connec2012s[,i]<-Weidist[,i]*racArea
}
for (i in 1:length(racArea)) {
  connec2012s[i,i]<-0
}
connec2012s<-apply(connec2012s,2,sum)
connec2012s<-cbind("patche_ID"=attr(connec2012s,"names"),
                   "connec2012s"=connec2012s)

#let's compute the connectivity index in fall 2013
coor2013<-patche_info[!is.na(patche_info$PA_2013),]
coor2013<-coor2013[coor2013$PA_2013==1,]
coor2013<-coor2013[order(coor2013$ID),]
#xy<-as.matrix(data2013pop@other$xy)
xy<-as.matrix(coor2013[,c("Longitude","Latitude")])
dimnames(xy)[[1]]<-coor2013$ID
Weidist<-vegdist(xy,method="euclidean",diag=TRUE, upper=TRUE)
#'alpha' parameter of Laine and Hanski 2006, here we put 1km (take care that 
#the same unit is used in the Weidist matrix)
alpha<-1
Weidist<-as.matrix(Weidist/1000) #conversion from m to km
#conversion of area in square meter to the square root of the area in square km
racArea<-sqrt(coor2013$fallPLM2_2001_2008/1000000)
#if some area informations are lacking, we assume that their area is very small
racArea[is.na(racArea)]<-0 
Weidist<-exp(-alpha*Weidist)
connec2013<-matrix(nrow=(dim(Weidist)[1]),ncol=(dim(Weidist)[2]))
dimnames(connec2013)<-dimnames(Weidist)
for (i in 1:length(racArea)) {
  connec2013[,i]<-Weidist[,i]*racArea
}
for (i in 1:length(racArea)) {
  connec2013[i,i]<-0
}
connec2013<-apply(connec2013,2,sum)
connec2013<-cbind("patche_ID"=attr(connec2013,"names"),"connec2013"=connec2013)

#let's compute the connectivity index in spring 2013
coor2013s<-patche_info[!is.na(patche_info$PA_S2013),]
coor2013s<-coor2013s[coor2013s$PA_S2013==1,]
coor2013s<-coor2013s[order(coor2013s$ID),]
#xy<-as.matrix(data2013pop@other$xy)
xy<-as.matrix(coor2013s[,c("Longitude","Latitude")])
dimnames(xy)[[1]]<-coor2013s$ID
Weidist<-vegdist(xy,method="euclidean",diag=TRUE, upper=TRUE)
#'alpha' parameter of Laine and Hanski 2006, here we put 1km (take care that 
#the same unit is used in the Weidist matrix)
alpha<-1
Weidist<-as.matrix(Weidist/1000) #conversion from m to km
#conversion of area in square meter to the square root of the area in square km
racArea<-sqrt(coor2013s$fallPLM2_2001_2008/1000000)
#if some area informations are lacking, we assume that their area is very small
racArea[is.na(racArea)]<-0
Weidist<-exp(-alpha*Weidist)
connec2013s<-matrix(nrow=(dim(Weidist)[1]),ncol=(dim(Weidist)[2]))
dimnames(connec2013s)<-dimnames(Weidist)
for (i in 1:length(racArea)) {
  connec2013s[,i]<-Weidist[,i]*racArea
}
for (i in 1:length(racArea)) {
  connec2013s[i,i]<-0
}
connec2013s<-apply(connec2013s,2,sum)
connec2013s<-cbind("patche_ID"=attr(connec2013s,"names"),
                   "connec2013s"=connec2013s)


###############################################################################
#computation of the different indice of presence/absence 
###############################################################################

#first we design a function to count the uninterupted cumulated number of years 
#of presence of powdery mildew 
cumul<-function(x) {
  if (is.na(x[length(x)])) {
    n<-"NA"} else {
      if (x[length(x)]==0) {
        n<-0} else {
          n<-1
          for (i in length(x):2) {
            if (is.na(x[i]-x[(i-1)])) {
              break} else {
                if (x[i]-x[(i-1)]==0) {
                  n<-n+1} else {
                    break
                  }
              }
          }
        }
    }
  return(n)
}

#computation for 2013####

#here you have the matrix of presence/absence for considered years, last column
#have to be the sampling year
PresAbs<-patche_info[,17:29] 
#counting the number of years where presence of powdery mildew was signaled
prez_sum<-cbind("patche_ID"=patche_info$ID,
                "presence_sum"=rowSums(PresAbs,na.rm=TRUE))
#mean number of occurence of mildew (very similar to the previous index)
prez_mean<-cbind("patche_ID"=patche_info$ID,
                 "presence_mean"=rowMeans(PresAbs,na.rm=TRUE))
#cumulative number of presence years before sampling
cumu_sum<-cbind("patche_ID"=patche_info$ID,
                "cumulative_sum"=apply(PresAbs,1,cumul))
#we merge the different indices of presence/absence
presabs_index13<-cbind(prez_sum,"presence_mean"=prez_mean[,2],
                       "cumulative_sum"=cumu_sum[,2])

#computation for 2012####

#here you have the matrix of presence/absence for considered years, last column
#have to be the sampling year
PresAbs<-patche_info[,17:28]
#counting the number of years where presence of powdery mildew was signaled
prez_sum<-cbind("patche_ID"=patche_info$ID,
                "presence_sum"=rowSums(PresAbs,na.rm=TRUE))
#mean number of occurence of mildew (very similar to the previous index)
prez_mean<-cbind("patche_ID"=patche_info$ID,
                 "presence_mean"=rowMeans(PresAbs,na.rm=TRUE))
#cumulative number of presence years before sampling
cumu_sum<-cbind("patche_ID"=patche_info$ID,
                "cumulative_sum"=apply(PresAbs,1,cumul))
#we merge the different indices of presence/absence
presabs_index12<-cbind(prez_sum,"presence_mean"=prez_mean[,2],
                       "cumulative_sum"=cumu_sum[,2])

#computation for 2011####

#here you have the matrix of presence/absence for considered years, last column
#have to be the sampling year
PresAbs<-patche_info[,17:27]
#counting the number of years where presence of powdery mildew was signaled
prez_sum<-cbind("patche_ID"=patche_info$ID,
                "presence_sum"=rowSums(PresAbs,na.rm=TRUE))
#mean number of occurence of mildew (very similar to the previous index)
prez_mean<-cbind("patche_ID"=patche_info$ID,
                 "presence_mean"=rowMeans(PresAbs,na.rm=TRUE))
#cumulative number of presence years before sampling
cumu_sum<-cbind("patche_ID"=patche_info$ID,
                "cumulative_sum"=apply(PresAbs,1,cumul))
#we merge the different indices of presence/absence
presabs_index11<-cbind(prez_sum,"presence_mean"=prez_mean[,2],
                       "cumulative_sum"=cumu_sum[,2])

#computation for 2010####

#here you have the matrix of presence/absence for considered years, last column
#have to be the sampling year
PresAbs<-patche_info[,17:26]
#counting the number of years where presence of powdery mildew was signaled
prez_sum<-cbind("patche_ID"=patche_info$ID,
                "presence_sum"=rowSums(PresAbs,na.rm=TRUE))
#mean number of occurence of mildew (very similar to the previous index)
prez_mean<-cbind("patche_ID"=patche_info$ID,
                 "presence_mean"=rowMeans(PresAbs,na.rm=TRUE))
#cumulative number of presence years before sampling
cumu_sum<-cbind("patche_ID"=patche_info$ID,
                "cumulative_sum"=apply(PresAbs,1,cumul))
#we merge the different indices of presence/absence
presabs_index10<-cbind(prez_sum,"presence_mean"=prez_mean[,2],
                       "cumulative_sum"=cumu_sum[,2])


###############################################################################
#computing some statistics about the MLG composition of the populations
###############################################################################

genotype<-geno_hom
#splitting the dataset into one dataset for each year, at the same time we 
#remove all dubious individuals (suspected of co-infection or with too many 
#missing data, ie more than one here)
geno2010<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS10" & 
                     genotype$nb_snp_het<1 & genotype$nb_missing<1 & 
                     genotype$nb_snp_typed>18 & genotype$duplicate!="DUP" & 
                     genotype$leave_ID=="a",]
geno2011<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS11" & 
                     genotype$nb_snp_het<1 & genotype$nb_missing<1 & 
                     genotype$nb_snp_typed>18 & genotype$duplicate!="DUP" & 
                     genotype$leave_ID=="a",]
geno2012<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS12" & 
                     genotype$nb_snp_het<1 & genotype$nb_missing<1 & 
                     genotype$nb_snp_typed>18 & genotype$duplicate!="DUP" & 
                     genotype$leave_ID=="a",]
geno2013<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS13" & 
                     genotype$nb_snp_het<1 & genotype$nb_missing<1 & 
                     genotype$nb_snp_typed>18 & genotype$duplicate!="DUP" & 
                     genotype$leave_ID=="a",]

#we add a column which identified the new genotypes and the genotypes already
#present in 2012, 2011 or 2010
geno2013<-cbind(geno2013,"new_old2012"=ifelse((as.factor(geno2013$MLG)) 
                      %in% (levels(as.factor(geno2012$MLG))),"OLD_2012","NEW"))
geno2013<-cbind(geno2013,"new_old2011"=ifelse((as.factor(geno2013$MLG)) 
                      %in% (levels(as.factor(geno2011$MLG))),"OLD_2011","NEW"))
geno2013<-cbind(geno2013,"new_old2010"=ifelse((as.factor(geno2013$MLG)) 
                      %in% (levels(as.factor(geno2010$MLG))),"OLD_2010","NEW"))

#we add a column which identified the new genotypes and the genotypes already
#present in 2011 or 2010
geno2012<-cbind(geno2012,"new_old2011"=ifelse((as.factor(geno2012$MLG)) 
                      %in% (levels(as.factor(geno2011$MLG))),"OLD_2011","NEW"))
geno2012<-cbind(geno2012,"new_old2010"=ifelse((as.factor(geno2012$MLG)) 
                      %in% (levels(as.factor(geno2010$MLG))),"OLD_2010","NEW"))

#we add a column which identified the new genotypes and the genotypes already
#present in 2010
geno2011<-cbind(geno2011,"new_old2010"=ifelse((as.factor(geno2011$MLG)) 
                      %in% (levels(as.factor(geno2010$MLG))),"OLD_2010","NEW"))

#adding a column for new genotype at the SIN level
tt<-c()
ttt<-c()
for (i in 1:86) {
  tt<-cbind("sample_ID"=as.character(geno2012$sample_ID[geno2012$SIN_86==i & 
      !is.na(geno2012$SIN_86)]),
      "new_old_SINlev"=((ifelse((as.factor(geno2012$MLG[geno2012$SIN_86==i & 
      !is.na(geno2012$SIN_86)])) 
      %in% (levels(as.factor(geno2011$MLG[geno2011$SIN_86==i & 
      !is.na(geno2011$SIN_86)]))),"OLD_2011","NEW"))))
  ttt<-rbind(ttt,tt)
}
ttt<-as.data.frame(ttt)
geno2012<-merge(geno2012,ttt,by="sample_ID",all.x=TRUE)

tt<-c()
ttt<-c()
for (i in 1:86) {
  tt<-cbind("sample_ID"=as.character(geno2013$sample_ID[geno2013$SIN_86==i & 
      !is.na(geno2013$SIN_86)]),
      "new_old_SINlev"=((ifelse((as.factor(geno2013$MLG[geno2013$SIN_86==i & 
      !is.na(geno2013$SIN_86)])) 
      %in% (levels(as.factor(geno2012$MLG[geno2012$SIN_86==i & 
      !is.na(geno2012$SIN_86)]))),"OLD_2012","NEW"))))
  ttt<-rbind(ttt,tt)
}
ttt<-as.data.frame(ttt)
geno2013<-merge(geno2013,ttt,by="sample_ID",all.x=TRUE)

#relatedness between coinfection strains
coinr2010<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS10" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a" & 
                    genotype$nb_snp_het>0,]
coinr2010<-drop.levels(coinr2010)
coinr2010<-tapply(coinr2010$nb_snp_het,INDEX=coinr2010$patche_ID,FUN=mean)

coinr2011<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS11" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a" & 
                    genotype$nb_snp_het>0,]
coinr2011<-drop.levels(coinr2011)
coinr2011<-tapply(coinr2011$nb_snp_het,INDEX=coinr2011$patche_ID,FUN=mean)

coinr2012<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS12" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a" & 
                    genotype$nb_snp_het>0,]
coinr2012<-drop.levels(coinr2012)
coinr2012<-tapply(coinr2012$nb_snp_het,INDEX=coinr2012$patche_ID,FUN=mean)

coinr2013<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS13" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a" & 
                    genotype$nb_snp_het>0,]
coinr2013<-drop.levels(coinr2013)
coinr2013<-tapply(coinr2013$nb_snp_het,INDEX=coinr2013$patche_ID,FUN=mean)

#to investigate the coinfection rate, we transform the nb_snp_het column
geno2010t<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS10" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2010t[(geno2010t$nb_snp_het>0),"nb_snp_het"]<-"COINF"
geno2010t[(geno2010t$nb_snp_het==0),"nb_snp_het"]<-"PURE"
table(geno2010t$nb_snp_het,geno2010t$patche_ID)
geno2011t<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS11" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2011t[(geno2011t$nb_snp_het>0),"nb_snp_het"]<-"COINF"
geno2011t[(geno2011t$nb_snp_het==0),"nb_snp_het"]<-"PURE"
geno2012t<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS12" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2012t[(geno2012t$nb_snp_het>0),"nb_snp_het"]<-"COINF"
geno2012t[(geno2012t$nb_snp_het==0),"nb_snp_het"]<-"PURE"
geno2013t<-genotype[!is.na(genotype$leave_ID) & genotype$survey=="MS13" & 
                    genotype$nb_missing<1 & genotype$nb_snp_typed>18 & 
                    genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2013t[(geno2013t$nb_snp_het>0),"nb_snp_het"]<-"COINF"
geno2013t[(geno2013t$nb_snp_het==0),"nb_snp_het"]<-"PURE"

op<-par(mfrow=c(2,2))
plot(Aland,col="grey90")
title(main="2010 genotyped patches")
points(geno2010$Longitude,geno2010$Latitude,cex=1,bg="red",pch=21)
plot(Aland,col="grey90")
title(main="2011 genotyped patches")
points(geno2011$Longitude,geno2011$Latitude,cex=1,bg="blue",pch=21)
plot(Aland,col="grey90")
title(main="2012 genotyped patches")
points(geno2012$Longitude,geno2012$Latitude,cex=1,bg="green",pch=21)
plot(Aland,col="grey90")
title(main="2013 genotyped patches")
points(geno2013$Longitude,geno2013$Latitude,cex=1,bg="yellow",pch=21)
par(op)

#we can compute a contengy table for genotype by patche for each year
table(geno2010$MLG,geno2010$patche_ID)
#number of genotyped individuals per patch
colSums(table(geno2010$MLG,geno2010$patche_ID))
#number of different MLG per patch
colSums(ifelse(table(geno2010$MLG,geno2010$patche_ID)==0,0,1))
nb_ind_MLG2010<-cbind("number_coinf"=table(geno2010t$nb_snp_het,
                                           geno2010t$patche_ID)[1,],
                      "number_pure"=table(geno2010t$nb_snp_het,
                                          geno2010t$patche_ID)[2,],
                      "number_genotyped"=table(geno2010t$patche_ID))
nb_ind_MLG2010<-cbind("patche_ID"=row.names(nb_ind_MLG2010),nb_ind_MLG2010)
temp<-cbind("number_MLG" = colSums(ifelse(table(geno2010$MLG,
                                                geno2010$patche_ID)==0,0,1)))
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2010<-merge(nb_ind_MLG2010,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("patche_ID"=row.names(coinr2010),coinr2010)
nb_ind_MLG2010<-merge(nb_ind_MLG2010,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)

table(geno2011$MLG,geno2011$patche_ID)
#number of genotyped individuals per patch
colSums(table(geno2011$MLG,geno2011$patche_ID))
#number of different MLG per patch
colSums(ifelse(table(geno2011$MLG,geno2011$patche_ID)==0,0,1))
nb_ind_MLG2011<-cbind("number_coinf"=table(geno2011t$nb_snp_het,
                                           geno2011t$patche_ID)[1,],
                      "number_pure"=table(geno2011t$nb_snp_het,
                                          geno2011t$patche_ID)[2,],
                      "number_genotyped"=table(geno2011t$patche_ID))
nb_ind_MLG2011<-cbind("patche_ID"=row.names(nb_ind_MLG2011),nb_ind_MLG2011)
temp<-cbind("number_MLG" = colSums(ifelse(table(geno2011$MLG,
                                                geno2011$patche_ID)==0,0,1)))
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2011<-merge(nb_ind_MLG2011,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("number_pure_new"=table(geno2011$new_old2010,
                                    geno2011$patche_ID)[1,])
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2011<-merge(nb_ind_MLG2011,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("patche_ID"=row.names(coinr2011),coinr2011)
nb_ind_MLG2011<-merge(nb_ind_MLG2011,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)

table(geno2012$MLG,geno2012$patche_ID)
#number of genotyped individuals per patch
colSums(table(geno2012$MLG,geno2012$patche_ID))
#number of different MLG per patch
colSums(ifelse(table(geno2012$MLG,geno2012$patche_ID)==0,0,1))
nb_ind_MLG2012<-cbind("number_coinf"=table(geno2012t$nb_snp_het,
                                           geno2012t$patche_ID)[1,],
                      "number_pure"=table(geno2012t$nb_snp_het,
                                          geno2012t$patche_ID)[2,],
                      "number_genotyped"=table(geno2012t$patche_ID))
nb_ind_MLG2012<-cbind("patche_ID"=row.names(nb_ind_MLG2012),nb_ind_MLG2012)
temp<-cbind("number_MLG" = colSums(ifelse(table(geno2012$MLG,
                                                geno2012$patche_ID)==0,0,1)))
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2012<-merge(nb_ind_MLG2012,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("number_pure_new"=table(geno2012$new_old2011,
                                    geno2012$patche_ID)[1,])
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2012<-merge(nb_ind_MLG2012,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("number_pure_new_SIN"=table(geno2012$new_old_SINlev,
                                        geno2012$patche_ID)[1,])
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2012<-merge(nb_ind_MLG2012,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("patche_ID"=row.names(coinr2012),coinr2012)
nb_ind_MLG2012<-merge(nb_ind_MLG2012,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)

table(geno2013$MLG,geno2013$patche_ID)
#number of genotyped individuals per patch
colSums(table(geno2013$MLG,geno2013$patche_ID))
#number of different MLG per patch
colSums(ifelse(table(geno2013$MLG,geno2013$patche_ID)==0,0,1)) 
nb_ind_MLG2013<-cbind("number_coinf"=table(geno2013t$nb_snp_het,
                                           geno2013t$patche_ID)[1,],
                      "number_pure"=table(geno2013t$nb_snp_het,
                                          geno2013t$patche_ID)[2,],
                      "number_genotyped"=table(geno2013t$patche_ID))
nb_ind_MLG2013<-cbind("patche_ID"=row.names(nb_ind_MLG2013),nb_ind_MLG2013)
temp<-cbind("number_MLG" = colSums(ifelse(table(geno2013$MLG,
                                                geno2013$patche_ID)==0,0,1)))
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2013<-merge(nb_ind_MLG2013,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("number_pure_new"=table(geno2013$new_old2012,
                                    geno2013$patche_ID)[1,])
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2013<-merge(nb_ind_MLG2013,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("number_pure_new_SIN"=table(geno2013$new_old_SINlev,
                                        geno2013$patche_ID)[1,])
temp<-cbind("patche_ID"=row.names(temp),temp)
nb_ind_MLG2013<-merge(nb_ind_MLG2013,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)
temp<-cbind("patche_ID"=row.names(coinr2013),coinr2013)
nb_ind_MLG2013<-merge(nb_ind_MLG2013,temp,by.x="patche_ID",by.y="patche_ID",
                      all=TRUE)

stats_patch2010<-nb_ind_MLG2010
stats_patch2011<-nb_ind_MLG2011
stats_patch2012<-nb_ind_MLG2012
stats_patch2013<-nb_ind_MLG2013


###############################################################################
#combinning stats, presabs and connectivity indices
###############################################################################

temp<-merge(stats_patch2010,presabs_index10,by.x="patche_ID",by.y="patche_ID")
temp<-merge(temp,connec2010,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2011,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012s,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2013,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
nb_ind_MLG2010<-merge(temp,connec2013s,by.x="patche_ID",by.y="patche_ID",
                      all.x=TRUE)

temp<-merge(stats_patch2011,presabs_index11,by.x="patche_ID",by.y="patche_ID")
temp<-merge(temp,connec2010,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2011,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012s,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2013,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
nb_ind_MLG2011<-merge(temp,connec2013s,by.x="patche_ID",by.y="patche_ID",
                      all.x=TRUE)

temp<-merge(stats_patch2012,presabs_index12,by.x="patche_ID",by.y="patche_ID")
temp<-merge(temp,connec2010,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2011,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012s,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2013,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
nb_ind_MLG2012<-merge(temp,connec2013s,by.x="patche_ID",by.y="patche_ID",
                      all.x=TRUE)

temp<-merge(stats_patch2013,presabs_index13,by.x="patche_ID",by.y="patche_ID")
temp<-merge(temp,connec2010,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2011,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2012s,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
temp<-merge(temp,connec2013,by.x="patche_ID",by.y="patche_ID",all.x=TRUE)
nb_ind_MLG2013<-merge(temp,connec2013s,by.x="patche_ID",by.y="patche_ID",
                      all.x=TRUE)

nb_ind_MLG2010<-merge(nb_ind_MLG2010,patche_info,by.x="patche_ID",by.y="ID")
nb_ind_MLG2011<-merge(nb_ind_MLG2011,patche_info,by.x="patche_ID",by.y="ID")
nb_ind_MLG2012<-merge(nb_ind_MLG2012,patche_info,by.x="patche_ID",by.y="ID")
nb_ind_MLG2013<-merge(nb_ind_MLG2013,patche_info,by.x="patche_ID",by.y="ID")
stat_patch2010<-nb_ind_MLG2010[nb_ind_MLG2010$number_genotyped!="0",]
stat_patch2011<-nb_ind_MLG2011[nb_ind_MLG2011$number_genotyped!="0",]
stat_patch2012<-nb_ind_MLG2012[nb_ind_MLG2012$number_genotyped!="0",]
stat_patch2013<-nb_ind_MLG2013[nb_ind_MLG2013$number_genotyped!="0",]


###############################################################################
#export the table
###############################################################################

write.table(stat_patch2010,file="stat_patch2010.txt",row.names=FALSE,
            quote=FALSE,sep="\t")
write.table(stat_patch2011,file="stat_patch2011.txt",row.names=FALSE,
            quote=FALSE,sep="\t")
write.table(stat_patch2012,file="stat_patch2012.txt",row.names=FALSE,
            quote=FALSE,sep="\t")
write.table(stat_patch2013,file="stat_patch2013.txt",row.names=FALSE,
            quote=FALSE,sep="\t")


