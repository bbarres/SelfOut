###############################################################################
###############################################################################
#Importing and preparing the data for analysis
###############################################################################
###############################################################################

library(raster)
library(maptools)
library(vegan)
library(mapplots)
library(RColorBrewer)
library(rgdal)
library(gdata)
library(adegenet)
library(combinat)

#the following code aims at importing the raw data (Aland shapefile for the 
#map, information on host patches and genotype information for individual
#samples). The information on patches and genotypes are then combined 
#and split into four different files, one per year. Several indices are
#computed and added to the files: several connectivity indices between 
#patches, the prevalence of coinfection, the number of different 
#Multilocus Genotypes, several genetic diversity indices, different indices
#of the prevalence of powdery mildew in the patches. Then all these 
#information is combined in one table for each studied years (2010, 2011, 
#2012 and 2013) and exported so that all these indices can be used without 
#having to process the raw data each time. The code can be a bit tedious and 
#can certainly be improved (both for clarity and efficiency). However, this 
#code is stored and available on github for us to keep track of what we have 
#done and in case someone would be interested in some details of the data 
#processing


###############################################################################
#loading the different datasets
###############################################################################

#the code to deal with shapefile instead of RData file for GIS purpose
Aland<-readShapePoly("data/RECNO.shp",proj4string=CRS("+init=epsg:2393"))
Aland<-spTransform(Aland,CRS("+init=epsg:3067"))
plot(Aland,col=grey(0.7),axes=FALSE,lty=0)
patchshape<-readShapePoly("All patches.shp",proj4string=CRS("+init=epsg:3067"))

plot(Aland,col=grey(0.85),lty=0)
plot(patchshape,col="red",lty=0,add=TRUE)

#loading patches informations
patche_info<-read.table("data/coord_all_patch12.txt",header=TRUE, sep="\t",dec=".")
#we reorganize the coord data file
patche_info<-patche_info[,c(1:12,14:dim(patche_info)[2])]
points(patche_info$Longitude,patche_info$Latitude,cex=1,bg="white",pch=21)

#loading the sample information data
sample_info<-read.table("data/sample_info4.txt", header=TRUE, sep="\t", dec=".")

#loading genotypes data
geno_hom<-read.table("data/geno_hom10_13.txt",header=TRUE,sep="\t",
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

#we add a column with the multilocus genotype (MLG). This is the combination 
#of the 19 individual SNPs information. 
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
#computing allelic richness and genotypic richness
###############################################################################

#definition of an function to compute allelic richness using a genind object 
#for infile

AllRich<-function(data) #data is a genind object (see adegenet for details)
{
  #Conversion from 'genind' object to 'genpop' object
  datapop<-genind2genpop(data, process.other=TRUE, other.action=mean)
  #First, determining the smaller number of allele across sampled population
  matloc<-t(matrix(data=datapop@loc.fac,nrow=(dim(datapop@tab)[2]), 
                   ncol=(dim(datapop@tab)[1])))
  matpop<-matrix(data=rownames(datapop@tab), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]))
#   matpop<-matrix(data=datapop@pop.names, nrow=(dim(datapop@tab)[1]), 
#                  ncol=(dim(datapop@tab)[2]))
  conf<-list(matpop, matloc)
  effN<-(tapply(datapop@tab, conf, sum))
  effN<-effN[order(as.numeric(rownames(effN))),]
  colnames(effN)<-locNames(data)
  echMin<-min(effN)
  
  #Second, build of the matrix of total number of sampled allele 
  truc<-t(as.matrix(table(datapop@loc.fac)))
  x<-matrix(nrow=(dim(effN)[1]), ncol=(dim(effN)[2]), data=truc,byrow=TRUE)
  effTot<-matrix(rep(t(effN),t(x)), nrow=(dim(datapop@tab)[1]), 
                 ncol=(dim(datapop@tab)[2]), byrow=TRUE)
  
  #Third, compute the matrix of Ar for each population/loci combination
  #(see El Mousadik and Petit 1996 for details)
  CoMat<-matrix(nrow=(dim(datapop@tab)[1]),ncol=(dim(datapop@tab)[2]))
  for (i in 1:(dim(datapop@tab)[1])) {
    for (j in 1:(dim(datapop@tab)[2])) {
      CoMat[i,j]<-(1-(nCm(effTot[i,j]-datapop@tab[i,j],echMin)/
                        nCm(effTot[i,j],echMin)))
    }
  }
  
  #Allelic richness in each population, for each LOCUS
  ArLOC<-(tapply(CoMat, conf, sum))
  ArLOC<-ArLOC[order(as.numeric(dimnames(ArLOC)[[1]])),]
  colnames(ArLOC)<-locNames(data)
  ##determining mean Allelic Richness across site and loci
  #determining mean Allelic Richness across loci
  #Ar<-(apply(ArLOC,1,mean))
  rez<-list("Minimum Sampling Size"=echMin,"Allelic Richness Matrix"=ArLOC)
  return(rez)
}

#Let's have a fresh start and load again the original genotype information 
genotype<-geno_hom
#In order to compute allelic or genotypic richness, it is important to set 
#a minimum sample size. The AllRich function set automaticaly the minimum 
#number of samples by identifying the smallest population size in the file. 
#Our problem here is that there are populations (patches) with only one 
#individual. This would lead to a very imprecise evaluation of the Allelic 
#richness. Therefore, we are limiting the analysed dataset to patches with 
#at least a certain number (min_eff) of individuals
min_eff<-2

#Allelic richness of 2010 populations ##############
#formating dataset
geno2010<-genotype[genotype$survey=="MS10" & genotype$nb_snp_het==0 
                   & genotype$nb_missing<1 & genotype$nb_snp_typed==19 
                   & genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
#list of population with at least 2 individuals
higheff<-names(table(geno2010$patche_ID))[table(geno2010$patche_ID)>min_eff] 
#here we subset the file with the list of chosen patches
geno2010<-geno2010[geno2010$patche_ID %in% higheff,] 
geno2010<-drop.levels(geno2010)
#then we convert the file to an genind format of the adegenet package
geno2010ade<-replace(geno2010,geno2010=="AA",1)
geno2010ade<-replace(geno2010ade,geno2010ade=="CC",2)
geno2010ade<-replace(geno2010ade,geno2010ade=="GG",3)
geno2010ade<-replace(geno2010ade,geno2010ade=="TT",4)
loc_name<-colnames(geno2010ade[,3:(nb_SNP+3-1)])
data2010<-df2genind(geno2010ade[,loc_name],ploidy=1,
                    ind.names=geno2010ade$sample_ID,
                    loc.names=loc_name,pop=geno2010ade$patche_ID)
data2010@other$xy<-geno2010ade[,c("Longitude","Latitude")]
data2010@other$area<-geno2010ade[,c("Area_real")]
data2010pop<-genind2genpop(data2010,process.other=TRUE)

#Now we compute the allelic richness in each selected patches
AllRich(data2010)
#if we want to know the number of individual used in the rarefaction procedure 
AllRich(data2010)[1]
Ar2010<-AllRich(data2010)[[2]]
Ar2010<-apply(Ar2010,1,mean)
plot(Ar2010)
Ar2010<-cbind("patche_ID"=attr(Ar2010,"names"),"Ar2010"=Ar2010)
temp<-merge(Ar2010,presabs_index10,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Ar)))
boxplot(as.numeric(as.character(temp$Ar))~as.numeric(temp$cumulative_sum))

#computation of genotypic richness in 2010
data2010geno<-df2genind(cbind(geno2010ade[,loc_name],
                        as.matrix(as.numeric(as.factor(geno2010ade$MLG)))),
                        ploidy=1,ind.names=geno2010ade$sample_ID,
                        loc.names=c(loc_name,"MLG"),pop=geno2010ade$patche_ID)
data2010geno@other$xy<-geno2010ade[,c("Longitude","Latitude")]
data2010geno@other$area<-geno2010ade[,c("Area_real")]
Gr2010<-AllRich(data2010geno)[[2]]
Gr2010<-replace(Gr2010,Gr2010<1,1)
Gr2010<-cbind("patche_ID"=rownames(Gr2010),"Gr2010"=Gr2010[,"MLG"])
Gr2010
temp<-merge(Gr2010,presabs_index10,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Gr)))

#production of the final table for 2010
temp<-merge(Gr2010,Ar2010,by="patche_ID")
dimnames(coor2010)[[2]][1]<-"patche_ID"
Div_Patch2010<-merge(temp,nb_ind_MLG2010[nb_ind_MLG2010$number_genotyped!="0",],
                     by="patche_ID",all.y=TRUE)
Div_Patch2010<-cbind(Div_Patch2010,
              "S_0_1"=as.numeric(as.character(Div_Patch2010$connec2010)>0.005))
Div_Patch2010<-cbind(Div_Patch2010,
              "S_class"=as.numeric(as.character(Div_Patch2010$connec2010)>0.1))
Div_Patch2010$S_class[as.numeric(as.character(Div_Patch2010$connec2010)<0.005)]<-0
Div_Patch2010$S_class[as.numeric(as.character(Div_Patch2010$connec2010)<0.010) & 
                     as.numeric(as.character(Div_Patch2010$connec2010)>=0.005)]<-1
Div_Patch2010$S_class[as.numeric(as.character(Div_Patch2010$connec2010)<0.015) & 
                     as.numeric(as.character(Div_Patch2010$connec2010)>=0.010)]<-2
Div_Patch2010$S_class[as.numeric(as.character(Div_Patch2010$connec2010)<0.020) & 
                     as.numeric(as.character(Div_Patch2010$connec2010)>=0.015)]<-3
Div_Patch2010$S_class[as.numeric(as.character(Div_Patch2010$connec2010))>=0.025]<-4
#we replace NA values with 1, because the powdery mildew was at least present 
#the sampling year
Div_Patch2010[is.na(as.numeric(as.character(Div_Patch2010$cumulative_sum))),
              "cumulative_sum"]<-1
Div_Patch2010<-cbind(Div_Patch2010,
            "Age_0_1"=as.numeric(as.character(Div_Patch2010$cumulative_sum)>1))
#some examples of plot with the Div_Patch2010
boxplot(as.numeric(as.character(Div_Patch2010$Gr2010))~Div_Patch2010$S_0_1)
boxplot(as.numeric(as.character(Div_Patch2010$Gr2010))~Div_Patch2010$S_class)
boxplot(as.numeric(as.character(Div_Patch2010$Gr2010))~Div_Patch2010$Age_0_1)

#Allelic richness of 2011 populations ##############
#formating dataset
geno2011<-genotype[genotype$survey=="MS11" & genotype$nb_snp_het==0 & 
                     genotype$nb_missing<1 & genotype$nb_snp_typed==19 & 
                     genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
#list of population with more than 2 individuals
higheff<-names(table(geno2011$patche_ID))[table(geno2011$patche_ID)>min_eff]
#we remove these pop because for one SNP, they have a size of less than 3 
#(because of missing data)
higheff<-higheff[higheff!="817" & higheff!="3285"] 
geno2011<-geno2011[geno2011$patche_ID %in% higheff,]
geno2011<-drop.levels(geno2011)
geno2011ade<-replace(geno2011,geno2011=="AA",1)
geno2011ade<-replace(geno2011ade,geno2011ade=="CC",2)
geno2011ade<-replace(geno2011ade,geno2011ade=="GG",3)
geno2011ade<-replace(geno2011ade,geno2011ade=="TT",4)
loc_name<-colnames(geno2011ade[,3:(nb_SNP+3-1)])
data2011<-df2genind(geno2011ade[,loc_name],ploidy=1,
                    ind.names=geno2011ade$sample_ID,
                    loc.names=loc_name,pop=geno2011ade$patche_ID)
data2011@other$xy<-geno2011ade[,c("Longitude","Latitude")]
data2011@other$area<-geno2011ade[,c("Area_real")]
data2011pop<-genind2genpop(data2011,process.other=TRUE)

#Now we compute the allelic richness in each selected patches
AllRich(data2011)
#if we want to know the number of individual used in the rarefaction procedure
AllRich(data2011)[1]
Ar2011<-AllRich(data2011)[[2]]
Ar2011<-apply(Ar2011,1,mean)
plot(Ar2011)
Ar2011<-cbind("patche_ID"=attr(Ar2011,"names"),"Ar2011"=Ar2011)
temp<-merge(Ar2011,presabs_index11,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Ar)))
boxplot(as.numeric(as.character(temp$Ar))~as.numeric(temp$presence_sum))

#computation of genotypic richness in 2011
data2011geno<-df2genind(cbind(geno2011ade[,loc_name],
                        as.matrix(as.numeric(as.factor(geno2011ade$MLG)))),
                        ploidy=1,ind.names=geno2011ade$sample_ID,
                        loc.names=c(loc_name,"MLG"),pop=geno2011ade$patche_ID)
data2011geno@other$xy<-geno2011ade[,c("Longitude","Latitude")]
data2011geno@other$area<-geno2011ade[,c("Area_real")]
Gr2011<-AllRich(data2011geno)[[2]]
Gr2011<-replace(Gr2011,Gr2011<1,1)
Gr2011<-cbind("patche_ID"=rownames(Gr2011),"Gr2011"=Gr2011[,"MLG"])
Gr2011
temp<-merge(Gr2011,presabs_index11,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Gr)))

#production of the final table for 2011
temp<-merge(Gr2011,Ar2011,by="patche_ID")
dimnames(coor2011)[[2]][1]<-"patche_ID"
Div_Patch2011<-merge(temp,nb_ind_MLG2011[nb_ind_MLG2011$number_genotyped!="0",],
                     by="patche_ID",all.y=TRUE)
Div_Patch2011<-cbind(Div_Patch2011,
              "S_0_1"=as.numeric(as.character(Div_Patch2011$connec2011)>0.01))
Div_Patch2011<-cbind(Div_Patch2011,
             "S_class"=as.numeric(as.character(Div_Patch2011$connec2011)>0.01))
Div_Patch2011$S_class[as.numeric(as.character(Div_Patch2011$connec2011)<0.01)]<-0
Div_Patch2011$S_class[as.numeric(as.character(Div_Patch2011$connec2011)<0.02) & 
                   as.numeric(as.character(Div_Patch2011$connec2011)>=0.01)]<-1
Div_Patch2011$S_class[as.numeric(as.character(Div_Patch2011$connec2011)<0.03) & 
                   as.numeric(as.character(Div_Patch2011$connec2011)>=0.02)]<-2
Div_Patch2011$S_class[as.numeric(as.character(Div_Patch2011$connec2011)<0.04) & 
                   as.numeric(as.character(Div_Patch2011$connec2011)>=0.03)]<-3
Div_Patch2011$S_class[as.numeric(as.character(Div_Patch2011$connec2011))>=0.04]<-4
Div_Patch2011[is.na(as.numeric(as.character(Div_Patch2011$cumulative_sum))),
              "cumulative_sum"]<-1 
#we replace NA values with 1, because the powdery mildew was at least present 
#the sampling year
Div_Patch2011<-cbind(Div_Patch2011,
            "Age_0_1"=as.numeric(as.character(Div_Patch2011$cumulative_sum)>1))
#some examples of plot with the Div_Patch2011
boxplot(as.numeric(as.character(Div_Patch2011$Gr2011))~Div_Patch2011$S_0_1)
boxplot(as.numeric(as.character(Div_Patch2011$Gr2011))~Div_Patch2011$S_class)
boxplot(as.numeric(as.character(Div_Patch2011$Gr2011))~Div_Patch2011$Age_0_1)

#Allelic richness of 2012 populations ##############
#formating dataset
geno2012<-genotype[genotype$survey=="MS12" & genotype$nb_snp_het==0 & 
                     genotype$nb_missing<1 & genotype$nb_snp_typed==19 & 
                     genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
#list of population with more than 2 individuals
higheff<-names(table(geno2012$patche_ID))[table(geno2012$patche_ID)>min_eff]
#here we subset the file with the list of chosen patches
geno2012<-geno2012[geno2012$patche_ID %in% higheff,] 
geno2012<-drop.levels(geno2012)
#then we convert the file to an genind format of the adegenet package
geno2012ade<-replace(geno2012,geno2012=="AA",1)
geno2012ade<-replace(geno2012ade,geno2012ade=="CC",2)
geno2012ade<-replace(geno2012ade,geno2012ade=="GG",3)
geno2012ade<-replace(geno2012ade,geno2012ade=="TT",4)
loc_name<-colnames(geno2012ade[,3:(nb_SNP+3-1)])
data2012<-df2genind(geno2012ade[,loc_name],ploidy=1,
                    ind.names=geno2012ade$sample_ID,
                    loc.names=loc_name,pop=geno2012ade$patche_ID)
data2012@other$xy<-geno2012ade[,c("Longitude","Latitude")]
data2012@other$area<-geno2012ade[,c("Area_real")]
data2012pop<-genind2genpop(data2012,process.other=TRUE)

#Now we compute the allelic richness in each selected patches
AllRich(data2012)
#if we want to know the number of individual used in the rarefaction procedure
AllRich(data2012)[1]
Ar2012<-AllRich(data2012)[[2]]
Ar2012<-apply(Ar2012,1,mean)
plot(Ar2012)
Ar2012<-cbind("patche_ID"=attr(Ar2012,"names"),"Ar2012"=Ar2012)
plot(Ar2012)
temp<-merge(Ar2012,presabs_index12,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Ar)))
boxplot(as.numeric(as.character(temp$Ar))~as.numeric(temp$cumulative_sum))

#computation of genotypic richness in 2012
data2012geno<-df2genind(cbind(geno2012ade[,loc_name],
                        as.matrix(as.numeric(as.factor(geno2012ade$MLG)))),
                        ploidy=1,ind.names=geno2012ade$sample_ID,
                        loc.names=c(loc_name,"MLG"),pop=geno2012ade$patche_ID)
data2012geno@other$xy<-geno2012ade[,c("Longitude","Latitude")]
data2012geno@other$area<-geno2012ade[,c("Area_real")]
Gr2012<-AllRich(data2012geno)[[2]]
Gr2012<-cbind("patche_ID"=rownames(Gr2012),"Gr2012"=Gr2012[,"MLG"])
Gr2012
temp<-merge(Gr2012,presabs_index12,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Gr)))

#production of the final table for 2012
temp<-merge(Gr2012,Ar2012,by="patche_ID")
dimnames(coor2012)[[2]][1]<-"patche_ID"
Div_Patch2012<-merge(temp,nb_ind_MLG2012[nb_ind_MLG2012$number_genotyped!="0",],
                     by="patche_ID",all.y=TRUE)
Div_Patch2012<-cbind(Div_Patch2012,
              "S_0_1"=as.numeric(as.character(Div_Patch2012$connec2012)>0.005))
Div_Patch2012<-cbind(Div_Patch2012,
              "S_class"=as.numeric(as.character(Div_Patch2012$connec2012)>0.1))
Div_Patch2012$S_class[as.numeric(as.character(Div_Patch2012$connec2012)<0.005)]<-0
Div_Patch2012$S_class[as.numeric(as.character(Div_Patch2012$connec2012)<0.01) & 
                  as.numeric(as.character(Div_Patch2012$connec2012)>=0.005)]<-1
Div_Patch2012$S_class[as.numeric(as.character(Div_Patch2012$connec2012)<0.015) & 
                  as.numeric(as.character(Div_Patch2012$connec2012)>=0.010)]<-2
Div_Patch2012$S_class[as.numeric(as.character(Div_Patch2012$connec2012)<0.02) & 
                  as.numeric(as.character(Div_Patch2012$connec2012)>=0.015)]<-3
Div_Patch2012$S_class[as.numeric(as.character(Div_Patch2012$connec2012))>=0.025]<-4
Div_Patch2012[is.na(as.numeric(as.character(Div_Patch2012$cumulative_sum))),
              "cumulative_sum"]<-1 
#we replace NA values with 1, because the powdery mildew was at least present 
#the sampling year
Div_Patch2012<-cbind(Div_Patch2012,
            "Age_0_1"=as.numeric(as.character(Div_Patch2012$cumulative_sum)>1))
#some examples of plot with the Div_Patch2012
boxplot(as.numeric(as.character(Div_Patch2012$Gr2012))~Div_Patch2012$S_0_1)
boxplot(as.numeric(as.character(Div_Patch2012$Gr2012))~Div_Patch2012$S_class)
boxplot(as.numeric(as.character(Div_Patch2012$Gr2012))~Div_Patch2012$Age_0_1)

#Allelic richness of 2013 populations ##############
#formating dataset
geno2013<-genotype[genotype$survey=="MS13" & genotype$nb_snp_het==0 & 
                     genotype$nb_missing<1 & genotype$nb_snp_typed==19 & 
                     genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
#list of population with more than 2 individuals
higheff<-names(table(geno2013$patche_ID))[table(geno2013$patche_ID)>min_eff] 
#here we subset the file with the list of chosen patches
geno2013<-geno2013[geno2013$patche_ID %in% higheff,] 
geno2013<-drop.levels(geno2013)
#then we convert the file to an genind format of the adegenet package
geno2013ade<-replace(geno2013,geno2013=="AA",1)
geno2013ade<-replace(geno2013ade,geno2013ade=="CC",2)
geno2013ade<-replace(geno2013ade,geno2013ade=="GG",3)
geno2013ade<-replace(geno2013ade,geno2013ade=="TT",4)
loc_name<-colnames(geno2013ade[,3:(nb_SNP+3-1)])
data2013<-df2genind(geno2013ade[,loc_name],ploidy=1,
                    ind.names=geno2013ade$sample_ID,
                    loc.names=loc_name,pop=geno2013ade$patche_ID)
data2013@other$xy<-geno2013ade[,c("Longitude","Latitude")]
data2013@other$area<-geno2013ade[,c("Area_real")]
data2013pop<-genind2genpop(data2013,process.other=TRUE)

#Now we compute the allelic richness in each selected patches
AllRich(data2013)
#if we want to know the number of individual used in the rarefaction procedure
AllRich(data2013)[1]
Ar2013<-AllRich(data2013)[[2]]
Ar2013<-apply(Ar2013,1,mean)
plot(Ar2013)
Ar2013<-cbind("patche_ID"=attr(Ar2013,"names"),"Ar2013"=Ar2013)
plot(Ar2013)
temp<-merge(Ar2013,presabs_index13,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Ar)))
boxplot(as.numeric(as.character(temp$Ar))~as.numeric(temp$cumulative_sum))

#computation of genotypic richness in 2012
data2013geno<-df2genind(cbind(geno2013ade[,loc_name],
                        as.matrix(as.numeric(as.factor(geno2013ade$MLG)))),
                        ploidy=1,ind.names=geno2013ade$sample_ID,
                        loc.names=c(loc_name,"MLG"),pop=geno2013ade$patche_ID)
data2013geno@other$xy<-geno2013ade[,c("Longitude","Latitude")]
data2013geno@other$area<-geno2013ade[,c("Area_real")]
Gr2013<-AllRich(data2013geno)[[2]]
Gr2013<-cbind("patche_ID"=rownames(Gr2013),"Gr2013"=Gr2013[,"MLG"])
Gr2013
temp<-merge(Gr2013,presabs_index13,by="patche_ID")
plot(as.numeric(as.character(temp$cumulative_sum)),
     as.numeric(as.character(temp$Gr)))

#production of the final table for 2012
temp<-merge(Gr2013,Ar2013,by="patche_ID")
dimnames(coor2013)[[2]][1]<-"patche_ID"
Div_Patch2013<-merge(temp,nb_ind_MLG2013[nb_ind_MLG2013$number_genotyped!="0",]
                     ,by="patche_ID",all.y=TRUE)
Div_Patch2013<-cbind(Div_Patch2013,
              "S_0_1"=as.numeric(as.character(Div_Patch2013$connec2013)>0.005))
Div_Patch2013<-cbind(Div_Patch2013,
              "S_class"=as.numeric(as.character(Div_Patch2013$connec2013)>0.1))
Div_Patch2013$S_class[as.numeric(as.character(Div_Patch2013$connec2013)<0.005)]<-0
Div_Patch2013$S_class[as.numeric(as.character(Div_Patch2013$connec2013)<0.01) & 
                  as.numeric(as.character(Div_Patch2013$connec2013)>=0.005)]<-1
Div_Patch2013$S_class[as.numeric(as.character(Div_Patch2013$connec2013)<0.015) & 
                  as.numeric(as.character(Div_Patch2013$connec2013)>=0.010)]<-2
Div_Patch2013$S_class[as.numeric(as.character(Div_Patch2013$connec2013)<0.020) & 
                  as.numeric(as.character(Div_Patch2013$connec2013)>=0.015)]<-3
Div_Patch2013$S_class[as.numeric(as.character(Div_Patch2013$connec2013))>=0.025]<-4
Div_Patch2013[is.na(as.numeric(as.character(Div_Patch2013$cumulative_sum))),
              "cumulative_sum"]<-1 
#we replace NA values with 1, because the powdery mildew was at least present 
#the sampling year
Div_Patch2013<-cbind(Div_Patch2013,
            "Age_0_1"=as.numeric(as.character(Div_Patch2013$cumulative_sum)>1))
#some examples of plot with the Div_Patch2013
boxplot(as.numeric(as.character(Div_Patch2013$Gr2013))~Div_Patch2013$S_0_1)
boxplot(as.numeric(as.character(Div_Patch2013$Gr2013))~Div_Patch2013$S_class)
boxplot(as.numeric(as.character(Div_Patch2013$Gr2013))~Div_Patch2013$Age_0_1)


###############################################################################
#export the tables and clean the environment
###############################################################################

write.table(Div_Patch2010,file="data/stat_patch2010.txt",row.names=FALSE,
            quote=FALSE,sep="\t")
write.table(Div_Patch2011,file="data/stat_patch2011.txt",row.names=FALSE,
            quote=FALSE,sep="\t")
write.table(Div_Patch2012,file="data/stat_patch2012.txt",row.names=FALSE,
            quote=FALSE,sep="\t")
write.table(Div_Patch2013,file="data/stat_patch2013.txt",row.names=FALSE,
            quote=FALSE,sep="\t")

rm(list=ls())


###############################################################################
#END
###############################################################################