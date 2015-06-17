###############################################################################
###############################################################################
#Building the dataset per commuune and SIN_86 for LD computation
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first

setwd("~/work/Rfichiers/Githuber/SelfOut_data")

#In order to compute Linkage Disequilibrium (LD), we use the GENEPOP on the web 
#software. The computation of the LD is done at the Aland level as well as 
#two different sublevels: commune and semi-independant network (SIN). Because of 
#the large number of these sublevels, we have prepared the dataset using the 
#following code, which produces files that can be easily used as input for 
#GENEPOP on the web after a few editions

library(plyr)
library(gdata)
#Let's have a fresh start and load again the original genotype information 
genotype<-geno_hom

#formating dataset: 2010####
geno2010ade<-genotype[genotype$survey=="MS10" & genotype$nb_snp_het==0 
                   & genotype$nb_missing<1 & genotype$nb_snp_typed==19 
                   & genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2010ade<-geno2010ade[!is.na(geno2010ade$patche_ID),]
geno2010ade<-drop.levels(geno2010ade)
#we convert genotypes wuth letters into genotypes with number
geno2010ade<-replace(geno2010ade,geno2010ade=="AA",1)
geno2010ade<-replace(geno2010ade,geno2010ade=="CC",2)
geno2010ade<-replace(geno2010ade,geno2010ade=="GG",3)
geno2010ade<-replace(geno2010ade,geno2010ade=="TT",4)
#we then identify the different MLG in the different sublevel
geno2010Comm<-ddply(geno2010ade,.(Commune,MLG),summarise,nb_geno=length(MLG))
geno2010SIN86<-ddply(geno2010ade,.(SIN_86,MLG),summarise,nb_geno=length(MLG))
#exporting the file to a text format
write.table(geno2010Comm,file="geno2010Comm.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
write.table(geno2010SIN86,file="geno2010SIN86.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
#cleaning the environment
rm(geno2010Comm,geno2010SIN86,geno2010ade)

#formating dataset: 2011####
geno2011ade<-genotype[genotype$survey=="MS11" & genotype$nb_snp_het==0 
                      & genotype$nb_missing<1 & genotype$nb_snp_typed==19 
                      & genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2011ade<-geno2011ade[!is.na(geno2011ade$patche_ID),]
geno2011ade<-drop.levels(geno2011ade)
#we convert genotypes wuth letters into genotypes with number
geno2011ade<-replace(geno2011ade,geno2011ade=="AA",1)
geno2011ade<-replace(geno2011ade,geno2011ade=="CC",2)
geno2011ade<-replace(geno2011ade,geno2011ade=="GG",3)
geno2011ade<-replace(geno2011ade,geno2011ade=="TT",4)
#we then identify the different MLG in the different sublevel
geno2011Comm<-ddply(geno2011ade,.(Commune,MLG),summarise,nb_geno=length(MLG))
geno2011SIN86<-ddply(geno2011ade,.(SIN_86,MLG),summarise,nb_geno=length(MLG))
#exporting the file to a text format
write.table(geno2011Comm,file="geno2011Comm.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
write.table(geno2011SIN86,file="geno2011SIN86.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
#cleaning the environment
rm(geno2011Comm,geno2011SIN86,geno2011ade)

#formating dataset: 2012####
geno2012ade<-genotype[genotype$survey=="MS12" & genotype$nb_snp_het==0 
                      & genotype$nb_missing<1 & genotype$nb_snp_typed==19 
                      & genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2012ade<-geno2012ade[!is.na(geno2012ade$patche_ID),]
geno2012ade<-drop.levels(geno2012ade)
#we convert genotypes wuth letters into genotypes with number
geno2012ade<-replace(geno2012ade,geno2012ade=="AA",1)
geno2012ade<-replace(geno2012ade,geno2012ade=="CC",2)
geno2012ade<-replace(geno2012ade,geno2012ade=="GG",3)
geno2012ade<-replace(geno2012ade,geno2012ade=="TT",4)
#we then identify the different MLG in the different sublevel
geno2012Comm<-ddply(geno2012ade,.(Commune,MLG),summarise,nb_geno=length(MLG))
geno2012SIN86<-ddply(geno2012ade,.(SIN_86,MLG),summarise,nb_geno=length(MLG))
#exporting the file to a text format
write.table(geno2012Comm,file="geno2012Comm.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
write.table(geno2012SIN86,file="geno2012SIN86.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
#cleaning the environment
rm(geno2012Comm,geno2012SIN86,geno2012ade)

#formating dataset: 2013####
geno2013ade<-genotype[genotype$survey=="MS13" & genotype$nb_snp_het==0 
                      & genotype$nb_missing<1 & genotype$nb_snp_typed==19 
                      & genotype$duplicate!="DUP" & genotype$leave_ID=="a",]
geno2013ade<-geno2013ade[!is.na(geno2013ade$patche_ID),]
geno2013ade<-drop.levels(geno2013ade)
#we convert genotypes wuth letters into genotypes with number
geno2013ade<-replace(geno2013ade,geno2013ade=="AA",1)
geno2013ade<-replace(geno2013ade,geno2013ade=="CC",2)
geno2013ade<-replace(geno2013ade,geno2013ade=="GG",3)
geno2013ade<-replace(geno2013ade,geno2013ade=="TT",4)
#we then identify the different MLG in the different sublevel
geno2013Comm<-ddply(geno2013ade,.(Commune,MLG),summarise,nb_geno=length(MLG))
geno2013SIN86<-ddply(geno2013ade,.(SIN_86,MLG),summarise,nb_geno=length(MLG))
#exporting the file to a text format
write.table(geno2013Comm,file="geno2013Comm.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
write.table(geno2013SIN86,file="geno2013SIN86.txt",sep="\t",
            quote=FALSE,row.names=FALSE)
#cleaning the environment
rm(geno2013Comm,geno2013SIN86,geno2013ade)
