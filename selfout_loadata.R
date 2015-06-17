###############################################################################
###############################################################################
#Load the data for analysis
###############################################################################
###############################################################################

setwd("~/work/Rfichiers/Githuber/SelfOut_data")
#before running this code, you need to import and transform the dataset using 
#'selfout_importdata.R' first


###############################################################################
#loading the different datasets
###############################################################################

#the code to deal with shapefile instead of RData file for GIS purpose
Aland<-readShapePoly("RECNO.shp",proj4string=CRS("+init=epsg:2393"))
Aland<-spTransform(Aland,CRS("+init=epsg:3067"))
patchshape<-readShapePoly("All patches.shp",proj4string=CRS("+init=epsg:3067"))

#loading patches informations
patche_info<-read.table("coord_all_patch12.txt",header=TRUE, sep="\t",dec=".")
#we reorganize the coord data file
patche_info<-patche_info[,c(1:12,14:dim(patche_info)[2])]

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
rm(i,MLGmat,MLG)

#we load the datatable which combine patches, genotypes information and 
#indices that has been prepared using selfout_importdata.R
coinf2010<-read.table("stat_patch2010.txt",header=TRUE,sep="\t",
                      colClasses=c("factor","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","factor",
                                   "character","character","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","numeric","factor",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","factor",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric",
                                   "factor","factor","factor")
)


coinf2011<-read.table("stat_patch2011.txt",header=TRUE,sep="\t",
                      colClasses=c("factor","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric","factor",
                                   "character","character","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","numeric","factor",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","factor",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric",
                                   "factor","factor","factor")
)


coinf2012<-read.table("stat_patch2012.txt",header=TRUE,sep="\t",
                      colClasses=c("factor","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric","factor",
                                   "character","character","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","numeric","factor",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","factor",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric",
                                   "factor","factor","factor")
)

coinf2013<-read.table("stat_patch2013.txt",header=TRUE,sep="\t",
                      colClasses=c("factor","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric","factor",
                                   "character","character","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","factor","factor",
                                   "factor","factor","factor","factor","numeric","factor",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","factor","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","factor",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric","numeric",
                                   "numeric","factor","numeric","numeric",
                                   "numeric","numeric","numeric","numeric","numeric",
                                   "factor","factor","factor")
)



coinf<-coinf2012 #here you can change the file you want to use for the analysis

#adding a percentage of coinfection column to the data table 
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)

#modification of the number of MLG in case of coinfection
coinfcorr<-coinf
coinfcorr[coinfcorr$number_coinf>0 & coinfcorr$number_MLG<2,]$number_MLG<-2

#dataset for P/A of coinfection instead of prevalence of coinfection
coinfYN<-coinf
coinfYN[coinfYN$number_coinf>0,]$number_coinf<-1

#dataset with at least four sample for overwintering of the population analysis
coinf4corr<-coinfcorr
coinf4corr<-coinf4corr[coinf4corr$number_genotyped>3,]

#dataset for impact of coinfection on prevalence of the disease (using evolution of 
#prevalence between spring and fall of relative abundance of powdery mildew)
evolinf<-coinf
evolinf<-evolinf[evolinf$RA_S2012>0,]
evolinf<-evolinf[!is.na(evolinf$RA_S2012),]

#dataset for impact of coinfection on prevalence of the disease (using evolution of 
#prevalence between spring and fall of relative abundance of powdery mildew) with at leat 4 samples
evolinf4<-coinfcorr
evolinf4<-evolinf4[evolinf4$RA_S2012>0,]
evolinf4<-evolinf4[evolinf4$number_genotyped>3,]
evolinf4<-evolinf4[!is.na(evolinf4$AA_S2012),]
evolinf4<-evolinf4[!is.na(evolinf4$AA_F2012),]
evolinf4<-cbind(evolinf4,"YNcoinf"=cut(evolinf4$percoinf,breaks=100*(-1:1),labels=c(0,1)))
evolinf4<-cbind(evolinf4,"coinfCAT"=cut(evolinf4$percoinf,breaks=c(-1,1,24,49,74,100),labels=c(0,1,2,3,4),
                                        ordered_result=TRUE))


##########################################################################################
#formating dataset for production of chasmotecia
##########################################################################################


chasmo<-read.table("chasmo.txt",header=TRUE,sep="\t",
                      colClasses=c("factor","numeric","factor","numeric","numeric","numeric",
                                   "factor","numeric","factor")
)

