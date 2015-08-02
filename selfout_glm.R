###############################################################################
###############################################################################
#Explore effect of coinfection on the disease dynamic
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")

library(lme4)
library(MASS)
library(car)
library(visreg)

#this code contains the details of the models used for the analysis of the 
#impact of coinfection on the overwintering success, the production of the 
#overwintering structures (chasmotecia) and the production of new Multilocus 
#Genotypes (MLG) the following year


###############################################################################
#effect of coinfection on the overwintering of the disease at the patch level
###############################################################################

#we first investigate the effect of the binary variable for presence/absence 
#of coinfection the previous year on overwintering success

#2013 preparing the datafile (coinfection 0/1)
coinf<-coinf2013
#patches with coinfected sampled should have at list 2 different MLG, 
#we therefore apply a correction to the numnber of MLG for patches 
#with coinfection but with less than 2 single identified MLG
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable: with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-coinf[!is.na(coinf$AA_F2013),]
#we test simple models taking into account the presence/absence of coinfection 
#and the absolute abundance of the disease before the winter
OvWi13_1<-glm(PA_S2014~coinfYN+AA_F2013,family=binomial,data=coinf)
summary(OvWi13_1)
OvWi13_2<-glm(PA_S2014~AA_F2013,family=binomial,data=coinf)
summary(OvWi13_2)
OvWi13_3<-glm(PA_S2014~coinfYN,family=binomial,data=coinf)
summary(OvWi13_3)
anova(OvWi13_2,OvWi13_1,test="Chisq")
summary(OvWi13_1)
drop1(OvWi13_1,test="Chisq")
#keeping only the model of interest
rm(OvWi13_2,OvWi13_3)
#2013 an example of plot
op<-par(mfrow=c(1,2))
visreg(OvWi13_1,"coinfYN",scale="response",overlay=TRUE,
       xlab="No coinf (0)/ coinf (1)",ylab="P(P/A 2014)")
visreg(OvWi13_1,"AA_F2013",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Abundance Fall",ylab="P(P/A 2014)")
par(op)
#2013 other visualisation
plot(coinf$PA_S2014~coinf$coinfYN)

#2012 preparing the datafile (coinfection 0/1)
coinf<-coinf2012
#patches with coinfected sampled should have at list 2 different MLG, 
#we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-coinf[!is.na(coinf$AA_F2012),]
#we test simple models taking into account the presence/absence of coinfection 
#and the absolute abundance of the disease before the winter
OvWi12_1<-glm(PA_S2013~coinfYN+AA_F2012,family=binomial,data=coinf)
summary(OvWi12_1)
OvWi12_2<-glm(PA_S2013~AA_F2012,family=binomial,data=coinf)
summary(OvWi12_2)
OvWi12_3<-glm(PA_S2013~coinfYN,family=binomial,data=coinf)
summary(OvWi12_3)
anova(OvWi12_2,OvWi12_1,test="Chisq")
summary(OvWi12_1)
drop1(OvWi12_1,test="Chisq")
#keeping only the model of interest
rm(OvWi12_2,OvWi12_3)
#2012 an example of plot
op<-par(mfrow=c(1,2))
visreg(OvWi12_1,"coinfYN",scale="response",overlay=TRUE,
       xlab="No coinf (0)/ coinf (1)",ylab="P(P/A 2013)")
visreg(OvWi12_1,"AA_F2012",rug=2,scale="response",jitter=TRUE,by="coinfYN",
       overlay=TRUE,partial=FALSE,xlab="Abundance Fall",ylab="P(P/A 2013)")
par(op)
#2012 other visualisation
plot(coinf$PA_S2013~coinf$coinfYN)


#we then investigate the effect of the percentage of coinfection on the  
#overwintering success. In order to do that, we limit the dataset to 
#patches with at least 3 individuals, so that the estimation of the percentage 
#of coinfection is relevant

#2013 preparing the datafile (coinfection %)
coinf<-coinf2013
#patches with coinfected sampled should have at list 2 different MLG, 
#we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)
coinf<-coinf[!is.na(coinf$AA_F2013),]
#we select patches with at least 3 samples
coinf<-coinf[coinf$number_genotyped>2,]
#we test simple models taking into account the percentage of coinfection 
#and the absolute abundance of the disease before the winter
pOvWi13_1<-glm(PA_S2014~percoinf+AA_F2013,family=binomial,data=coinf)
summary(pOvWi13_1)
pOvWi13_2<-glm(PA_S2014~AA_F2013,family=binomial,data=coinf)
summary(pOvWi13_2)
pOvWi13_3<-glm(PA_S2014~percoinf,family=binomial,data=coinf)
summary(pOvWi13_3)
anova(pOvWi13_2,pOvWi13_1,test="Chisq")
summary(pOvWi13_1)
drop1(pOvWi13_1,test="Chisq")
#keeping only the model of interest
rm(pOvWi13_2,pOvWi13_3)

#2012 preparing the datafile (coinfection %)
coinf<-coinf2012
#patches with coinfected sampled should have at list 2 different MLG, 
#we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)
coinf<-coinf[!is.na(coinf$AA_F2012),]
#we select patches with at least 3 samples
coinf<-coinf[coinf$number_genotyped>2,]
#we test simple models taking into account the percentage of coinfection 
#and the absolute abundance of the disease before the winter
pOvWi12_1<-glm(PA_S2013~percoinf+AA_F2012,family=binomial,data=coinf)
summary(pOvWi12_1)
pOvWi12_2<-glm(PA_S2013~AA_F2012,family=binomial,data=coinf)
summary(pOvWi12_2)
pOvWi12_3<-glm(PA_S2013~percoinf,family=binomial,data=coinf)
summary(pOvWi12_3)
anova(pOvWi12_2,pOvWi12_1,test="Chisq")
summary(pOvWi12_1)
drop1(pOvWi12_1,test="Chisq")
#keeping only the model of interest
rm(pOvWi12_2,pOvWi12_3)


###############################################################################
#Effect of coinfection on the production of chasmotecia
###############################################################################

#P/A of chasmo against P/A coinf####

#2013 with P/A chasmo and P/A coinf#######
coinf<-coinf2013
#patches with coinfected sampled should have at list 2 different MLG, 
#we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-data.frame(coinf,"chasmoYN"=coinf$chasmo_f_2013)
coinf$chasmoYN[(coinf$chasmoYN)>0]<-1
coinf$chasmoYN<-as.factor(coinf$chasmoYN)
coinf<-coinf[!is.na(coinf$AA_F2013),]
#we build models to test if the presence/absence of coinfection have an 
#impact on the presence/absence of chasmotecia
chasmo13_1<-glm(chasmoYN~coinfYN+AA_F2013,family=binomial,data=coinf)
summary(chasmo13_1)
chasmo13_2<-glm(chasmoYN~AA_F2013,family=binomial,data=coinf)
summary(chasmo13_2)
chasmo13_3<-glm(chasmoYN~coinfYN,family=binomial,data=coinf)
summary(chasmo13_3)
anova(chasmo13_2,chasmo13_1,test="Chisq")
#keeping only the model of interest
rm(chasmo13_2,chasmo13_3)

#2012 with P/A chasmo and P/A coinf#######
coinf<-coinf2012
#patches with coinfected sampled should have at list 2 different MLG, 
#we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-data.frame(coinf,"chasmoYN"=coinf$chasmo_f_2012)
coinf$chasmoYN[(coinf$chasmoYN)>0]<-1
coinf$chasmoYN<-as.factor(coinf$chasmoYN)
coinf<-coinf[!is.na(coinf$AA_F2012),]
coinf<-coinf[!is.na(coinf$chasmoYN),]
#we build models to test if the presence/absence of coinfection have an 
#impact on the presence/absence of chasmotecia
chasmo12_1<-glm(chasmoYN~coinfYN+AA_F2012,family=binomial,data=coinf)
summary(chasmo12_1)
chasmo12_2<-glm(chasmoYN~AA_F2012,family=binomial,data=coinf)
summary(chasmo12_2)
chasmo12_3<-glm(chasmoYN~coinfYN,family=binomial,data=coinf)
summary(chasmo12_3)
anova(chasmo12_2,chasmo12_1,test="Chisq")
rm(chasmo12_2,chasmo12_3)


###############################################################################
#Effect of coinfection on the production of new genotypes the following year
###############################################################################

#we first merge the data of 2012 and 2013
datanewgeno<-merge(coinf2012,coinf2013,by="patche_ID")
datanewgeno<-datanewgeno[datanewgeno$number_genotyped.x>2,]
datanewgeno<-cbind(datanewgeno,"percoinf12"=datanewgeno$number_coinf.x*100/
                     datanewgeno$number_genotyped.x)

#we start to explore the model with the following explanatory variables: the 
#percentage of coinfection of the previous year, the prevalence of chasmothecia 
#in the previous year and the pathogen connectivity in the previous year. 
#the response response variable is the number of MLG observed in the patch 
#that were not present in the previous year
modnewgeno<-glm(number_pure_new.y~percoinf12+chasmoplant2012.x+connec2012.x,
                  family=poisson,data=datanewgeno)
summary(modnewgeno)
modnewgeno<-glm(number_pure_new.y~percoinf12+connec2012.x,
                family=poisson,data=datanewgeno)
summary(modnewgeno)

#we just keep the percentage of coinfection and the connectivity in the model
modnewgeno_1<-glm(number_pure_new.y~percoinf12+connec2012.x,
                family=poisson,data=datanewgeno)
summary(modnewgeno_1)
modnewgeno_2<-glm(number_pure_new.y~percoinf12,
                  family=poisson,data=datanewgeno)
summary(modnewgeno_2)
modnewgeno_3<-glm(number_pure_new.y~connec2012.x,
                  family=poisson,data=datanewgeno)
summary(modnewgeno_3)
anova(modnewgeno_3,modnewgeno_1,test="Chisq")
anova(modnewgeno_2,modnewgeno_1,test="Chisq")

#a plot of the variation of the number of new MLG per patch as a function of 
#the percentage of coinfection in the same patch the previous year
visreg(modnewgeno_1,"percoinf12",scale="response",lwd=10,font.lab=2,font=2,
       font.axis=2,xlab="% of coinfection in 2012",line.par=list(col="black"),
       ylab="Number of new MLG in 2013",bty="l",rug=FALSE,xlim=c(0,100),axes=FALSE)
box(bty="l",lwd=3)
axis(1,lwd=3,font=2)
axis(2,lwd=3,font=2,las=1)
rug(jitter(datanewgeno$percoinf12,amount=1),col="violetred")
#transparent colors don't seem to work with 'rug'. In case this is fixed at 
#some point, here is an example of transparent color: 
#rgb(208,32,144,alpha=0.3,maxColorValue=255)

#export to a pdf file of 12 x 8 inches

#another figure
breaks<-c(-1,10,20,30,40,50)
freq.cut<-cut(datanewgeno$percoinf12,breaks)
boxplot(datanewgeno$number_pure_new.y~freq.cut)

#cleaning the environment
rm(breaks,freq.cut,modnewgeno_2,modnewgeno_3)


###############################################################################
#END
###############################################################################