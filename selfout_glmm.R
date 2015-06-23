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
#on overwintering success. In order to do that, we limit the dataset to 
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
#keeping only the model of interest
rm(pOvWi12_2,pOvWi12_3)


###############################################################################
#Effect of coinfection on the production of chasmotecia
###############################################################################

#P/A of chasmo against P/A coinf####

#2013 with P/A chasmo and P/A coinf#######
coinf<-coinf2013
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
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
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
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

chasmo12<-glmer(chasmoYN~coinfYN+AA_F2012+Gdiv+connec2012+cumulative_sum+(1|SIN_86),family=binomial,data=coinf)
summary(chasmo12) #problem due to the random factor which capture all the variance
chasmo12<-glm(chasmoYN~coinfYN+AA_F2012+Gdiv+connec2012+cumulative_sum,family=binomial,data=coinf)
summary(chasmo12)
chasmo12_1<-glm(chasmoYN~coinfYN+AA_F2012,family=binomial,data=coinf)
summary(chasmo12_1)
chasmo12_2<-glm(chasmoYN~AA_F2012,family=binomial,data=coinf)
summary(chasmo12_2)
chasmo12_3<-glm(chasmoYN~coinfYN,family=binomial,data=coinf)
summary(chasmo12_3)
anova(chasmo12_2,chasmo12_1,test="Chisq")

#level of chasmothecia against the percentage of coinfection####

#2013#######
coinf<-coinf2013
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)
coinf<-coinf[!is.na(coinf$AA_F2013),]
coinf<-coinf[coinf$number_genotyped>2,]

chasmo13<-glmer(chasmo_f_2013~percoinf+AA_F2013+Gdiv+(1|SIN_86),family=poisson,data=coinf)
summary(chasmo13)
chasmo13<-glm(chasmoplant2013~coinfYN+AA_F2013,family=gaussian,data=coinf)
summary(chasmo13)

#2012#######
coinf<-coinf2012
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)
coinf<-coinf[!is.na(coinf$AA_F2012),]
coinf<-coinf[coinf$number_genotyped>2,]

chasmo12<-glmer(chasmo_f_2012~percoinf+AA_F2012+Gdiv+(1|SIN_86),family=poisson,data=coinf)
summary(chasmo12)


###############################################################################
#Is the number of MLG affected by connectivity?
###############################################################################

plot(table(coinf4corr$number_MLG))
nbMLG<-glm(number_MLG~percoinf+connec2012+Area_real+cumulative_sum,family=poisson,data=coinf4corr)
summary(nbMLG)
nbMLG<-glm(number_MLG~percoinf+Area_real,family=poisson,data=coinf4corr)
summary(nbMLG)
visreg(nbMLG,"Area_real",scale="response",jitter=TRUE,
       xlab="Patche's Area",ylab="Number MLG")
visreg(nbMLG,"percoinf",scale="response",jitter=TRUE,
       xlab="Percentage of coinfection",ylab="Number MLG")

plot(table(coinf4corr$number_pure_new))
nbnewMLG<-glm(number_pure_new~percoinf+Area_real+connec2012+cumulative_sum,family=poisson,data=coinf4corr)
summary(nbnewMLG)
nbnewMLG<-glm(number_pure_new~percoinf+Area_real+cumulative_sum,family=poisson,data=coinf4corr)
summary(nbnewMLG)
visreg(nbnewMLG,"Area_real",scale="response",jitter=TRUE,
       xlab="Patche's Area",ylab="Number MLG")
visreg(nbnewMLG,"percoinf",scale="response",jitter=TRUE,
       xlab="Percentage of coinfection",ylab="Number MLG")
visreg(nbnewMLG,"cumulative_sum",scale="response",jitter=TRUE,
       xlab="Consecutive years of mildew",ylab="Number MLG")


#is this the right way of doing it: 
nbMLG<-glm(number_MLG~percoinf+connec2012+Area_real+cumulative_sum+offset(log(number_genotyped)) ,family=poisson,data=coinf4corr)
summary(nbMLG)




###############################################################################
#Effect of coinfection on the production of new genotypes the following year
###############################################################################

datanewgeno<-merge(coinf2012,coinf2013,by="patche_ID")
datanewgeno<-datanewgeno[datanewgeno$number_genotyped.x>2,]
datanewgeno<-cbind(datanewgeno,"percoinf12"=datanewgeno$number_coinf.x*100/datanewgeno$number_genotyped.x)
datanewgeno2<-datanewgeno[!is.na(datanewgeno$coinr2012),]

modnewgeno<-glm(number_pure_new.y~percoinf12+chasmoplant2012.x+connec2012.x+coinr2012,
                family=poisson,data=datanewgeno)
summary(modnewgeno)
modnewgeno<-glm(number_pure_new.y~percoinf12+connec2012.x+coinr2012,
                family=poisson,data=datanewgeno)
summary(modnewgeno)
modnewgeno<-glm(number_pure_new.y~percoinf12+coinr2012,
                family=poisson,data=datanewgeno)
summary(modnewgeno)

modnewgeno_1<-glm(number_pure_new.y~percoinf12+connec2012.x,
                family=poisson,data=datanewgeno)
summary(modnewgeno_1)
modnewgeno_2<-glm(number_pure_new.y~percoinf12,
                  family=poisson,data=datanewgeno)
summary(modnewgeno_2)
modnewgeno_3<-glm(number_pure_new.y~connec2012.x,
                  family=poisson,data=datanewgeno)
summary(modnewgeno_3)
anova(modnewgeno_2,modnewgeno_1,test="Chisq")


#New MLG define at the SIN level (instead of the metapopulation level)
modnewgeno<-glm(number_pure_new_SIN.y~percoinf12,
                family=poisson,data=datanewgeno)
summary(modnewgeno)

visreg(modnewgeno,"percoinf12",scale="response",jitter=TRUE,
       xlab="Percentage of coinfection",ylab="New MLG")

breaks<-c(-1,10,20,30,40,50)
freq.cut<-cut(datanewgeno$percoinf12,breaks)
boxplot(datanewgeno$number_pure_new.y~freq.cut)


###############################################################################
#Effect of new MLG on the evolution of intra-epidemic dynamic
###############################################################################

model.ratioperc<-lm((AA_F2013-AA_S2013)~number_pure_new+connec2013+percoinf,data=evolinf4)
summary(model.ratioperc)
anova(model.ratioperc)
model.ratioperc<-lm((AA_F2013)~number_pure_new+connec2013+percoinf,data=evolinf4)
summary(model.ratioperc)

###############################################################################
#factors affecting coinfection in space
###############################################################################

#2013#######
coinf<-coinf2013
#patches with coinfected sampled should have at least 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)

prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2013+AA_F2013,family=binomial,data=coinf)
summary(prevalcoinf)
#if we replace the raw number of genotype by the Genotypic diversity
prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2013+AA_F2013,family=binomial,data=coinf)
summary(prevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2013+AA_F2013+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2013+AA_F2013+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)

op<-par(mfrow=c(1,4))
visreg(prevalcoinf,"cumulative_sum",scale="response",overlay=TRUE,
       xlab="Population age",ylab="Proba Coinfection")
visreg(prevalcoinf,"AA_F2013",rug=2,scale="response",jitter=TRUE,
       overlay=TRUE,partial=FALSE,xlab="Abundance Fall",ylab="Proba Coinfection")
visreg(prevalcoinf,"Gdiv",rug=2,scale="response",jitter=TRUE,
       overlay=TRUE,partial=FALSE,xlab="Genotypic diversity",ylab="Proba Coinfection")
visreg(prevalcoinf,"connec2013",rug=2,scale="response",jitter=TRUE,
       overlay=TRUE,partial=FALSE,xlab="Connectivity",ylab="Proba Coinfection")
par(op)


#2012#######
coinf<-coinf2012
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)

prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2012+AA_F2012,family=binomial,data=coinf)
summary(prevalcoinf)
#if we replace the raw number of genotype by the Genotypic diversity
prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2012+AA_F2012,family=binomial,data=coinf)
summary(prevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2012+AA_F2012+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2012+AA_F2012+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)


#2011#######
coinf<-coinf2011
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)

prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2011+AA_F2011,family=binomial,data=coinf)
summary(prevalcoinf)
#if we replace the raw number of genotype by the Genotypic diversity
prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2011+AA_F2011,family=binomial,data=coinf)
summary(prevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2011+AA_F2011+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2011+AA_F2011+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)


#2010#######
coinf<-coinf2010
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)

prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2010+AA_F2010,family=binomial,data=coinf)
summary(prevalcoinf)
#if we replace the raw number of genotype by the Genotypic diversity
prevalcoinf<-glm(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2010+AA_F2010,family=binomial,data=coinf)
summary(prevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+number_MLG+connec2010+AA_F2010+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)
mixprevalcoinf<-glmer(cbind(number_coinf,number_genotyped-number_coinf)~cumulative_sum+Gdiv+connec2010+AA_F2010+(1|SIN_86),family=binomial,data=coinf)
summary(mixprevalcoinf)


###############################################################################
#effect of coinfection on the prevalence of the disease and dynamic of the epidemic
###############################################################################

#####evolution of disease abundance 2012 and 2013 ####

#2013#######
coinf<-coinf2013
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
#if we want to limit to populations with at least 3 individuals genotyped
coinf<-coinf[coinf$RA_S2013>0,]
coinf<-coinf[coinf$number_genotyped>2,]
coinf<-coinf[!is.na(coinf$AA_S2013),]
coinf<-coinf[!is.na(coinf$AA_F2013),]
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)

#variation of the absolute abundance between spring and fall with random SIN effect
DynAA13.S<-lmer(AA_F2013-AA_S2013~connec2013+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(DynAA13.S)
DynAA13.coin<-lmer(AA_F2013-AA_S2013~percoinf+(1|SIN_86),
                   data=coinf,REML=FALSE)
summary(DynAA13.coin)
DynAA13.gdiv<-lmer(AA_F2013-AA_S2013~Gdiv+(1|SIN_86),
                   data=coinf,REML=FALSE)
summary(DynAA13.gdiv)
DynAA13.coinS<-lmer(AA_F2013-AA_S2013~percoinf+connec2013+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(DynAA13.coinS)
DynAA13.coingdiv<-lmer(AA_F2013-AA_S2013~percoinf+Gdiv+(1|SIN_86),
                       data=coinf,REML=FALSE)
summary(DynAA13.coingdiv)
DynAA13.Sgdiv<-lmer(AA_F2013-AA_S2013~connec2013+Gdiv+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(DynAA13.Sgdiv)
DynAA13.full<-lmer(AA_F2013-AA_S2013~percoinf*Gdiv+(1|SIN_86),
                   data=coinf,REML=FALSE)
summary(DynAA13.full)

anova(DynAA13.coinS,DynAA13.S)
anova(DynAA13.coinS,DynAA13.coin)
anova(DynAA13.Sgdiv,DynAA13.S)
anova(DynAA13.Sgdiv,DynAA13.gdiv)
anova(DynAA13.coingdiv,DynAA13.gdiv)
anova(DynAA13.coingdiv,DynAA13.coin)
anova(DynAA13.full)
anova(DynAA13.coingdiv,DynAA13.full)


#2012#######
coinf<-coinf2012
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
#if we want to limit to populations with at least 3 individuals genotyped
coinf<-coinf[coinf$RA_S2012>0,]
coinf<-coinf[coinf$number_genotyped>2,]
coinf<-coinf[!is.na(coinf$AA_S2012),]
coinf<-coinf[!is.na(coinf$AA_F2012),]
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)

#variation of the absolute abundance between spring and fall with random SIN effect
DynAA12.S<-lmer(AA_F2012-AA_S2012~connec2012+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(DynAA12.S)
DynAA12.coin<-lmer(AA_F2012-AA_S2012~percoinf+(1|SIN_86),
                   data=coinf,REML=FALSE)
summary(DynAA12.coin)
DynAA12.gdiv<-lmer(AA_F2012-AA_S2012~Gdiv+(1|SIN_86),
                   data=coinf,REML=FALSE)
summary(DynAA12.gdiv)
DynAA12.coinS<-lmer(AA_F2012-AA_S2012~percoinf+connec2012+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(DynAA12.coinS)
DynAA12.coingdiv<-lmer(AA_F2012-AA_S2012~percoinf+Gdiv+(1|SIN_86),
                       data=coinf,REML=FALSE)
summary(DynAA12.coingdiv)
DynAA12.Sgdiv<-lmer(AA_F2012-AA_S2012~connec2012+Gdiv+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(DynAA12.Sgdiv)
DynAA12.full<-lmer(AA_F2012-AA_S2012~percoinf*Gdiv+(1|SIN_86),
                   data=coinf,REML=FALSE)
summary(DynAA12.full)

anova(DynAA12.coinS,DynAA12.S)
anova(DynAA12.coinS,DynAA12.coin)
anova(DynAA12.coingdiv,DynAA12.gdiv)
anova(DynAA12.Sgdiv,DynAA12.S)
anova(DynAA12.full)
anova(DynAA12.coingdiv,DynAA12.full)


#####disease abundance in autumn 2010, 2011, 2012 and 2013 ####

#2013#######
coinf<-coinf2013
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
#if we want to limit to populations with at least 3 individuals genotyped
coinf<-coinf[coinf$RA_F2013>0,]
coinf<-coinf[coinf$number_genotyped>2,]
coinf<-coinf[!is.na(coinf$AA_F2013),]
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)

#effect on the absolute abundance between spring and fall with random SIN effect
AA13.S<-lmer(AA_F2013~connec2013+(1|SIN_86),
             data=coinf,REML=FALSE)
summary(AA13.S)
AA13.coin<-lmer(AA_F2013~percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA13.coin)
AA13.gdiv<-lmer(AA_F2013~Gdiv+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA13.gdiv)
AA13.coinS<-lmer(AA_F2013~percoinf+connec2013+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA13.coinS)
AA13.coingdiv<-lmer(AA_F2013~percoinf+Gdiv+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(AA13.coingdiv)
AA13.Sgdiv<-lmer(AA_F2013~connec2013+Gdiv+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA13.Sgdiv)
AA13.full<-lmer(AA_F2013~Gdiv*percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA13.full)

anova(AA13.coin)
anova(AA13.coinS,AA13.S)
anova(AA13.Sgdiv,AA13.S)
anova(AA13.Sgdiv,AA13.gdiv)
anova(AA13.coingdiv,AA13.gdiv)
anova(AA13.coingdiv,AA13.coin)
anova(AA13.coingdiv,AA13.S)
anova(AA13.coingdiv)
anova(AA13.coingdiv,AA13.full)
anova(AA13.full)
anova(AA13.gdiv,AA13.coin,AA13.coingdiv,AA13.full)

plot(table(coinf$AA_F2013))


#2012#######
coinf<-coinf2012
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
#if we want to limit to populations with at least 3 individuals genotyped
coinf<-coinf[coinf$RA_F2012>0,]
coinf<-coinf[coinf$number_genotyped>2,]
coinf<-coinf[!is.na(coinf$AA_F2012),]
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)

#effect on the absolute abundance between spring and fall with random SIN effect
AA12.S<-lmer(AA_F2012~connec2012+(1|SIN_86),
             data=coinf,REML=FALSE)
summary(AA12.S)
AA12.coin<-lmer(AA_F2012~percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA12.coin)
AA12.gdiv<-lmer(AA_F2012~Gdiv+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA12.gdiv)
AA12.coinS<-lmer(AA_F2012~percoinf+connec2012+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA12.coinS)
AA12.coingdiv<-lmer(AA_F2012~percoinf+Gdiv+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(AA12.coingdiv)
AA12.Sgdiv<-lmer(AA_F2012~connec2012+Gdiv+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA12.Sgdiv)
AA12.full<-lmer(AA_F2012~Gdiv*percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA12.full)

anova(AA12.coinS,AA12.S)
anova(AA12.Sgdiv,AA12.S)
anova(AA12.Sgdiv,AA12.gdiv)
anova(AA12.coingdiv,AA12.gdiv)
anova(AA12.coingdiv,AA12.coin)
anova(AA12.coingdiv,AA12.S)
anova(AA12.coingdiv)
anova(AA12.coingdiv,AA12.full)
anova(AA12.full)


#2011#######
coinf<-coinf2011
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
#if we want to limit to populations with at least 3 individuals genotyped
coinf<-coinf[coinf$RA_F2011>0,]
coinf<-coinf[coinf$number_genotyped>2,]
coinf<-coinf[!is.na(coinf$AA_F2011),]
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)

#effect on the absolute abundance between spring and fall with random SIN effect
AA11.S<-lmer(AA_F2011~connec2011+(1|SIN_86),
             data=coinf,REML=FALSE)
summary(AA11.S)
AA11.coin<-lmer(AA_F2011~percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA11.coin)
AA11.gdiv<-lmer(AA_F2011~Gdiv+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA11.gdiv)
AA11.coinS<-lmer(AA_F2011~percoinf+connec2011+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA11.coinS)
AA11.coingdiv<-lmer(AA_F2011~percoinf+Gdiv+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(AA11.coingdiv)
AA11.Sgdiv<-lmer(AA_F2011~connec2011+Gdiv+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA11.Sgdiv)
AA11.full<-lmer(AA_F2011~Gdiv*percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA11.full)

anova(AA11.coinS,AA11.S)
anova(AA11.Sgdiv,AA11.S)
anova(AA11.Sgdiv,AA11.gdiv)
anova(AA11.coingdiv,AA11.gdiv)
anova(AA11.coingdiv,AA11.coin)
anova(AA11.coingdiv,AA11.S)
anova(AA11.coingdiv)
anova(AA11.coingdiv,AA11.full)
anova(AA11.full)


#2010#######
coinf<-coinf2010
#patches with coinfected sampled should have at list 2 different MLG, we correct for that
coinf[coinf$number_coinf>0 & coinf$number_MLG<2,]$number_MLG<-2
#then we create a new variable of genotypic diversity index
coinf<-data.frame(coinf,"Gdiv"=coinf$number_MLG/coinf$number_genotyped)
#creation of a new binary variable with or without coinfection
coinf<-data.frame(coinf,"coinfYN"=coinf$number_coinf)
coinf$coinfYN[(coinf$coinfYN)>0]<-1
coinf$coinfYN<-as.factor(coinf$coinfYN)
#if we want to limit to populations with at least 3 individuals genotyped
coinf<-coinf[coinf$RA_F2010>0,]
coinf<-coinf[coinf$number_genotyped>2,]
coinf<-coinf[!is.na(coinf$AA_F2010),]
coinf<-cbind(coinf,"percoinf"=coinf$number_coinf*100/coinf$number_genotyped)

#effect on the absolute abundance between spring and fall with random SIN effect
AA10.S<-lmer(AA_F2010~connec2010+(1|SIN_86),
             data=coinf,REML=FALSE)
summary(AA10.S)
AA10.coin<-lmer(AA_F2010~percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA10.coin)
AA10.gdiv<-lmer(AA_F2010~Gdiv+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA10.gdiv)
AA10.coinS<-lmer(AA_F2010~percoinf+connec2010+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA10.coinS)
AA10.coingdiv<-lmer(AA_F2010~percoinf+Gdiv+(1|SIN_86),
                    data=coinf,REML=FALSE)
summary(AA10.coingdiv)
AA10.Sgdiv<-lmer(AA_F2010~connec2010+Gdiv+(1|SIN_86),
                 data=coinf,REML=FALSE)
summary(AA10.Sgdiv)
AA10.full<-lmer(AA_F2010~Gdiv*percoinf+(1|SIN_86),
                data=coinf,REML=FALSE)
summary(AA10.full)

anova(AA10.coinS,AA10.S)
anova(AA10.Sgdiv,AA10.S)
anova(AA10.Sgdiv,AA10.gdiv)
anova(AA10.coingdiv,AA10.gdiv)
anova(AA10.coingdiv,AA10.coin)
anova(AA10.coingdiv,AA10.S)
anova(AA10.coingdiv)
anova(AA10.coingdiv,AA10.full)
anova(AA10.full)




# #how to update a model
# model<-update(model, ~ . -)


