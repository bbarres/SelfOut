###############################################################################
###############################################################################
#Mono-MLG overwintering patches: is the same MLG found the following year?
###############################################################################
###############################################################################

#before using this code, you have to run 'selfout_loadata.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")


###############################################################################
#List the patches with one MLG that overwintered successfully
###############################################################################

#for the overwintering patches of 2011
mon2011<-coinf2011[coinf2011$number_coinf==0 & !is.na(coinf2011$number_coinf) & 
                   coinf2011$PA_S2012==1 & !is.na(coinf2011$PA_S2012) & 
                   coinf2011$number_MLG==1 & !is.na(coinf2011$number_MLG),
                   "patche_ID"]
mon2011<-as.character(mon2011)

#for the overwintering patches of 2012
mon2012<-coinf2012[coinf2012$number_coinf==0 & !is.na(coinf2012$number_coinf) & 
                    coinf2012$PA_S2013==1 & !is.na(coinf2012$PA_S2013) & 
                    coinf2012$number_MLG==1 & !is.na(coinf2012$number_MLG),
                   "patche_ID"]
mon2012<-as.character(mon2012)


###############################################################################
#List the patches with more than one MLG that overwintered successfully
###############################################################################

#for the overwintering patches of 2011
sev2011<-coinf2011[!is.na(coinf2011$number_coinf) & coinf2011$PA_S2012==1 & 
                     !is.na(coinf2011$PA_S2012) & coinf2011$number_MLG>1 & 
                     !is.na(coinf2011$number_MLG),
                   "patche_ID"]
sev2011<-as.character(sev2011)

#for the overwintering patches of 2012
sev2012<-coinf2012[!is.na(coinf2012$number_coinf) & coinf2012$PA_S2013==1 & 
                   !is.na(coinf2012$PA_S2013) & coinf2012$number_MLG>1 & 
                   !is.na(coinf2012$number_MLG),
                   "patche_ID"]
sev2012<-as.character(sev2012)


###############################################################################
#Function that compares the list of MLG from year to the next in the same patch
###############################################################################

#the following function compares a list of MLG sampled in one patch before and 
#after the winter. It returns a table with the patches ID, the number of 
#MLG found in both year and the number of MLG found before winter but not 
#after winter. 
#the arguments of the function are: the list of patches you want to compare, 
#the "year ID" of the year before winter ('year1') and after winter ('year2')
survivingMLG<-function(listpach,year1,year2){
  
  df<-data.frame("patch_ID"=character(),
                 "Survived"=logical(),
                 stringsAsFactors=FALSE)
  
  for (i in 1:length(listpach)){
    temp<-levels(as.data.frame(
             table(geno_hom[geno_hom$survey==year1 & !is.na(geno_hom$survey) & 
                   geno_hom$nb_missing==0 & geno_hom$nb_snp_het==0 &
                   geno_hom$patche_ID==listpach[i] & !is.na(geno_hom$patche_ID),
                   "MLG"]))[,1])  %in% 
          levels(as.data.frame(
             table(geno_hom[geno_hom$survey==year2 & !is.na(geno_hom$survey) & 
                   geno_hom$nb_missing==0 & geno_hom$nb_snp_het==0 &
                   geno_hom$patche_ID==listpach[i] & !is.na(geno_hom$patche_ID),
                   "MLG"]))[,1])
    temp<-cbind("patch_ID"=listpach[i],"persist"=temp)
    df<-rbind(df,temp)
  }
  
  return(table(df))
  
}


levels(as.data.frame(table(geno_hom[geno_hom$survey=="MS12" & !is.na(geno_hom$survey) & 
                           geno_hom$nb_missing==0 & geno_hom$nb_snp_het==0 &
                           geno_hom$patche_ID==mon2012[5] & 
                           !is.na(geno_hom$patche_ID),
                             "MLG"]))[,1]) %in% 
levels(as.data.frame(table(geno_hom[geno_hom$survey=="MS13" & !is.na(geno_hom$survey) & 
                           geno_hom$nb_missing==0 & geno_hom$nb_snp_het==0 &
                           geno_hom$patche_ID==mon2012[5] & 
                           !is.na(geno_hom$patche_ID),
                             "MLG"]))[,1])

levels(as.data.frame(table(geno_hom[geno_hom$survey=="MS12" & !is.na(geno_hom$survey) & 
                                    geno_hom$nb_missing==0 & geno_hom$nb_snp_het==0 &
                                    geno_hom$patche_ID==sev2012[5] & 
                                    !is.na(geno_hom$patche_ID),
                                    "MLG"]))[,1]) %in% 
levels(as.data.frame(table(geno_hom[geno_hom$survey=="MS13" & !is.na(geno_hom$survey) & 
                                    geno_hom$nb_missing==0 & geno_hom$nb_snp_het==0 &
                                    geno_hom$patche_ID==sev2012[5] & 
                                    !is.na(geno_hom$patche_ID),
                                    "MLG"]))[,1])

###############################################################################
#Runing the function
###############################################################################

temp<-survivingMLG(sev2012,"MS12","MS13")
colMeans(temp)
temp<-survivingMLG(mon2012,"MS12","MS13")
colMeans(temp)




