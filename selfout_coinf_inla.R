###############################################################################
###############################################################################
#Exploring the number of Coinfections in space using INLA models
###############################################################################
###############################################################################

#this code is adapted from the repository of Jussi Jousimo: 
#https://github.com/statguy/MultiLocusGenotype. It uses the SpaceTimeModels 
#package that allows an easier use of INLA to fit spatial and time models in 
#an R framework

#before using this code, you have to run 'selfout_loadata.R' first
setwd("~/work/Rfichiers/Githuber/SelfOut_data")

#the "SpaceTimeModels" package wasn't available on CRAN at the time this code 
#was written. So you have to install the 'devtools' package and install the 
#package straight from github
library(devtools)
install_github("statguy/SpaceTimeModels")
library(SpaceTimeModels)


###############################################################################
#Factors affecting the number of coinfections: 2012
###############################################################################

infections <- coinf2012
#patches with coinfected sampled should have at list 2 different MLG, 
#we therefore apply a correction to the numnber of MLG for patches 
#with coinfection but with less than 2 single identified MLG
infections[infections$number_coinf>0 & infections$number_MLG<2,]$number_MLG<-2

# Remove missing data and scale covariates
complete <- complete.cases(infections[,c("PLM2_Sept2012","number_coinf",
                                         "cumulative_sum","connec2012",
                                         "road_PA","Distance_to_shore",
                                         "PA_2011","number_MLG","AA_F2012")])
summary(complete)
infections <- infections[complete,]

infections$PA_2011 <- as.numeric(infections$PA_2011)
infections[,c("PLM2_Sept2012","number_coinf",
              "cumulative_sum","connec2012",
              "road_PA","Distance_to_shore",
              "PA_2011","AA_F2012")]<-scale(infections[,c("PLM2_Sept2012",
                                                    "number_coinf",
                                                    "cumulative_sum",
                                                    "connec2012","road_PA",
                                                    "Distance_to_shore",
                                                    "PA_2011","AA_F2012")])
coords <- sp::SpatialPoints(infections[,c("Longitude","Latitude")])
#infections$number_MLG <- infections$number_MLG - 1

#Construct estimation mesh
mesh <- SpaceTimeModels::NonConvexHullMesh$new(knots=coords,knotsScale=1e5, 
                                               cutoff=1e3,maxEdge=c(2.2e3,1e5),
                                               convex=0.1)
mesh$getINLAMesh()$n
mesh$plot()

#Setup spatial model
model<-ContinuousSpaceModel$new()
model$setSpatialMesh(mesh)
model$setSpatialPrior()
model$setLikelihood("binomial")

#Intercept-only model
model$setSmoothingModel()
model$addObservationStack(sp=coords, response=infections$number_coinf,
                          offset=infections$number_genotyped)
model$estimate()
model$summary() # WAIC = 1486.19
