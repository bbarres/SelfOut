###############################################################################
###############################################################################
#Exploring the number of Multilocus Genotype (MLG) in space using INLA models
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
#Factors affecting the number of MLG: 2012
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
              "PA_2011","AA_F2012")]<-infections[,c("PLM2_Sept2012",
                                                    "number_coinf",
                                                    "cumulative_sum",
                                                    "connec2012","road_PA",
                                                    "Distance_to_shore",
                                                    "PA_2011","AA_F2012")]
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
model$setLikelihood("nbinomial")

#Intercept-only model
model$setSmoothingModel()
model$addObservationStack(sp=coords, response=infections$number_MLG)
model$estimate()
model$summary() # WAIC = 1486.19

#Full model
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012 + Distance_to_shore + 
                           PA_2011 + cumulative_sum + connec2012 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1460.38

#Full model minus PA_2011
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012 + Distance_to_shore + 
                           cumulative_sum + connec2012 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1458.84

#Full model minus Distance_to_shore
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012 + 
                           PA_2011 + cumulative_sum + connec2012 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1458.90

#Full model minus cumulative_sum
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012 + Distance_to_shore + 
                           PA_2011 + connec2012 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1459.18

#Full model minus road_PA
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012 + Distance_to_shore + 
                           PA_2011 + cumulative_sum + connec2012, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1458.88

#Model with Absolute Abundance in the fall as fixed effect
model$setCovariatesModel(~ 1 + AA_F2012, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1457.46

#Model with Pathogen connectivity as fixed effect
model$setCovariatesModel(~ 1 + connec2012, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1485.47

#Model with Plantago coverage as fixed effect
model$setCovariatesModel(~ 1 + PLM2_Sept2012, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1480.94

#Model with Absolute Abundance of mildew, Plantago coverage and connectivity 
#as fixed effects. This is the best fitting model with nbinomial likelihood
#(lowest WAIC)
model$setCovariatesModel(~ 1 + AA_F2012 + PLM2_Sept2012 + connec2012, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1454.83


# There is no overdispersion, so poisson likelihood is enough
model$setLikelihood("poisson")
model$setCovariatesModel(~ 1 + AA_F2012 + PLM2_Sept2012 + connec2012, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1423.84

fittedResponse <- model$getFittedResponse()
table(infections$number_MLG)
table(trunc(fittedResponse$responseMean))

# Non-spatial models
nonSpatmod<-inla(number_MLG ~ 1, data=infections, family="poisson", 
                 control.predictor=list(compute=TRUE), 
                 control.compute=list(waic=TRUE))
summary(nonSpatmod) # WAIC = 1461.06
nonSpatmodBest<-inla(number_MLG ~ 1 + AA_F2012 + PLM2_Sept2012 + connec2012, 
                     data=infections, family="poisson", 
                     control.predictor=list(compute=TRUE), 
                     control.compute=list(waic=TRUE))
summary(nonSpatmodBest) # WAIC = 1423.98

# The spatial model appears not to improve the fit much
table(round(nonSpatmodBest$summary.fitted.values$mean))


###############################################################################
#Factors affecting the number of MLG: 2013
###############################################################################

infections <- coinf2013
#patches with coinfected sampled should have at list 2 different MLG, 
#we therefore apply a correction to the numnber of MLG for patches 
#with coinfection but with less than 2 single identified MLG
infections[infections$number_coinf>0 & infections$number_MLG<2,]$number_MLG<-2

# Remove missing data and scale covariates
complete <- complete.cases(infections[,c("PLM2_Sept2013","number_coinf",
                                         "cumulative_sum","connec2013",
                                         "road_PA","Distance_to_shore",
                                         "PA_2012","number_MLG","AA_F2013")])
summary(complete)
infections <- infections[complete,]

infections$PA_2012 <- as.numeric(infections$PA_2012)
infections[,c("PLM2_Sept2013","number_coinf",
              "cumulative_sum","connec2013",
              "road_PA","Distance_to_shore",
              "PA_2012","AA_F2013")]<-infections[,c("PLM2_Sept2013",
                                                    "number_coinf",
                                                    "cumulative_sum",
                                                    "connec2013","road_PA",
                                                    "Distance_to_shore",
                                                    "PA_2012","AA_F2013")]
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
model$setLikelihood("nbinomial")

#Intercept-only model
model$setSmoothingModel()
model$addObservationStack(sp=coords, response=infections$number_MLG)
model$estimate()
model$summary() # WAIC = 1756.63

#Full model
model$setCovariatesModel(~ 1 + PLM2_Sept2013 + AA_F2013 + Distance_to_shore + 
                           PA_2012 + cumulative_sum + connec2013 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1720.98

#Full model minus PA_2012
model$setCovariatesModel(~ 1 + PLM2_Sept2013 + AA_F2013 + Distance_to_shore + 
                           cumulative_sum + connec2013 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1720.56

#Full model minus Distance_to_shore
model$setCovariatesModel(~ 1 + PLM2_Sept2013 + AA_F2013 + 
                           PA_2012 + cumulative_sum + connec2013 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1719.90

#Full model minus cumulative_sum
model$setCovariatesModel(~ 1 + PLM2_Sept2013 + AA_F2013 + Distance_to_shore + 
                           PA_2012 + connec2013 + road_PA, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1719.64

#Full model minus road_PA
model$setCovariatesModel(~ 1 + PLM2_Sept2013 + AA_F2013 + Distance_to_shore + 
                           PA_2012 + cumulative_sum + connec2013, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1719.87

#Model with Absolute Abundance in the fall as fixed effect
model$setCovariatesModel(~ 1 + AA_F2013, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1717.08

#Model with Pathogen connectivity as fixed effect
model$setCovariatesModel(~ 1 + connec2013, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1757.32 BAD

#Model with Plantago coverage as fixed effect
model$setCovariatesModel(~ 1 + PLM2_Sept2013, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1757.54 BAD

#Model with Absolute Abundance of mildew and connectivity as fixed effects. 
#This is the best fitting model with nbinomial likelihood (lowest WAIC)
model$setCovariatesModel(~ 1 + AA_F2013 + connec2013, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1716.91


# There is no overdispersion, so poisson likelihood is enough
model$setLikelihood("poisson")
model$setCovariatesModel(~ 1 + AA_F2013 + connec2013, 
                         covariates=infections)
model$clearStack()$addObservationStack(sp=coords,
                                       response=infections$number_MLG,
                                       covariates=infections)
model$estimate()
model$summary() # WAIC = 1679.34

fittedResponse <- model$getFittedResponse()
table(infections$number_MLG)
table(round(fittedResponse$responseMean))

# Non-spatial models
nonSpatmod<-inla(number_MLG ~ 1, data=infections, family="poisson", 
                 control.predictor=list(compute=TRUE), 
                 control.compute=list(waic=TRUE))
summary(nonSpatmod) # WAIC = 1720.61
nonSpatmodBest<-inla(number_MLG ~ 1 + AA_F2013 + connec2013, 
                     data=infections, family="poisson", 
                     control.predictor=list(compute=TRUE), 
                     control.compute=list(waic=TRUE))
summary(nonSpatmodBest) # WAIC = 1678.58

# The spatial model appears not to improve the fit much
table(round(nonSpatmodBest$summary.fitted.values$mean))

