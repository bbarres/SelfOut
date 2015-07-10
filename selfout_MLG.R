###############################################################################
###############################################################################
#Exploring the number of Multilocus Genotype (MLG) using INLA models
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
#Factors affecting the number of MLG
###############################################################################

infections <- coinf2012

# Remove missing data and scale covariates
complete <- complete.cases(infections[,c("PLM2_Sept2012","number_coinf",
                                         "cumulative_sum","connec2012",
                                         "road_PA","Distance_to_shore",
                                         "PA_2011","number_MLG","AA_F2012")])
length(complete)
infections <- infections[complete,]

infections$PA_2011 <- as.numeric(infections$PA_2011)
infections[,c("PLM2_Sept2012","number_coinf",
              "cumulative_sum","connec2012",
              "road_PA","Distance_to_shore",
              "PA_2011","AA_F2012")] <- infections[,c("PLM2_Sept2012","number_coinf",
                                                      "cumulative_sum","connec2012",
                                                      "road_PA","Distance_to_shore",
                                                      "PA_2011","AA_F2012")]
coords <- sp::SpatialPoints(infections[,c("Longitude","Latitude")])
#infections$number_MLG <- infections$number_MLG - 1

# Construct estimation mesh
mesh <- SpaceTimeModels::NonConvexHullMesh$new(knots=coords, knotsScale=1e5, cutoff=1e3, maxEdge=c(2.2e3, 1e5), convex=0.1)
mesh$getINLAMesh()$n
mesh$plot()

# Setup spatial model
model <- ContinuousSpaceModel$new()
model$setSpatialMesh(mesh)
model$setSpatialPrior()
model$setLikelihood("nbinomial")

# Intercept-only model
model$setSmoothingModel()
model$addObservationStack(sp=coords, response=infections$number_MLG)
model$estimate()
model$summary() # WAIC = 1486.21

# Full model
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012 + Distance_to_shore + PA_2011, covariates=infections)
model$clearStack()$addObservationStack(sp=coords, response=infections$number_MLG, covariates=infections)
model$estimate()
model$summary() # WAIC = 1458.55

model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012 + PA_2011, covariates=infections)
model$clearStack()$addObservationStack(sp=coords, response=infections$number_MLG, covariates=infections)
model$estimate()
model$summary() # WAIC = 1457.05

model$setCovariatesModel(~ 1 + AA_F2012 + PA_2011, covariates=infections)
model$clearStack()$addObservationStack(sp=coords, response=infections$number_MLG, covariates=infections)
model$estimate()
model$summary() # WAIC = 1458.94

# Best fitting model with nbinomial likelihood (lowest WAIC)
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012, covariates=infections)
model$clearStack()$addObservationStack(sp=coords, response=infections$number_MLG, covariates=infections)
model$estimate()
model$summary() # WAIC = 1455.63
fittedResponse <- model$getFittedResponse()
table(infections$number_MLG)
table(round(fittedResponse$responseMean))

model$setCovariatesModel(~ 1 + AA_F2012, covariates=infections)
model$clearStack()$addObservationStack(sp=coords, response=infections$number_MLG, covariates=infections)
model$estimate()
model$summary() # WAIC = 1457.46

# There is no overdispersion, so poisson likelihood is enough
model$setLikelihood("poisson")
model$setCovariatesModel(~ 1 + PLM2_Sept2012 + AA_F2012, covariates=infections)
model$clearStack()$addObservationStack(sp=coords, response=infections$number_MLG, covariates=infections)
model$estimate()
model$summary() # WAIC = 1424.59

fittedResponse <- model$getFittedResponse()
table(infections$number_MLG)
table(round(fittedResponse$responseMean))

# Non-spatial models
summary(inla(number_MLG ~ 1, data=infections, family="poisson", control.predictor=list(compute=TRUE), control.compute=list(waic=TRUE)))
# WAIC = 1461.06
fit <- inla(number_MLG ~ 1 + PLM2_Sept2012 + AA_F2012, data=infections, family="poisson", control.predictor=list(compute=TRUE), control.compute=list(waic=TRUE))
summary(fit)
# WAIC = 1425.23
# The spatial model appears not to improve the fit much
table(round(fit$summary.fitted.values$mean))