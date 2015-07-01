# SelfOut
The R code for the article on recombination of Podosphaera plantaginis in the Aland Islands metapopulation. 

##List of the different scripts

  * **selfout_importdata.R:** script that is used to import and format the raw data. This script is for "internal" use only. It creates the input datafile that are used in the 'selfout_loadata.R' script
  * **selfout_loadata.R:** script you have to run before any other scripts (except for 'selfout_importdata.R')
  * **selfout_colonization.R:** script for analyzing colonisation of patches during the epidemic (*i.e.* intra epidemic season)
  * **selfout_LD_subsample.R:** script to split the data set in different subpopulations so that LD tests can be done taking into account the possible geographical structuration of the populations

