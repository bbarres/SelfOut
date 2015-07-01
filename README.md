# SelfOut


The R code for the article on recombination of *Podosphaera plantaginis* in the Åland Islands metapopulation. 

---

![alt text](http://www.nature.com/ncomms/2015/150108/ncomms6975/images/ncomms6975-f3.jpg "A map of coinfection in Aland in 2012")


---


>This map is extracted from the 2015 article of **_Susi et al_**: [Co-infection alters population dynamics of infectious disease](http://www.nature.com/ncomms/2015/150108/ncomms6975/full/ncomms6975.html) published in [*Nature Communications*](http://www.nature.com/ncomms/index.html). *(This will be updated as soon as a newer figure is available)* 



##List of the different scripts

  * **selfout_importdata.R:** script that is used to import and format the raw data. This script is for "internal" use only. It creates the input datafile that are used in the 'selfout_loadata.R' script
  * **selfout_loadata.R:** script you have to run before any other scripts (except for 'selfout_importdata.R')
  * **selfout_colonization.R:** script for analyzing colonization of patches during the epidemic (*i.e.* intra epidemic season)
  * **selfout_LD_subsample.R:** script to split the data set in different subpopulations so that LD tests can be done taking into account the possible geographical structuration of the populations
  * **selfout_samplingmap.R:** script to plot the sampling maps for each year of the survey (2010, 2011, 2012 and 2013)
  * **selfout_glm.R:** script for the generalized linear model used to analysis the impact of coinfection on different life traits of *P. plantaginis*

