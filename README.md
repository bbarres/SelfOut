[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.34573.svg)](http://dx.doi.org/10.5281/zenodo.34573)

# SelfOut

The R and the Matlab code for the article on recombination of a plant pathogenic fungus, *Podosphaera plantaginis*, in the Åland Islands metapopulation. 

---
 

       2012                |           2013          
:-------------------------:|:-------------------------:
![alt text](https://q77cda.db.files.1drv.com/y4mLortDGyBMAXBZal2NOntFHOssZmA_dfBG8w1CkEJig4DGM6Y_gEk8nDwCmcPswBT8Siyc7mB40FVB3YGWMglc45u7J9P_XSTadkfZ2n4GB9ILG2ueO3mZb1RMrtdUeBaLmqQsYpFukxMp2m599Pl4lIWGyBse6R2bgAQ1jAH5xEH4c_IqzSaREQw2yp1qZKxAaNU0Z1CFEqUMvEBxF1AsQ "A map of coinfected and single infected samples in 2012")  |  ![alt text](https://q77dda.db.files.1drv.com/y4muIZV4ZQ3-m1NSS82eTtX0f_e0J6zalHgGmgId6RCCLxEfszhu5_Y8u4HtY5qoysWiwiPNKGKLqjNt97jebjDZm9deQ7DxrHj5WkyzJs2I9sRH1G4QPhNvoWd8yE0-Q9ixOmcFi-SGXeagz5mOtQdxANcYPKbwRN4m260kxsVv3xIJBG9zS7bIEZWwhcOht2gxQm6LBuXCwqlmHgAn2I4zA "A map of coinfected and single infected samples in 2013")



>These maps show the samples that were found to be coinfected (red) or single infected (green) in 2012 and 2013. This study is a follow up of the 2015 article of **_Susi et al_**: [Co-infection alters population dynamics of infectious disease](http://www.nature.com/ncomms/2015/150108/ncomms6975/full/ncomms6975.html). 

# 


## List of the different scripts

  * **coinfModel:** This folder contains all the Matlab scripts used for running the ABC analysis of the relative frequency of single infections vs co-infectionson in the host metapopulation.These scripts were developped by Jukka P. Sirén. 
  * **selfout_importdata.R:** script that is used to import and format the raw data. This script is for "internal" use only. It creates the input datafiles that are used in the `selfout_loadata.R` script.
  * **selfout_loadata.R:** script you have to run before any other scripts (except for `selfout_importdata.R`).
  * **selfout_glm.R:** script for the generalized linear model used to analysis the impact of coinfection on different life traits of *P. plantaginis*.
  * **selfout_colonization.R:** script for analyzing colonization of patches during the epidemic (*i.e.* intra epidemic season).
  * **selfout_monoMLG.R:** script to identify patches with a unique MLG that overwintered successfully and to check if the MLG is found in the same patch the next year. The patches with several MLG are also investigated: we check which MLG overwintered successfully and which didn't. 
  * **selfout_LD_subsample.R:** script to split the data set in different subpopulations so that LD tests can be done taking into account the possible geographical structuration of the populations.
  * **selfout_samplingmap.R:** script to plot the sampling maps for each year of the survey (2010, 2011, 2012 and 2013). 
  * **selfout_barplot.R:** script for the barplot figures that show several aspect of the consequence of coinfection on the disease dynamic of the ribwort powdery mildew. 
  * **selfout_newMLGplot.R:** script for the plot figure that show the number of new MLG as a function of the proportion of coinfection in patch the previous year. 
  * **selfout_distintrapatch.R:** script to compute the mean distances between samples within patch, for every patch. 
  * **selfout_diversitymap.R:** script for producing map of the Allelic richness and Genotypic richness in Åland patches. 
  * **selfout_MLG_inla.R:** script to explore the variation of the number of MLG across space. This code was adapted from a [repository](https://github.com/statguy/MultiLocusGenotype) of [Jussi Jousimo](https://github.com/statguy). It uses the `SpaceTimeModels` package that allows an easier use of INLA to fit spatial and time models in an R framework. You can install this package using [this](https://github.com/statguy/SpaceTimeModels) github repository. 
  * **selfout_coinfectionmap.R:** script to produce the map of the relative risk surface of coinfection vs single infection in 2012 and 2013. See below for an example of the output that can be produced. 
  * **selfout_GenoRichmap.R:** script to produce the map of the genotypic richness in 2012 and 2013. The genotypic richness was evaluated on a minimum sample size of 3 successfully genotyped pure individuals per patch. A spatial smoothing of the genotypic richness values was performed with the same bandwidth used for the coinfection map. 


## Some map examples

       2012               |           2013          
:-------------------------:|:-------------------------:
![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/CoinRisk12.png "Relative risk maps of coinfection vs infection in 2012. This map was obtained using 'spatstat' R package")  |  ![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/CoinRisk13.png "Relative risk maps of coinfection vs infection in 2013. This map was obtained using 'spatstat' R package")

>These maps show the relative risk to find a coinfected sample in the metapopulation of the plant pathogen *P. plantaginis*. The darker the red, the higher the probability to find a plant with multiple strains on a leaf. These maps show that coinfection is found everywhere in Åland Islands, that there is a lot of variation of the prevalence of coinfection across the system and that there is also some temporal variation of the distribution of coinfection from one year to another. 


---

![alt text](http://drive.google.com/drive/folders/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/Plantago_plant_model2.png "A drawing of a Plantago lanceolata")


---



