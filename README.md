#SelfOut

The R code for the article on recombination of a plant pathogenic fungus, *Podosphaera plantaginis*, in the Åland Islands metapopulation. 

---


       2012                |           2013          
:-------------------------:|:-------------------------:
![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/Samp2012.png "A map of coinfected and single infected samples in 2012")  |  ![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/Samp2013.png "A map of coinfected and single infected samples in 2013")


---


>These maps show the samples that were found to be coinfected (red) or single infected (green) in 2012 and 2013. This study is a follow up of the 2015 article of **_Susi et al_**: [Co-infection alters population dynamics of infectious disease](http://www.nature.com/ncomms/2015/150108/ncomms6975/full/ncomms6975.html). 



##List of the different scripts

  * **selfout_importdata.R:** script that is used to import and format the raw data. This script is for "internal" use only. It creates the input datafiles that are used in the `selfout_loadata.R` script.
  * **selfout_loadata.R:** script you have to run before any other scripts (except for `selfout_importdata.R`).
  * **selfout_glm.R:** script for the generalized linear model used to analysis the impact of coinfection on different life traits of *P. plantaginis*.
  * **selfout_colonization.R:** script for analyzing colonization of patches during the epidemic (*i.e.* intra epidemic season).
  * **selfout_monoMLG.R:** script to identify patches with a unique MLG that overwintered successfully and to check if the MLG is found in the same patch the next year. The patches with several MLG are also investigated: we check which MLG overwintered successfully and which didn't. 
  * **selfout_LD_subsample.R:** script to split the data set in different subpopulations so that LD tests can be done taking into account the possible geographical structuration of the populations.
  * **selfout_samplingmap.R:** script to plot the sampling maps for each year of the survey (2010, 2011, 2012 and 2013). 
  * **selfout_barplot.R:** script for the barplot figures that show several aspect of the consequence of coinfection on the disease dynamic of the ribwort powdery mildew. 
  * **selfout_distintrapatch.R:** script to compute the mean distances between samples within patch, for every patch.  
  * **selfout_diversitymap.R:** script for producing map of the Allelic richness and Genotypic richness in Åland patches. 
  * **selfout_MLG_INLA.R:** script to explore the variation of the number of MLG across space. This code was adapted from a [repository](https://github.com/statguy/MultiLocusGenotype) of [Jussi Jousimo](https://github.com/statguy). It uses the `SpaceTimeModels` package that allows an easier use of INLA to fit spatial and time models in an R framework. You can install this package using [this](https://github.com/statguy/SpaceTimeModels) github repository. 


##Some map examples

       2012               |           2013          
:-------------------------:|:-------------------------:
![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/CoinRisk12.png "Relative risk maps of coinfection vs infection in 2012. This map was obtained using 'spatstat' R package")  |  ![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/CoinRisk13.png "Relative risk maps of coinfection vs infection in 2013. This map was obtained using 'spatstat' R package")

>These maps show the relative risk to find a coinfected sample in the metapopulation of the plant pathogen *P. plantaginis*. The darker the red, the higher the probability to find a plant with multiple strains on a leaf. These maps show that coinfection is found everywhere in Åland Islands, that there is a lot of variation of the prevalence of coinfection across the system and that there is also some temporal variation of the distribution of coinfection from one year to another. 


---

![alt text](http://googledrive.com/host/0B-FIusWb7o6PfjdhbUJncm1mdjM1NnQ1TWl6MHhZUnNRZjd6RkUtUVo5WlFsVURTV0lvQjA/Plantago_plant_model2.png "A drawing of a Plantago lanceolata")


---



