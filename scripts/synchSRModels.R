#*************************************************************************************
# synchSRModels.R
# Date revised: May 26, 2018
# Inputs: stock-recruit and environmental data
# Explainer: 1) Fits basic stock recruit models,
#			 2) Incorporates environmental covariates and IDs top model
# 			 3) Tests whether synchrony patters dissipate in residuals of top model 
#               (i.e. is there evidence that synch is associated with environmental drivers)
#*************************************************************************************

# setwd("C:/github/synchSalmon/")
setwd("/Users/cam/github/synchSalmon") #Cam's Mac wd

require(here); require(dplyr); require(tidyr)

source(here("scripts/synchFunctions.R"))

## Import covariates
sstPca <- read.table(here("/data/env/cleanCovar.csv")) #cleaned environmental and catch covariates

