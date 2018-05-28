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

require(here); require(dplyr); require(tidyr); require(lme4); require(MuMIn)

source(here("scripts/synchFunctions.R"))

recDat <- read.csv(here("/data/scaledFullData.csv")) #cleaned environmental and catch covariates
recDat$stk <- as.factor(recDat$stk)

#_________________________________________________________________________
### Run model selection comparing performance of subset of covariates

## Check base ricker model behaves
ric1 <- lmer(prod ~ ets + (1 + ets | stk), data = recDat)
summary(ric1)
coef(ric1) 

allVars <- names(recDat[6:16])
catchVars <- names(recDat[6:7])
envVars <- names(recDat[8:16])
modScoreFull <- NULL
modNameFull <- NULL

for (j in seq_along(allVars)) { #one covariate models
	mod <- lmer(prod ~ ets + recDat[, allVars[j]] + (1 + ets | stk), REML = FALSE, 
				data = recDat)
	modName <- paste(allVars[j])
	modScore <- AICc(mod)[1]
	modNameFull <- c(modNameFull, modName)
	modScoreFull <- c(modScoreFull, modScore)
}

for (i in seq_along(catchVars)) { #two covariate models
	for (h in seq_along(envVars)) {
		mod <- lmer(prod ~ ets + recDat[, catchVars[i]] + recDat[, envVars[h]] + 
					(1 + ets | stk), REML = FALSE, data = recDat) #assume covariate effects are shared
		modName <- paste(catchVars[i], envVars[h], sep = "_")
		modScore <- AICc(mod)[1]
		modNameFull <- c(modNameFull, modName)
		modScoreFull <- c(modScoreFull, modScore)
	}
}

aicScores <- data.frame(mod = modNameFull,
						aicc = modScoreFull)
aicScores <- with(aicScores, aicScores[order(aicc), ])

# Calculate residuals from mixed effects and multiple models
topRicMixed <- lmer(prod ~ ets + pinkCatch + sstPrime + (1 + ets | stk), REML = TRUE, 
			   data = recDat)
mixedResid <- resid(topRicMixed)

stkIDs <- unique(recDat$stk)
residVec <- NULL
for (i in seq_along(stkIDs)) {
	d <- recDat[recDat$stk == stkIDs[i], ]
	singleRic <- lm(prod ~ ets + pinkCatch + sstPrime, data = d)
	residVec <- c(residVec, resid(singleRic))
}

recDat$mixedModResid <- mixedResid
recDat$modResid <- residVec