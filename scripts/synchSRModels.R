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


#_________________________________________________________________________
### Compare synchrony in productivity through time w/ residuals
require(synchrony); require(zoo)


ts <- recDat %>% group_by(stk) %>% summarise(tsLength=length(!is.na(ets)), firstYr=min(yr))
selStks <- ts[which(ts$tsLength == "61"), ]$stk
recDatTrim <- recDat[recDat$stk %in% selStks, ]

wideProd <- spread(recDatTrim[, c("yr", "stk", "prod")], stk, prod)
wideResid <- spread(recDatTrim[, c("yr", "stk", "modResid")], stk, modResid)
prodMat <- as.matrix(wideProd[, -1])
residMat <- as.matrix(wideResid[, -1])
yrs <- unique(wideProd$yr)

prodWtdCV <- rollapplyr(prodMat, width = 10, function(x) wtdCV(x), fill = NA,
						by.column = FALSE)
prodSynch <- rollapplyr(prodMat, width = 10, function(x) community.sync(x)$obs, 
						fill = NA, by.column = FALSE)
prodAgCV <- rollapplyr(prodMat, width = 10, function(x) cvAgg(x), fill = NA, 
						by.column = FALSE)
prodCorr <- rollapplyr(prodMat, width = 10, function(x) meancorr(x)$obs, fill = NA, 
						by.column = FALSE)
par(mfrow=c(2,2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
plot(prodWtdCV ~ yrs)
plot(prodSynch ~ yrs)
plot(prodAgCV ~ yrs)
plot(prodCorr ~ yrs)
mtext(side=3, "Fraser Sockeye Metapopulation Dynamics", outer=TRUE)


residWtdCV <- rollapplyr(residMat, width = 10, function(x) wtdCV(x), fill = NA,
						by.column = FALSE)
residSynch <- rollapplyr(residMat, width = 10, function(x) community.sync(x)$obs, 
						fill = NA, by.column = FALSE)
residAgCV <- rollapplyr(residMat, width = 10, function(x) cvAgg(x), fill = NA, 
						by.column = FALSE)
residCorr <- rollapplyr(residMat, width = 10, function(x) meancorr(x)$obs, fill = NA, 
						by.column = FALSE)
par(mfrow=c(2,2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
plot(residWtdCV ~ yrs)
plot(residSynch ~ yrs)
plot(residAgCV ~ yrs)
plot(residCorr ~ yrs)
mtext(side=3, "Fraser Sockeye Metapopulation Dynamics", outer=TRUE)
