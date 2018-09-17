#*************************************************************************************
# synchExplore.R
# Date revised: April 20, 2018; ONGOING
# Inputs: stock-recruit data csv
# Outputs: figure pdfs
# Explainer: Uses models in Thibaut and Connolly 2013 plus code from S. Anderson to 
#           first explore how component variability and synchrony contribute to aggregate 
#           variability in Fraser sockeye, then identify covariates and track changes 
#			in residuals
#*************************************************************************************

setwd("C:/github/synchSalmon/")
# setwd("/Users/cam/github/synchSalmon") #Cam's Mac wd

require(here); require(synchrony); require(zoo); require(ggplot2); require(dplyr); require(tidyr)

source(here("scripts/synchFunctions.R"))

## Data clean
recDat1 <- read.csv(here("/data/sox/fraserRecDatEFF.csv"), stringsAsFactors = FALSE)
recDat2 <- read.csv(here("/data/sox/fraserRecDatTrim.csv"), stringsAsFactors = FALSE) 
recDat <- merge(recDat1, recDat2[, c("stk", "yr", "ets")], by = c("stk", "yr")) #combine ets and eff estimates
recDat <- with(recDat, recDat[order(stk, yr),])
recDat$prod <- log(recDat$rec/recDat$ets)
for (i in 1:nrow(recDat)) { #add top-fitting SR model
  stk <- recDat$stk[i]
  if (stk == 1 | stk == 2 | stk == 6 | stk == 8 | stk == 9) {
    recDat$model[i] <- "larkin"
  }	else {
    recDat$model[i] <- "ricker"
  } 
}
# Add lagged ETS abundances for larkin models
recDat$ets1 <- NA
recDat$ets2 <- NA
recDat$ets3 <- NA
stkIndex <- unique(recDat$stk)
for(j in seq_along(stkIndex)) {
  d <- subset(recDat, stk == stkIndex[j])
  for (i in 1:nrow(d)) { #add top-fitting SR model
    d$ets1[i] <- ifelse(i < 1, NA, d$ets[i - 1])
    d$ets2[i] <- ifelse(i < 2, NA, d$ets[i - 2])
    d$ets3[i] <- ifelse(i < 3, NA, d$ets[i - 3])
  }
  recDat[recDat$stk == stkIndex[j], c("ets1", "ets2", "ets3")] <- d[, c("ets1", "ets2", "ets3")]
}
# Trim and convert to matrix
ts <- recDat %>% 
  group_by(stk) %>% 
  summarise(tsLength = length(!is.na(prod)), firstYr = min(yr), lastYr = max(yr))
selectedStks <- c(1, seq(from=3, to=10, by=1), 18, 19) #stocks w/ time series from 1948
selectedStks2 <- seq(from = 1, to = 19, by =1 )[-c(11, 15, 17)] #stocks w/ time series from 1973
recDatTrim1 <- subset(recDat, !is.na(prod) & !is.na(eff3) & !yr == "2012")
recDatTrim <- recDatTrim1[recDatTrim1$stk %in% selectedStks,]
recDatTrim2 <- recDatTrim1[recDatTrim1$stk %in% selectedStks2,]
recDatTrim2 <- recDatTrim2[recDatTrim2$yr > 1972, ]
# recDatTrim2 <- subset(recDat, !is.na(prod) & !is.na(eff3) & !yr == "2012")
# write.csv(recDatTrim, file = here("data/recDatLongTS.csv"), row.names = FALSE)

wideRec <- spread(recDatTrim[,c("stk", "yr", "ets")], stk, ets)
recMat <- as.matrix(wideRec[,-1])
wideProd <- spread(recDatTrim[,c("stk", "yr", "prod")], stk, prod)
prodMat <- as.matrix(wideProd[,-1])
yrs <- unique(wideProd$yr)
wideProd2 <- spread(recDatTrim2[,c("stk", "yr", "prod")], stk, prod)
prodMat2 <- as.matrix(wideProd2[,-1])


#_________________________________________________________________________
## Frequentist framework
rollWtdCVS <- rollapplyr(recMat, width=10, function(x) wtdCV(x), fill=NA, by.column=FALSE)
rollSynchS <- rollapplyr(recMat, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
rollAgCVS <- rollapplyr(recMat, width=10, function(x) cvAgg(x), fill=NA, by.column=FALSE)
rollCorrS <- rollapplyr(recMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
rollWtdCV <- rollapplyr(prodMat, width=10, function(x) wtdCV(x, recMat = recMat), fill=NA, by.column=FALSE)
rollSynch <- rollapplyr(prodMat, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
rollAgCV <- rollapplyr(prodMat, width=10, function(x) cvAgg(x, recMat = recMat), fill=NA, by.column=FALSE)
rollCorr <- rollapplyr(prodMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
rollS <- rollapplyr(recMat, width=10, function(x) wtdMean(x), fill=NA, by.column=FALSE)
rollProd <- rollapplyr(prodMat, width=10, function(x) wtdMean(x, recMat = recMat), fill=NA, by.column=FALSE)


## Fit Ricker and Larkin models to each CU and examine trends in residuals (accounts for 
## density dependence and changes in exploitation rates through time)
stkIndex <- unique(recDatTrim$stk)
residVec <- NULL
# par(mfrow = c(3, 3))
for(j in seq_along(stkIndex)) {
  d <- subset(recDatTrim, stk == stkIndex[j])
  mod <- unique(d$model)
  if (mod == "ricker") {
    srMod <- lm(prod ~ eff, data = d)
  }
  if (mod == "larkin") {
    srMod <- lm(prod ~ eff + eff1 + eff2 + eff3, data = d)
  }
  residVec <- c(residVec, resid(srMod))
  # acf(resid(srMod), main = stkIndex[j])
}
recDatTrim$logResid <- residVec
recDatTrim$modResid <- exp(residVec)

wideResid <- spread(recDatTrim[ , c("stk", "yr", "modResid")], stk, modResid)
residMat <- as.matrix(wideResid[ , -1])
wideLogResid <- spread(recDatTrim[ , c("stk", "yr", "logResid")], stk, logResid)
logResidMat <- as.matrix(wideLogResid[ , -1])


rollWtdCVR <- rollapplyr(residMat, width=10, function(x) wtdCV(x, recMat = recMat), fill = NA, by.column = FALSE)
rollSynchR <- rollapplyr(residMat, width=10, function(x) community.sync(x)$obs, fill = NA, by.column = FALSE)
rollAgCVR <- rollapplyr(residMat, width=10, function(x) cvAgg(x, recMat = recMat), fill = NA, by.column = FALSE)
rollLogResid <- rollapplyr(logResidMat, width=10, function(x) wtdMean(x, recMat = recMat), fill=NA, by.column=FALSE)

## Plot trends
pdf(here("outputs/expFigs/varTrends.pdf"), height = 6, width = 6)
par(mfrow=c(2, 2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
plot(rollWtdCVS ~ yrs, type = "l", ylab = "Weighted Mean Component CV")
plot(rollSynchS ~ yrs, type = "l", ylab = "Synchrony Index")
plot(rollAgCVS ~ yrs, type = "l", ylab = "Aggregate CV")
plot(rollS2 ~ yrs, type = "l", ylab = "Mean")
mtext(side=3, "Variability in Recruit Abundance", outer=TRUE)
plot(rollWtdCV ~ yrs, type = "l", ylab = "Weighted Mean Component CV")
plot(rollSynch ~ yrs, type = "l", ylab = "Synchrony Index")
plot(rollAgCV ~ yrs, type = "l", ylab = "Aggregate CV")
plot(rollProd ~ yrs, type = "l", ylab = "Mean")
mtext(side=3, "Variability in log(R/S)", outer=TRUE)
plot(rollWtdCVR ~ yrs, type = "l", ylab = "Weighted Mean Component CV")
plot(rollSynchR ~ yrs, type = "l", ylab = "Synchrony Index")
plot(rollAgCVR ~ yrs, type = "l", ylab = "Aggregate CV")
plot(rollLogResid ~ yrs, type = "l", ylab = "Mean")
mtext(side=3, "Variability in Model Residuals", outer=TRUE)
dev.off()

rollCorrR <- rollapplyr(residMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
rollCorrLogR <- rollapplyr(logResidMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
plot(rollCorrR ~ yrs, type = "l", ylab = "Mean Correlation Residuals")
plot(rollCorrLogR ~ yrs, type = "l", ylab = "Mean Correlation Log Residuals")


#_________________________________________________________________________
## Look at metrics from Yamane et al. 2018
rollCurrentCV <- rollapplyr(recMat, width = 10, function(x) calcCV(x, current = TRUE), fill = NA, by.column = FALSE)
rollNullCV <- rollapplyr(recMat, width = 10, function(x) calcCV(x, current = FALSE), fill = NA, by.column = FALSE)
rollDivDefS <- rollCurrentCV - rollNullCV

plot(rollDivDefS ~ yrs, type = "l", ylab = "Diversity Deficit (null)")  
plot(rollDivDefS ~ rollSynchS)

rollCurrentCV <- rollapplyr(prodMat, width = 10, function(x) calcCV(x, current = TRUE), fill = NA, by.column = FALSE)
rollNullCV <- rollapplyr(prodMat, width = 10, function(x) calcCV(x, current = FALSE), fill = NA, by.column = FALSE)
rollDivDefProd <- rollCurrentCV - rollNullCV

plot(rollDivDefProd ~ yrs, type = "l", ylab = "Diversity Deficit (null)")  
plot(rollDivDefProd ~ rollSynch)
abline(a = 0, b = 1)
cor(rollDivDefProd[10:39], rollSynch[10:39])


#_________________________________________________________________________
## Look at CU-specific time series of residuals and productivity
pdf(here("outputs/expFigs/prodByCU.pdf"), height = 6, width = 6)
par(mfrow=c(4, 3), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
for(i in seq_along(stkIndex)) {
  d <- prodMat[ , i]
  plot(d ~ yrs, type = "l", ylab = "log R/S")
}
par(mfrow=c(4, 3), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
for(i in seq_along(stkIndex)) {
  d <- logResidMat[ , i]
  plot(d ~ yrs, type = "l", ylab = "Residuals")
  abline(h = 0)
}
par(mfrow=c(4, 3), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
for(i in seq_along(stkIndex)) {
  d <- logResidMat[ , i]
  hist(d, ylab = "Residuals")
  abline(v = 0, col = "red")
}
dev.off()

