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
recDat$prod <- log(recDat$rec/recDat$eff)
for (i in 1:nrow(recDat)) { #add top-fitting SR model
  stk <- recDat$stk[i]
  if (stk == 1 | stk == 2 | stk == 6 | stk == 8 | stk == 9) {
    recDat$model[i] <- "larkin"
  }	else {
    recDat$model[i] <- "ricker"
  } 
}
# Add lagged eff abundances for larkin models
recDat$eff1 <- NA
recDat$eff2 <- NA
recDat$eff3 <- NA
stkIndex <- unique(recDat$stk)
for(j in seq_along(stkIndex)) {
  d <- subset(recDat, stk == stkIndex[j])
  for (i in 1:nrow(d)) { #add top-fitting SR model
    d$eff1[i] <- ifelse(i < 1, NA, d$eff[i - 1])
    d$eff2[i] <- ifelse(i < 2, NA, d$eff[i - 2])
    d$eff3[i] <- ifelse(i < 3, NA, d$eff[i - 3])
  }
  recDat[recDat$stk == stkIndex[j], c("eff1", "eff2", "eff3")] <- d[, c("eff1", "eff2", "eff3")]
}
# Trim and convert to matrix
ts <- recDatTrim %>% 
  group_by(stk) %>% 
  summarise(tsLength = length(!is.na(prod)), firstYr = min(yr), lastYr = max(yr))
selectedStks <- c(1, seq(from=3, to=10, by=1), 18, 19)
recDatTrim <- recDat[recDat$stk %in% selectedStks,]
recDatTrim <- subset(recDatTrim, !is.na(prod) & !is.na(eff3) & !yr == "2012")


wideRec <- spread(recDatTrim[,c("stk", "yr", "ets")], stk, ets)
recMat <- as.matrix(wideRec[,-1])


#_________________________________________________________________________
## Frequentist framework
rollWtdCV <- rollapplyr(recMat, width=10, function(x) wtdCV(x), fill=NA, by.column=FALSE)
rollSynch <- rollapplyr(recMat, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
rollAgCV <- rollapplyr(recMat, width=10, function(x) cvAgg(x), fill=NA, by.column=FALSE)
rollCorr <- rollapplyr(recMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
yrs <- unique(wideRec$yr)

## Raw trends
pdf(here("outputs/expFigs/spwnrVar.pdf"), height = 5, width = 6)
par(mfrow=c(2,2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
plot(rollWtdCV ~ yrs, type = "l")
plot(rollSynch ~ yrs, type = "l")
plot(rollAgCV ~ yrs, type = "l")
plot(rollCorr ~ yrs, type = "l")
mtext(side=3, "Variability in Spawner Abundance", outer=TRUE)
dev.off()

cor(y = rollAgCV[10:61], x = rollWtdCV[10:61])
cor(y = rollAgCV[10:61], x = rollSynch[10:61])

#_________________________________________________________________________
## Fit Ricker and Larkin models to each and examine trends in residuals (accounts for 
## density dependence and changes in exploitation rates through time)
# stkIndex <- unique(recDatTrim$stk)
# residVec <- NULL
# for(j in seq_along(stkIndex)) {
#   d <- subset(recDatTrim, stk == stkIndex[j])
#   mod <- unique(d$model)
#   if (mod == "ricker") {
#     srMod <- lm(prod ~ eff, data = d)
#   }
#   if (mod == "larkin") {
#     srMod <- lm(prod ~ eff + eff1 + eff2 + eff3, data = d)
#   }
#   residVec <- c(residVec, resid(srMod))
# }
# recDatTrim$modResid <- residVec

# NOTE: results in nonsensical variance estimates because residuals are centered on 0, focus
# on log(R/S) instead

wideProd <- spread(recDatTrim[,c("stk", "yr", "prod")], stk, prod)
prodMat <- as.matrix(wideProd[,-1])

rollWtdCV <- rollapplyr(prodMat, width=10, function(x) wtdCV(x), fill=NA, by.column=FALSE)
rollSynch <- rollapplyr(prodMat, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
rollAgCV <- rollapplyr(prodMat, width=10, function(x) cvAgg(x), fill=NA, by.column=FALSE)
rollCorr <- rollapplyr(prodMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
yrs <- unique(wideResid$yr)

## Raw trends
pdf(here("outputs/expFigs/prodVar.pdf"), height = 5, width = 6)
par(mfrow=c(2,2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
plot(rollWtdCV ~ yrs, type = "l")
plot(rollSynch ~ yrs, type = "l")
plot(rollAgCV ~ yrs, type = "l")
plot(rollCorr ~ yrs, type = "l")
mtext(side=3, "Variability in log(R/S)", outer=TRUE)
dev.off()

cor(y = rollAgCV[10:61], x = rollWtdCV[10:61])
cor(y = rollAgCV[10:61], x = rollSynch[10:61])
