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

# setwd("C:/github/synchSalmon/")
setwd("/Users/cam/github/synchSalmon") #Cam's Mac wd

require(here); require(synchrony); require(zoo); require(ggplot2); require(dplyr); require(tidyr)

source(here("scripts/synchFunctions.R"))

## Convert long dataframe to wide matrix
recDat <- read.csv(here("data/sox/fraserRecDatTrim.csv"), stringsAsFactors=F)
recDat <- with(recDat, recDat[order(stk, yr),])
ts <- recDat %>% group_by(stk) %>% summarise(tsLength=length(!is.na(ets)), firstYr=min(yr))
selectedStks <- c(1, seq(from=3, to=11, by=1), 18, 19)
recDatTrim <- recDat[recDat$stk %in% selectedStks,]
# recDat <- recDat[!recDat$stk == "15",] #remove stk with short TS
# recDat <- recDat[!recDat$yr < "1973",] #remove earlier years to ensure time series balanced

wideRec <- spread(recDatTrim[,c(1:3)], stk, ets)
recMat <- as.matrix(wideRec[,-1])


#_________________________________________________________________________
## Frequentist framework
rollWtdCV <- rollapplyr(recMat, width=10, function(x) wtdCV(x), fill=NA, by.column=FALSE)
rollSynch <- rollapplyr(recMat, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
rollAgCV <- rollapplyr(recMat, width=10, function(x) cvAgg(x), fill=NA, by.column=FALSE)
rollCorr <- rollapplyr(recMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
yrs <- unique(wideRec$yr)

## Raw trends
par(mfrow=c(2,2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
plot(rollWtdCV ~ yrs)
plot(rollSynch ~ yrs)
plot(rollAgCV ~ yrs)
plot(rollCorr ~ yrs)
mtext(side=3, "Fraser Sockeye Metapopulation Dynamics", outer=TRUE)



