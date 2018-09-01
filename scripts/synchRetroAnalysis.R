#*************************************************************************************
# synchRetroAnalysis.R
# Date revised: Aug 8, 2018; ONGOING
# Inputs: stock-recruit data csv
# Outputs: figure pdfs
# Explainer: Spin off from synchExplore.R that's trimmed down to focus only on relevant 
# metrics
#*************************************************************************************

require(here); require(synchrony); require(zoo); require(ggplot2); require(dplyr); require(tidyr); require(viridis)

source(here("scripts/synchFunctions.R"))

## Data clean
# SR data
recDat1 <- read.csv(here("/data/sox/fraserRecDatEFF.csv"), stringsAsFactors = FALSE)
recDat2 <- read.csv(here("/data/sox/fraserRecDatTrim.csv"), stringsAsFactors = FALSE) 
recDat <- merge(recDat1, recDat2[, c("stk", "yr", "ets")], by = c("stk", "yr")) #combine ets and eff estimates
recDat <- with(recDat, recDat[order(stk, yr),])
recDat <- recDat %>%
  mutate(prod = log(rec/eff), 
         prodZ = as.numeric(scale(prod, center = TRUE, scale = TRUE)),
         spwnRetYr = yr + 4)
for (i in 1:nrow(recDat)) { #add top-fitting SR model
  stk <- recDat$stk[i]
  if (stk == 1 | stk == 2 | stk == 6 | stk == 8 | stk == 9) {
    recDat$model[i] <- "larkin"
  }	else {
    recDat$model[i] <- "ricker"
  } 
}

# Add lagged EFF abundances for larkin models
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
ts <- recDat %>% 
  group_by(stk) %>% 
  summarise(tsLength = length(!is.na(prod)), firstYr = min(yr), lastYr = max(yr))
selectedStks <- c(1, seq(from=3, to=10, by=1), 18, 19) #stocks w/ time series from 1948
selectedStksShort <- seq(from = 1, to = 19, by =1 )[-c(11, 15, 17)] #stocks w/ time series from 1973
recDatTrim1 <- subset(recDat, !is.na(prod) & !is.na(eff3) & !yr == "2012")
recDatTrim <- recDatTrim1[recDatTrim1$stk %in% selectedStks,]
recDatTrimS <- recDatTrim1[recDatTrim1$stk %in% selectedStksShort,]
recDatTrimS <- recDatTrimS[recDatTrimS$yr > 1972, ]

wideRec <- spread(recDatTrim[,c("stk", "yr", "ets")], stk, ets)
recMat <- as.matrix(wideRec[,-1])
wideProd <- spread(recDatTrim[,c("stk", "yr", "prod")], stk, prod)
prodMat <- as.matrix(wideProd[,-1])
yrs <- unique(wideProd$yr)
wideRecS <- spread(recDatTrimS[,c("stk", "yr", "ets")], stk, ets)
recMatS <- as.matrix(wideRecS[,-1])
wideProdS <- spread(recDatTrimS[,c("stk", "yr", "prod")], stk, prod)
prodMatS <- as.matrix(wideProdS[,-1])
yrsS <- unique(wideProdS$yr)

# Calculate aggregate abundance
aggRecLong <- recDat[recDat$stk %in% selectedStks, ] %>% 
  group_by(spwnRetYr) %>% 
  summarize(aggRec = sum(rec))
aggRecShort <- recDat[recDat$stk %in% selectedStksShort, ] %>% 
  filter(yr > 1972) %>% 
  group_by(spwnRetYr) %>% 
  summarize(aggRec = sum(rec))


# Trim catch data
catchDat <- read.csv(here("/data/sox/fraserCatchDatTrim.csv"), stringsAsFactors = FALSE)
catchDatLong <- catchDat[catchDat$stk %in% selectedStks, ] %>% 
  group_by(yr) %>% 
  summarize(aggCatch = sum(totCatch))
catchDatShort <- catchDat[catchDat$stk %in% selectedStksShort, ] %>% 
  filter(yr > 1972) %>% 
  group_by(yr) %>% 
  summarize(aggCatch = sum(totCatch))


#_________________________________________________________________________
## Frequentist framework
rollCorr <- rollapplyr(prodMat, width=10, function(x) meancorr(x)$obs, fill=NA, by.column=FALSE)
rollWtdCV <- rollapplyr(prodMat, width=10, function(x) wtdCV(x, recMat = recMat), fill=NA, by.column=FALSE)
rollSynch <- rollapplyr(prodMat, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
rollAgCV <- rollapplyr(prodMat, width=10, function(x) cvAgg(x, recMat = recMat), fill=NA, by.column=FALSE)
rollWtdCVShort <- rollapplyr(prodMatS, width=10, function(x) wtdCV(x, recMat = recMatS), fill=NA, by.column=FALSE)
rollSynchShort <- rollapplyr(prodMatS, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
rollAgCVShort <- rollapplyr(prodMatS, width=10, function(x) cvAgg(x, recMat = recMatS), fill=NA, by.column=FALSE)
rollS <- rollapplyr(recMat, width=10, function(x) wtdMean(x), fill=NA, by.column=FALSE)


## Helper objects for plots
stks <- unique(recDatTrim1$stk)
meanP <- recDatTrim1 %>% #mean productivity
  group_by(yr) %>%
  summarise(logRS = mean(prod), zLogRS = mean(prodZ))
colPal <- viridis(n = length(stks), begin = 0, end = 1)

## Plot
pdf(here("figs/Fig1_RetroTrends.pdf"), height = 6, width = 9)
par(mfrow=c(2, 3), oma=c(0,0,0,0)+0.1, mar=c(2,4,1,1), cex.lab = 1.25)
usr <- par( "usr" )
plot(1, type="n", xlab="", ylab="Observed log(R/S)", xlim=range(recDatTrim1$yr), 
     ylim = c(min(recDatTrim1$prodZ), max(recDatTrim1$prodZ)))
for(i in seq_along(stks)) {
  d <- subset(recDatTrim1, stk == stks[i])
  lines(prodZ ~ yr, data= d, type = "l", ylab = "log RS", col = colPal[i])
}
lines(meanP$zLogRS ~ meanP$yr, lwd = 2)
plot(aggRec ~ spwnRetYr, data = aggRecLong, type = "l", ylab = "Aggregate Spawner Abundance", lwd = 1.5)
lines(aggRec ~ spwnRetYr, data = aggRecShort, col = "red")
plot(aggCatch ~ yr, data = catchDatLong, type = "l", ylab = "Aggregate Catch", 
     lwd = 1.5, ylim = c(0, max(catchDatShort$aggCatch)))
lines(aggCatch ~ yr, data = catchDatShort, col = "red")
plot(rollWtdCV ~ yrs, type = "l", ylab = "Weighted Mean Component CV", lwd = 1.5)
lines(rollWtdCVShort ~ yrsS, col = "red")
plot(rollSynch ~ yrs, type = "l", ylab = "Synchrony Index", lwd = 1.5)
lines(rollSynchShort ~ yrsS, col = "red")
plot(rollAgCV ~ yrs, type = "l", ylab = "Aggregate CV", lwd = 1.5)
lines(rollAgCVShort ~ yrsS, col = "red")
dev.off()


## Plot w/ rolling median recruit abundance instead of productivity
aggRecMat <- apply(recMat, 1, sum)
rollMedRec <- rollapplyr(aggRecMat, width=10, function(x) median(x), fill=NA, by.column=FALSE)
aggRecMatShort <- apply(recMatS, 1, sum)
rollMedRecShort <- rollapplyr(aggRecMatShort, width=10, function(x) median(x), fill=NA, by.column=FALSE)


pdf(here("figs/FigX_RetroTrends_wMedSpawners.pdf"), height = 6, width = 8)
par(mfrow=c(2, 2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
usr <- par( "usr" )
plot(rollMedRec ~ yrs, type = "l", ylab = "Median Aggregate Recruit Abundance", lwd = 1.5,
     ylim = c(0, max(rollMedRecShort, na.rm = TRUE)))
lines(rollMedRecShort ~ yrsS, col = "red")
legend("bottomleft", "A)", bty="n") 
plot(rollWtdCV ~ yrs, type = "l", ylab = "Weighted Mean Component CV", lwd = 1.5)
lines(rollWtdCVShort ~ yrsS, col = "red")
legend("bottomleft", "B)", bty="n") 
plot(rollSynch ~ yrs, type = "l", ylab = "Synchrony Index", lwd = 1.5)
lines(rollSynchShort ~ yrsS, col = "red")
legend("bottomleft", "C)", bty="n") 
plot(rollAgCV ~ yrs, type = "l", ylab = "Aggregate CV", lwd = 1.5)
lines(rollAgCVShort ~ yrsS, col = "red")
legend("bottomleft", "D)", bty="n") 
dev.off()


## Repeat analysis but exclude boom/bust years (2009 and 2010 returns)
# recDatTrimNoS <- recDatTrim %>% 
#   filter(!(yr == 2005 | yr == 2006))
# wideRec2 <- spread(recDatTrimNoS[,c("stk", "yr", "ets")], stk, ets)
# recMat2 <- as.matrix(wideRec2[,-1])
# wideProd2 <- spread(recDatTrimNoS[,c("stk", "yr", "prod")], stk, prod)
# prodMat2 <- as.matrix(wideProd2[,-1])
# yrs2 <- unique(wideProd2$yr)
# 
# stks <- unique(recDatTrimNoS$stk)
# meanP <- recDatTrimNoS %>% #mean productivity
#   group_by(yr) %>%
#   summarise(logRS = mean(prod))
# colPal <- viridis(n = length(stks), begin = 0, end = 1)
# 
# rollWtdCV2 <- rollapplyr(prodMat2, width=10, function(x) wtdCV(x, recMat = recMat2), fill=NA, by.column=FALSE)
# rollSynch2 <- rollapplyr(prodMat2, width=10, function(x) community.sync(x)$obs, fill=NA, by.column=FALSE)
# rollAgCV2 <- rollapplyr(prodMat2, width=10, function(x) cvAgg(x, recMat = recMat2), fill=NA, by.column=FALSE)
# 
# ## Plot
# pdf(here("figs/FigY_RetroTrends_noBoomBust.pdf"), height = 6, width = 8)
# par(mfrow=c(2, 2), oma=c(0,0,2,0)+0.1, mar=c(2,4,1,1))
# usr <- par( "usr" )
# plot(1, type="n", xlab="", ylab="Observed log(R/S)", xlim=range(recDatTrimNoS$yr), 
#      ylim = c(min(recDatTrimNoS$prod), max(recDatTrimNoS$prod)))
# for(i in seq_along(stks)) {
#   d <- subset(recDatTrimNoS, stk == stks[i])
#   lines(prod ~ yr, data= d, type = "l", ylab = "log RS", col = colPal[i])
# }
# lines(logRS ~ yr, data = meanP, lwd = 2)
# legend("bottomleft", "A)", bty="n") 
# plot(rollWtdCV2 ~ yrs2, type = "l", ylab = "Weighted Mean Component CV", lwd = 1.5)
# legend("bottomleft", "B)", bty="n") 
# plot(rollSynch2 ~ yrs2, type = "l", ylab = "Synchrony Index", lwd = 1.5)
# legend("bottomleft", "C)", bty="n") 
# plot(rollAgCV2 ~ yrs2, type = "l", ylab = "Aggregate CV", lwd = 1.5)
# legend("bottomleft", "D)", bty="n") 
# dev.off()
