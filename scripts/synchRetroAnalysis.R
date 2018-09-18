#*************************************************************************************
# synchRetroAnalysis.R
# Date revised: Aug 8, 2018; ONGOING
# Inputs: stock-recruit data csv
# Outputs: figure pdfs
# Explainer: Spin off from synchExplore.R that's trimmed down to focus only on relevant 
# metrics; adjusted Sep 12 to focus on exp(residuals) or non-logged R/S; use ETS to be
# consistent w/ forward simulation and FRSSI fitted SR models
#*************************************************************************************

require(here); require(synchrony); require(zoo); require(ggplot2); require(dplyr)
require(tidyr); require(viridis)

source(here("scripts/synchFunctions.R"))

## Data clean
# SR data
recDat1 <- read.csv(here("/data/sox/fraserRecDatEFF.csv"), stringsAsFactors = FALSE)
recDat2 <- read.csv(here("/data/sox/fraserRecDatTrim.csv"), stringsAsFactors = FALSE) 
recDat <- merge(recDat1, recDat2[, c("stk", "yr", "ets")], by = c("stk", "yr")) #combine ets and eff estimates
recDat <- with(recDat, recDat[order(stk, yr),])
recDat <- recDat %>%
  mutate(prod = (rec/ets),
         logProd = log(rec/ets),
         prodZ = as.numeric(scale(prod, center = TRUE, scale = TRUE)),
         spwnRetYr = yr + 4)
for (i in 1:nrow(recDat)) { #add top-fitting SR model
  stkN <- recDat$stk[i]
  if (stkN == 1 | stkN == 2 | stkN == 6 | stkN == 8 | stkN == 9) {
    recDat$model[i] <- "larkin"
  }	else {
    recDat$model[i] <- "ricker"
  } 
}

# Add lagged EFF abundances for larkin models
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
recDatTrim1 <- subset(recDat, !is.na(prod) & !is.na(ets3) & !yr == "2012")

saveRDS(recDatTrim1, file = here("data", "generated", "recDatTrim1.rds"))
# Fit single-stock SR models
stkIndex <- unique(recDatTrim1$stk)
residVec <- NULL
for(j in seq_along(stkIndex)) {
  d <- subset(recDatTrim1, stk == stkIndex[j])
  mod <- unique(d$model)
  if (mod == "ricker") {
    srMod <- lm(logProd ~ ets, data = d)
  }
  if (mod == "larkin") {
    srMod <- lm(logProd ~ ets + ets1 + ets2 + ets3, data = d)
  }
  residVec <- c(residVec, resid(srMod))
}
recDatTrim1$logResid <- residVec
recDatTrim1$modResid <- exp(residVec)

# Trim and convert to matrix
# recDatTrim1 %>% 
#   group_by(stk) %>% 
#   summarise(tsLength = length(!is.na(prod)), firstYr = min(yr), lastYr = max(yr))
selectedStks <- c(1, seq(from=3, to=10, by=1), 18, 19) #stocks w/ time series from 1948
recDatTrim <- recDatTrim1 %>% 
  filter(stk %in% selectedStks)
yrs <- unique(recDatTrim$yr)
spwnMat <- recDatTrim %>% 
  select(stk, yr, ets) %>% 
  spread(stk, ets) %>% 
  select(-yr) %>% 
  as.matrix()
prodMat <- recDatTrim %>% 
  select(stk, yr, prod) %>% 
  spread(stk, prod) %>% 
  select(-yr) %>% 
  as.matrix()
residMat <- recDatTrim %>% 
  select(stk, yr, modResid) %>% 
  spread(stk, modResid) %>% 
  select(-yr) %>% 
  as.matrix()
logResidMat <- recDatTrim %>% 
  select(stk, yr, logResid) %>% 
  spread(stk, logResid) %>% 
  select(-yr) %>% 
  as.matrix()
aggRecLong <- recDatTrim %>% 
  group_by(spwnRetYr) %>% 
  summarize(aggRec = sum(rec)) %>% 
  mutate(ts = "long")

saveRDS(prodMat, file = here("data", "generated", "prodMat.rds"))
prodDat <- data.frame(year = yrs,
                       wtdCV = rollapplyr(prodMat, width = 10, 
                                              function(x) wtdCV(x, recMat = spwnMat),
                                              fill = NA, by.column = FALSE),
                       synch = rollapplyr(prodMat, width = 10, 
                                              function(x) community.sync(x)$obs, 
                                              fill = NA, by.column = FALSE),
                       agCV = rollapplyr(prodMat, width = 10, 
                                             function(x) cvAgg(x, recMat = spwnMat),
                                             fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "prod", ts = "long") 
residDat <- data.frame(year = yrs, #repeat w/ resid data
                      wtdCV = rollapplyr(residMat, width = 10, 
                                         function(x) wtdCV(x, recMat = spwnMat),
                                         fill = NA, by.column = FALSE),
                      synch = rollapplyr(residMat, width = 10, 
                                         function(x) community.sync(x)$obs, 
                                         fill = NA, by.column = FALSE),
                      agCV = rollapplyr(residMat, width = 10, 
                                        function(x) cvAgg(x, recMat = spwnMat),
                                        fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "resid", ts = "long") 

# Changes in mean pairwise corrleations through time
cbind(yrs, 
      rollapplyr(prodMat, width = 10, function(x) meancorr(x)$obs, 
                 fill = NA, by.column = FALSE),
      rollapplyr(residMat, width = 10, function(x) meancorr(x)$obs, 
                 fill = NA, by.column = FALSE),
      rollapplyr(logResidMat, width = 10, function(x) meancorr(x)$obs, 
                 fill = NA, by.column = FALSE))


# repeat above, but for more stocks and shorter time period
selectedStksShort <- seq(from = 1, to = 19, by =1 )[-c(11, 15, 17)] 
recDatTrimS <- recDatTrim1 %>%  
  filter(stk %in% selectedStksShort, yr > 1972)
yrsShort <- unique(recDatTrimS$yr)
spwnMatShort <- recDatTrimS %>% 
  select(stk, yr, ets) %>% 
  spread(stk, ets) %>% 
  select(-yr) %>% 
  as.matrix()
prodMatShort <- recDatTrimS %>% 
  select(stk, yr, prod) %>% 
  spread(stk, prod) %>% 
  select(-yr) %>% 
  as.matrix()
residMatShort <- recDatTrimS %>% 
  select(stk, yr, modResid) %>% 
  spread(stk, modResid) %>% 
  select(-yr) %>% 
  as.matrix()
aggRecShort <- recDatTrimS %>% 
  group_by(spwnRetYr) %>% 
  summarize(aggRec = sum(rec)) %>% 
  mutate(ts = "short")

prodDatS <- data.frame(year = yrsShort,
                      wtdCV = rollapplyr(prodMatShort, width = 10, 
                                         function(x) wtdCV(x, recMat = spwnMatShort),
                                         fill = NA, by.column = FALSE),
                      synch = rollapplyr(prodMatShort, width = 10, 
                                         function(x) community.sync(x)$obs, 
                                         fill = NA, by.column = FALSE),
                      agCV = rollapplyr(prodMatShort, width = 10, 
                                        function(x) cvAgg(x, recMat = spwnMatShort),
                                        fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "prod", ts = "short") 
residDatS <- data.frame(year = yrsShort, #repeat w/ resid data
                       wtdCV = rollapplyr(residMatShort, width = 10, 
                                          function(x) wtdCV(x, recMat = spwnMatShort),
                                          fill = NA, by.column = FALSE),
                       synch = rollapplyr(residMatShort, width = 10, 
                                          function(x) community.sync(x)$obs, 
                                          fill = NA, by.column = FALSE),
                       agCV = rollapplyr(residMatShort, width = 10, 
                                         function(x) cvAgg(x, recMat = spwnMatShort),
                                         fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "resid", ts = "short")

cbind(yrsShort, 
      rollapplyr(prodMatShort, width = 10, function(x) meancorr(x)$obs, 
                 fill = NA, by.column = FALSE),
      rollapplyr(residMatShort, width = 10, function(x) meancorr(x)$obs, 
                 fill = NA, by.column = FALSE))

# Combine data
dat <- rbind(prodDat, residDat, prodDatS, residDatS) %>%
  mutate(var = as.factor(var), data = as.factor(data), ts = as.factor(ts)) %>% 
  mutate(var = recode(var, "rollAgCV" ="Aggregate CV", "rollSynch" = "Synchrony",
                      "rollWtdCV" = "Component CV", .default = levels(var)))
# saveRDS(dat, file = here("data/retroSynchData.Rdata"))


aggRec <- rbind(aggRecLong, aggRecShort) %>% 
  mutate(ts = as.factor(ts))

# Trim catch data
catchDat <- read.csv(here("/data/sox/fraserCatchDatTrim.csv"), 
                     stringsAsFactors = FALSE)
catchDatLong <- catchDat %>% 
  filter(stk %in% selectedStks) %>% 
  group_by(yr) %>% 
  summarize(catch = sum(totCatch)) %>% 
  mutate(ts = as.factor("long"))
catchDatShort <- catchDat %>% 
  filter(stk %in% selectedStksShort, yr > 1972) %>% 
  group_by(yr) %>% 
  summarize(catch = sum(totCatch)) %>% 
  mutate(ts = as.factor("short"))
aggCatch <- rbind(catchDatLong, catchDatShort)

# Trim productivity residual data
rawDat <- recDatTrim1 %>%
  select(stk, yr, prod, logProd, modResid, logResid) %>% 
  mutate(stk = as.factor(stk))

#_________________________________________________________________________
## Plot
xLab = "Year"
colPal <- viridis(n = length(unique(rawDat$stk)), begin = 0, end = 1)
axisSize = 12

rawProdPlot <- ggplot(rawDat, aes(x = yr, y = logResid, colour = stk)) + 
  labs(x = "", y = "Observed Productivitiy") + 
  geom_line(size = 1) +
  scale_color_manual(values = colPal, guide = FALSE) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", size = 1.75) +
  theme_sleekX()

aggSpwnPlot <- ggplot(aggRec, aes(x = spwnRetYr, y = aggRec, colour = ts, 
                                  size = ts)) + 
  geom_line() +
  scale_size_manual(values = c(1.25, 1), guide = FALSE) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Aggregate Spawner Abundance") 

aggCatchPlot <- ggplot(aggCatch, aes(x = yr, y = catch, colour = ts)) + 
  geom_line() +
  scale_size_manual(values = c(1.25, 1), guide = FALSE) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Aggregate Catch") 

pdf(here("figs/productivityTrends.pdf"), height = 6, width = 8)
ggplot(dat %>% filter(data == "prod"), aes(x = year, y = index, colour = ts)) + 
  geom_line() +
  scale_size_manual(values = c(1.25, 1), guide = FALSE) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Index", title = "R/S Trends") +
  facet_wrap(~var)
ggplot(dat %>% filter(data == "resid"), aes(x = year, y = index, colour = ts)) + 
  geom_line() +
  scale_size_manual(values = c(1.25, 1),  guide = FALSE) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Index", title = "Residual Trends") +
  facet_wrap(~var)
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



# wideRec <- spread(recDat[,c("stk", "yr", "ets")], stk, ets)
# recMat <- as.matrix(wideRec[,-1])
# wideProd <- spread(recDat[,c("stk", "yr", "prod")], stk, prod)
# prodMat <- as.matrix(wideProd[,-1])
# wideProd2 <- spread(recDat[,c("stk", "yr", "prod2")], stk, prod2)
# prodMat2 <- as.matrix(wideProd2[,-1])
# 
# temp <- recDat %>% filter(stk == 6) %>% select(yr, ets, eff, prod, prod2)
# rollapplyr(temp$prod, width = 10, function(x) 
#   sqrt(var(x, na.rm = TRUE))/ mean(x, na.rm = TRUE))
# rollapplyr(temp$prod2, width = 10, function(x) 
#   sqrt(var(x, na.rm = TRUE))/ mean(x, na.rm = TRUE))
# temp$prod - temp$prod2



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
