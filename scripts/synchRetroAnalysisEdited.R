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
selectedStks <- c(1, seq(from=3, to=10, by=1), 18, 19) #stocks w/ time series from 1948
recDatTrim <- recDatTrim1 %>% 
  filter(stk %in% selectedStks)
yrs <- unique(recDatTrim$yr)
spwnMat <- recDatTrim %>% 
  select(stk, yr, ets) %>% 
  spread(stk, ets) %>% 
  select(-yr) %>% 
  as.matrix()
recMat <- recDatTrim %>% 
  select(stk, yr, rec) %>% 
  spread(stk, rec) %>% 
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
aggRecLong <- recDatTrim %>% 
  group_by(spwnRetYr) %>% 
  summarize(aggRec = sum(rec)) %>% 
  mutate(ts = "long")
prodBySDat <- data.frame(year = yrs,
                       wtdCV = rollapplyr(prodMat, width = 12, 
                                              function(x) wtdCV(x, recMat = recMat),
                                              fill = NA, by.column = FALSE),
                       synch = rollapplyr(prodMat, width = 12, 
                                              function(x) community.sync(x)$obs, 
                                              fill = NA, by.column = FALSE),
                       agCV = rollapplyr(prodMat, width = 12, 
                                             function(x) cvAgg(x, recMat = recMat),
                                             fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "prod", ts = "long", weight = "s") 
residBySDat <- data.frame(year = yrs, #repeat w/ resid data
                      wtdCV = rollapplyr(residMat, width = 12, 
                                         function(x) wtdCV(x, recMat = recMat),
                                         fill = NA, by.column = FALSE),
                      synch = rollapplyr(residMat, width = 12, 
                                         function(x) community.sync(x)$obs, 
                                         fill = NA, by.column = FALSE),
                      agCV = rollapplyr(residMat, width = 12, 
                                        function(x) cvAgg(x, recMat = recMat),
                                        fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "resid", ts = "long", weight = "s")
prodByPDat <- data.frame(year = yrs,
                         wtdCV = rollapplyr(prodMat, width = 12, 
                                               function(x) wtdCV(x),
                                               fill = NA, by.column = FALSE),
                         synch = rollapplyr(prodMat, width = 12, 
                                            function(x) community.sync(x)$obs, 
                                            fill = NA, by.column = FALSE),
                         agCV = rollapplyr(prodMat, width = 12, 
                                              function(x) cvAgg(x),
                                              fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "prod", ts = "long", weight = "p") 
residByPDat <- data.frame(year = yrs, #repeat w/ resid data
                          wtdCV = rollapplyr(residMat, width = 12, 
                                                function(x) wtdCV(x),
                                                fill = NA, by.column = FALSE),
                          synch = rollapplyr(residMat, width = 12, 
                                             function(x) community.sync(x)$obs, 
                                             fill = NA, by.column = FALSE),
                          agCV = rollapplyr(residMat, width = 12, 
                                               function(x) cvAgg(x),
                                               fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "resid", ts = "long", weight = "p") 
prodDat <- data.frame(year = yrs,
                         wtdCV = rollapplyr(prodMat, width = 12, 
                                            function(x) wtdCV(x),
                                            fill = NA, by.column = FALSE),
                         synch = rollapplyr(prodMat, width = 12, 
                                            function(x) community.sync(x)$obs, 
                                            fill = NA, by.column = FALSE),
                         agCV = rollapplyr(prodMat, width = 12, 
                                           function(x) cvAgg(x),
                                           fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "prod", ts = "long", weight = "none") 
residDat <- data.frame(year = yrs, #repeat w/ resid data
                          wtdCV = rollapplyr(residMat, width = 12, 
                                             function(x) wtdCV(x, weight = FALSE),
                                             fill = NA, by.column = FALSE),
                          synch = rollapplyr(residMat, width = 12, 
                                             function(x) community.sync(x)$obs, 
                                             fill = NA, by.column = FALSE),
                          agCV = rollapplyr(residMat, width = 12, 
                                            function(x) cvAgg(x, weight = FALSE),
                                            fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "resid", ts = "long", weight = "none") 


# Changes in mean pairwise corrleations through time
# cbind(yrs, 
#       rollapplyr(prodMat, width = 12, function(x) meancorr(x)$obs, 
#                  fill = NA, by.column = FALSE),
#       rollapplyr(residMat, width = 12, function(x) meancorr(x)$obs, 
#                  fill = NA, by.column = FALSE),
#       rollapplyr(logResidMat, width = 12, function(x) meancorr(x)$obs, 
#                  fill = NA, by.column = FALSE))



# Combine data
dat <- rbind(prodBySDat, residBySDat, prodByPDat, residByPDat, 
             prodDat, residDat) %>%
  mutate(var = as.factor(var), data = as.factor(data), ts = as.factor(ts),
         weight = as.factor(weight))
# %>% 
#   mutate(var = recode(var, "agCV" ="Aggregate CV", "rollSynch" = "Synchrony",
#                       "wtdCV" = "Component CV", .default = levels(var)))
# saveRDS(dat, file = here("data/retroSynchData.Rdata"))


# Trim catch data
catchDat <- read.csv(here("/data/sox/fraserCatchDatTrim.csv"), 
                     stringsAsFactors = FALSE)
catchDatLong <- catchDat %>% 
  filter(stk %in% selectedStks) %>% 
  group_by(yr) %>% 
  summarize(catch = sum(totCatch)) %>% 
  mutate(ts = as.factor("long"))

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

aggSpwnPlot <- ggplot(aggRecLong, aes(x = spwnRetYr, y = aggRec)) + 
  geom_line() +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Aggregate Spawner Abundance") 

aggCatchPlot <- ggplot(aggCatch, aes(x = yr, y = catch)) + 
  geom_line() +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Aggregate Catch") 

pdf(here("figs/productivityTrends.pdf"), height = 6, width = 8)
ggplot(dat %>% filter(data == "prod"), aes(x = year, y = index, colour = weight)) + 
  geom_line() +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Index", title = "R/S Trends") +
  facet_wrap(~var)
ggplot(dat %>% filter(data == "resid"), aes(x = year, y = index, colour = weight)) + 
  geom_line() +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Index", title = "Residual Trends") +
  facet_wrap(~var)
dev.off()

t <- dat %>% filter(weight == "none", var == "wtdCV", data == "prod") %>% select(index)
t2 <- dat %>%  filter(weight == "p", var == "wtdCV", data == "prod") %>% select(index)
plot(unlist(t) ~ unlist(t2))
#______________________________________________________________________________
# repeat above, but for more stocks and shorter time period
# selectedStksShort <- seq(from = 1, to = 19, by =1 )[-c(11, 15, 17)] 
# recDatTrimS <- recDatTrim1 %>%  
#   filter(stk %in% selectedStksShort, yr > 1972)
# yrsShort <- unique(recDatTrimS$yr)
# spwnMatShort <- recDatTrimS %>% 
#   select(stk, yr, ets) %>% 
#   spread(stk, ets) %>% 
#   select(-yr) %>% 
#   as.matrix()
# prodMatShort <- recDatTrimS %>% 
#   select(stk, yr, prod) %>% 
#   spread(stk, prod) %>% 
#   select(-yr) %>% 
#   as.matrix()
# residMatShort <- recDatTrimS %>% 
#   select(stk, yr, modResid) %>% 
#   spread(stk, modResid) %>% 
#   select(-yr) %>% 
#   as.matrix()
# aggRecShort <- recDatTrimS %>% 
#   group_by(spwnRetYr) %>% 
#   summarize(aggRec = sum(rec)) %>% 
#   mutate(ts = "short")
# 
# prodDatS <- data.frame(year = yrsShort,
#                        wtdCV = rollapplyr(prodMatShort, width = 12, 
#                                           function(x) wtdCV(x, recMat = spwnMatShort),
#                                           fill = NA, by.column = FALSE),
#                        synch = rollapplyr(prodMatShort, width = 12, 
#                                           function(x) community.sync(x)$obs, 
#                                           fill = NA, by.column = FALSE),
#                        agCV = rollapplyr(prodMatShort, width = 12, 
#                                          function(x) cvAgg(x, recMat = spwnMatShort),
#                                          fill = NA, by.column = FALSE)) %>% 
#   gather(key = var, value = index, 2:4) %>% 
#   mutate(data = "prod", ts = "short") 
# residDatS <- data.frame(year = yrsShort, #repeat w/ resid data
#                         wtdCV = rollapplyr(residMatShort, width = 12, 
#                                            function(x) wtdCV(x, recMat = spwnMatShort),
#                                            fill = NA, by.column = FALSE),
#                         synch = rollapplyr(residMatShort, width = 12, 
#                                            function(x) community.sync(x)$obs, 
#                                            fill = NA, by.column = FALSE),
#                         agCV = rollapplyr(residMatShort, width = 12, 
#                                           function(x) cvAgg(x, recMat = spwnMatShort),
#                                           fill = NA, by.column = FALSE)) %>% 
#   gather(key = var, value = index, 2:4) %>% 
#   mutate(data = "resid", ts = "short")
# 
# # cbind(yrsShort, 
# #       rollapplyr(prodMatShort, width = 12, function(x) meancorr(x)$obs, 
# #                  fill = NA, by.column = FALSE),
# #       rollapplyr(residMatShort, width = 12, function(x) meancorr(x)$obs, 
# #                  fill = NA, by.column = FALSE))
catchDatShort <- catchDat %>% 
  filter(stk %in% selectedStksShort, yr > 1972) %>% 
  group_by(yr) %>% 
  summarize(catch = sum(totCatch)) %>% 
  mutate(ts = as.factor("short"))
aggCatch <- rbind(catchDatLong, catchDatShort)

