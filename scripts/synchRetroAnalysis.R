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
require(tidyr); require(viridis); require(ggpubr)

source(here("scripts/synchFunctions.R"))

select <- dplyr::select

## Data clean
# SR data
recDat1 <- read.csv(here("/data/sox/fraserRecDatEFF.csv"), stringsAsFactors = FALSE)
recDat2 <- read.csv(here("/data/sox/fraserRecDatTrim.csv"), stringsAsFactors = FALSE) 
recDat <- merge(recDat1, recDat2[, c("stk", "yr", "ets")], by = c("stk", "yr")) #combine ets and eff estimates
recDat <- with(recDat, recDat[order(stk, yr),])
recDat <- recDat %>%
  mutate(prod = (rec/ets),
         logProd = log(rec/ets),
         prodZ = NA,
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
    d$prodZ <- as.numeric(scale(d$logProd, center = TRUE, scale = TRUE))
  }
  recDat[recDat$stk == stkIndex[j], c("prodZ", "ets1", "ets2", "ets3")] <- d[, c("prodZ", "ets1", "ets2", "ets3")]
}
recDatTrim1 <- subset(recDat, !is.na(prod) & !is.na(ets3) & !yr == "2012")

# saveRDS(recDatTrim1, file = here("data", "generated", "recDatTrim1.rds"))
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

recDatTrim1 %>% 
  group_by(stk) %>% 
  summarize(mean(rec), mean(prod))

# Trim and convert to matrix
selectedStks <- c(1, seq(from=3, to=10, by=1), 18, 19) #stocks w/ time series from 1948
recDatTrim <- recDatTrim1 %>% 
  filter(stk %in% selectedStks)
yrs <- unique(recDatTrim$yr)
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
aggRec <- recDatTrim %>% 
  group_by(spwnRetYr) %>% 
  summarize(aggRec = sum(rec)) %>% 
  mutate(ts = "long")

saveRDS(recMat, file = here("data", "generated", "recMat.rds"))


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
  select(stk, yr, prodZ, logProd, modResid, logResid) %>% 
  mutate(stk = as.factor(stk))

# Import stan model outputs
m <- readRDS(here("data/generated/stanSynchModOut.rds")) #window size = 12
varNames <- c("cv_s", "phi", "cv_c")
parList <- lapply(seq_along(varNames), function(x) {
  d <- plyr::ldply(m, function(y) 
    quantile(y[[varNames[x]]], probs = c(0.1, 0.25, 0.5, 0.75, 0.9))) %>% 
    dplyr::rename(low = "10%", med = "50%", high = "90%") %>% 
    mutate(year = yrs)
})
names(parList) <- c("cvC", "synch", "cvA")

#_________________________________________________________________________
## Plot
colPal <- viridis(n = length(unique(rawDat$stk)), begin = 0, end = 1)
axisSize = 14

rawProdPlot <- ggplot(rawDat, aes(x = yr, y = logProd, colour = stk)) + 
  labs(x = "", y = "Standardized Productivity") + 
  geom_line(size = 0.75) +
  scale_color_manual(values = colPal, guide = FALSE) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", size = 1.25)  +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize),
        axis.title = element_text(size = axisSize))
aggRecPlot <- ggplot(aggRec, aes(x = spwnRetYr, y = aggRec)) +
  geom_line(size = 1.25) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize)) +
  labs(x = "", y = "Aggregate Recruit Abundance")
aggCatchPlot <- ggplot(catchDatLong, aes(x = yr, y = catch)) + 
  geom_line(size = 1.25) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize),
        axis.title = element_text(size = axisSize)) +
  labs(x = "", y = "Aggregate Catch") 
compCVPlot <- ggplot(parList$cvC, aes(x = year, y = med)) + 
  geom_line(size = 1.25) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha=0.2) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize),
        axis.title = element_text(size = axisSize)) +
  labs(x = "", y = "Component Variability") 
synchPlot <- ggplot(parList$synch, aes(x = year, y = med)) + 
  geom_line(size = 1.25) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha=0.2) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize),
        axis.title = element_text(size = axisSize)) +
  labs(x = "", y = "Synchrony")
agCVPlot <- ggplot(parList$cvA, aes(x = year, y = med)) + 
  geom_line(size = 1.25) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha=0.2) +
  theme_sleekX() +
  theme(axis.text = element_text(size = 0.9 * axisSize),
        axis.title = element_text(size = axisSize)) +
  labs(x = "", y = "Aggregate Variability")

png(here("figs/Fig1New_Retro.png"), height = 6, width = 7.5,
    units = "in", res = 200)
ggarrange(rawProdPlot, aggRecPlot, aggCatchPlot, compCVPlot, synchPlot, 
          agCVPlot, nrow = 2, ncol = 3)
dev.off()


# w/ faceting
# pdf(here("figs/productivityTrends.pdf"), height = 6, width = 8)
# ggplot(dat %>% filter(data == "prod"), aes(x = year, y = index)) + 
#   geom_line() +
#   theme_sleekX() +
#   theme(axis.text = element_text(size = 0.9 * axisSize)) +
#   labs(x = "", y = "Index", title = "R/S Trends") +
#   facet_wrap(~var)
# ggplot(dat %>% filter(data == "resid"), aes(x = year, y = index, colour = weight)) + 
#   geom_line() +
#   theme_sleekX() +
#   theme(axis.text = element_text(size = 0.9 * axisSize)) +
#   labs(x = "", y = "Index", title = "Residual Trends") +
#   facet_wrap(~var)
# dev.off()


#### Old frequentist version
prodDat <- data.frame(year = yrs,
                      wtdCV = rollapplyr(prodMat, width = 12, 
                                         function(x) wtdCV(x, weightMat = recMat),
                                         fill = NA, by.column = FALSE),
                      synch = rollapplyr(prodMat, width = 12, 
                                         function(x) community.sync(x)$obs, 
                                         fill = NA, by.column = FALSE),
                      agCV = rollapplyr(prodMat, width = 12, 
                                        function(x) cvAgg(x, weightMat = recMat),
                                        fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "prod", ts = "long", weight = "s") 
residDat <- data.frame(year = yrs, #repeat w/ resid data
                       wtdCV = rollapplyr(residMat, width = 12, 
                                          function(x) wtdCV(x, weightMat = recMat),
                                          fill = NA, by.column = FALSE),
                       synch = rollapplyr(residMat, width = 12, 
                                          function(x) community.sync(x)$obs, 
                                          fill = NA, by.column = FALSE),
                       agCV = rollapplyr(residMat, width = 12, 
                                         function(x) cvAgg(x, weightMat = recMat),
                                         fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "resid", ts = "long", weight = "s")
rDat <- data.frame(year = yrs,
                   wtdCV = rollapplyr(recMat, width = 12, 
                                      function(x) wtdCV(x),
                                      fill = NA, by.column = FALSE),
                   synch = rollapplyr(recMat, width = 12, 
                                      function(x) community.sync(x)$obs, 
                                      fill = NA, by.column = FALSE),
                   agCV = rollapplyr(recMat, width = 12, 
                                     function(x) cvAgg(x),
                                     fill = NA, by.column = FALSE)) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "rec", ts = "long", weight = "s")


# Combine data
dat <- rbind(rDat, prodDat, residDat) %>%
  mutate(var = as.factor(var), data = as.factor(data), ts = as.factor(ts),
         weight = as.factor(weight)) %>%
  mutate(var = recode(var, "agCV" ="Aggregate CV", "synch" = "Synchrony",
                      "wtdCV" = "Component CV", .default = levels(var)))

