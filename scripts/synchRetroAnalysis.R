#*******************************************************************************
# synchRetroAnalysis.R
# Date revised: Aug 8, 2018; ONGOING
# Inputs: stock-recruit data csv
# Outputs: figure pdfs
# Explainer: Spin off from synchExplore.R that's trimmed down to focus only on 
# relevantmetrics; adjusted Sep 12 to focus on exp(residuals) or non-logged R/S; 
# use ETS to be consistent w/ forward simulation and FRSSI fitted SR models
#*******************************************************************************

require(here); require(synchrony); require(zoo); require(ggplot2); 
require(tidyverse); require(viridis); require(ggpubr); require(samSim)

select <- dplyr::select

## Data clean
# SR data
recDat <- read.csv(here("/data/sox/fraserRecDat.csv"), 
                    stringsAsFactors = FALSE) 
recDat <- recDat %>%
  mutate(prod = (rec/ets),
         logProd = log(rec/ets),
         prodZ = NA,
         spwnRetYr = yr + 4)
for (i in 1:nrow(recDat)) { #add top-fitting SR model
  stkN <- recDat$stk[i]
  if (stkN == 1 | stkN == 2 | stkN == 6 | stkN == 815 | stkN == 9) {
    recDat$model[i] <- "larkin"
  }	else {
    recDat$model[i] <- "ricker"
  } 
}

# Add rolling estimates and lagged EFF abundances for larkin models
recDat$rollProd <- NA
recDat$ets1 <- NA
recDat$ets2 <- NA
recDat$ets3 <- NA
stkIndex <- unique(recDat$stk)
for(j in seq_along(stkIndex)) {
  d <- subset(recDat, stk == stkIndex[j])
  d$rollProd <- rollmeanr(d$logProd, 4, fill = NA)
  for (i in 1:nrow(d)) { #add top-fitting SR model
    d$ets1[i] <- ifelse(i < 1, NA, d$ets[i - 1])
    d$ets2[i] <- ifelse(i < 2, NA, d$ets[i - 2])
    d$ets3[i] <- ifelse(i < 3, NA, d$ets[i - 3])
    d$prodZ <- as.numeric(scale(d$logProd, center = TRUE, scale = TRUE))
  }
  recDat[recDat$stk == stkIndex[j], c("rollProd", "prodZ", "ets1", "ets2", 
                                      "ets3")] <- d[, c("rollProd", "prodZ", 
                                                        "ets1", "ets2", "ets3")]
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
  rollPVec <- c(rollPVec, resid(srMod))
}
recDatTrim1$logResid <- residVec
recDatTrim1$modResid <- exp(residVec)

# Trim and convert to matrix w/ stocks w/ time series from 1948
longStks <- recDatTrim1 %>% 
  group_by(stk) %>% 
  summarize(startDate = min(yr)) %>% 
  filter(!(startDate > 1951 | stk == "11")) %>% 
  select(stk) %>% 
  as.vector()
recDatTrim <- recDatTrim1 %>% 
  dplyr::filter(stk %in% longStks$stk)
# dplyr::filter(yr > 1975, !stk == "11", !stk == "815") #more stocks, shorter TS

yrs <- unique(recDatTrim$yr)
recMat <- recDatTrim %>% 
  select(stk, yr, rec) %>% 
  spread(stk, rec) %>% 
  select(-yr) %>% 
  as.matrix()
saveRDS(recMat, file = here("outputs", "generatedData", "recMat.rds"))

# Trim catch data
catchDat <- read.csv(here("/data/sox/fraserCatchDatTrim.csv"), 
                     stringsAsFactors = FALSE)
catchDatLong <- catchDat %>% 
  filter(stk %in% longStks$stk) %>% 
  group_by(yr) %>% 
  summarize(catch = sum(totalCatch), 
            ret = sum(totalRun)) %>% 
  mutate(rollCatch = rollmeanr(catch, 4, fill = NA),
         rollRet = rollmeanr(ret, 4, fill = NA))

# Trim productivity residual data
rawDat <- recDatTrim1 %>%
  select(stk, yr, rollProd, prodZ, logProd, modResid, logResid) %>% 
  filter(stk %in% longStks$stk) %>% 
  mutate(stk = as.factor(stk))

# Import TMB model outputs
modOut <- readRDS(here("outputs", "generatedData", "tmbSynchEst.rds"))
modOut$broodYr <- rep(yrs[12:61], times = 3)

#note that in Sean's output cv_s is "species" and cv_c is "community" following
#Thib and Connolly 2013
retroCVc <- modOut %>% 
  filter(term == "log_cv_s")  
retroPhi <- modOut %>% 
  filter(term == "logit_phi")
retroCVa <- modOut %>% 
  filter(term == "log_cv_c")

#_________________________________________________________________________
## Plot
stkN <- length(unique(rawDat$stk))

labDat <- data.frame(lab = c("a)", "b)", "c)", "d)", "e)", "f)"),
                     var = c("prod", "rec", "catch", "cv", "synch", "aggCV"))
plotYrs <- c(min(rawDat$yr), max(catchDatLong$yr))

rawProdPlot <- ggplot(rawDat, aes(x = yr, y = rollProd, col = stk)) + 
  labs(x = "", y = "log(Recruits/Spawner)") + 
  ylim(min(rawDat$rollProd),  max(rawDat$rollProd)) +
  geom_line(size = 0.45) +
  xlim(plotYrs) +
  scale_color_manual(values = rep("grey", length.out = stkN), guide = FALSE) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", size = 1.25)  +
  theme_sleekX(position = "top", axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "prod"),
            mapping = aes(x = min(plotYrs), y = max(rawDat$rollProd), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4)
aggRetPlot <- ggplot(catchDatLong, aes(x = yr, y = rollRet)) +
  geom_line(size = 1.25) +
  xlim(plotYrs) +
  theme_sleekX(position = "top", axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "rec"),
            mapping = aes(x = min(plotYrs), 
                          y = max(catchDatLong$rollRet, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Aggregate Return")
aggCatchPlot <- ggplot(catchDatLong, aes(x = yr, y = rollCatch)) + 
  geom_line(size = 1.15) +
  theme_sleekX(position = "top", axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  xlim(plotYrs) +
  geom_text(data = labDat %>% filter(var == "catch"),
            mapping = aes(x = min(plotYrs), 
                          y = max(catchDatLong$rollCatch, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Aggregate Catch") 
compCVPlot <- ggplot(retroCVc, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.15) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "cv"),
            mapping = aes(x = min(plotYrs), 
                          y = max(retroCVc$upr, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Component Variability") 
synchPlot <- ggplot(retroPhi, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.25) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "synch"),
            mapping = aes(x = min(plotYrs), 
                          y = max(retroPhi$upr, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Synchrony")
agCVPlot <- ggplot(retroCVa, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.25) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "aggCV"),
            mapping = aes(x = min(plotYrs), 
                          y = max(retroCVa$upr, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Aggregate Variability")

png(here("figs/Fig1New_Retro.png"), height = 4.5, width = 6.5,
    units = "in", res = 300)
ggarrange(rawProdPlot, aggRetPlot, aggCatchPlot, compCVPlot, synchPlot, 
          agCVPlot, nrow = 2, ncol = 3, heights = c(1, 1.1))
dev.off()


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
  mutate(data = "prod", ts = "long") 
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
  mutate(data = "resid", ts = "long")
rDat <- data.frame(year = yrs,
                   wtdCV = rollapplyr(recMat, width = 12, 
                                      function(x) wtdCV(x),
                                      fill = NA, by.column = FALSE),
                   synch = rollapplyr(recMat, width = 12, 
                                      function(x) community.sync(x)$obs, 
                                      fill = NA, by.column = FALSE)) %>%
  mutate(agCV = sqrt(synch) * wtdCV) %>% 
  gather(key = var, value = index, 2:4) %>% 
  mutate(data = "rec", ts = "long") %>% 
  filter(!is.na(index))

# Compare canned estimates to Sean's model's output
cvC1 <- rDat %>% 
  filter(var == "wtdCV")
cvC2 <- modOut %>% 
  filter(term == "log_cv_s")
plot(cvC1$index ~ cvC2$est)

cvA1 <- rDat %>% 
  filter(var == "agCV")
cvA2 <- modOut %>% 
  filter(term == "log_cv_c")
plot(cvA1$index ~ cvA2$est)

phi1 <- rDat %>% 
  filter(var == "synch")
phi2 <- modOut %>% 
  filter(term == "logit_phi")
plot(phi1$index ~ phi2$est)
## Match up exactly
