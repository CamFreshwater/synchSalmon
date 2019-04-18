#*******************************************************************************
# pbsSeminarFigures.R
# Date revised: April 10, 2019
# Explainer: Tweaks previously generated plots for use in April seminar at
# Pacific Biological Station.
#*******************************************************************************


require(here); require(synchrony); require(zoo); require(ggplot2); 
require(tidyverse); require(viridis); require(ggpubr); require(samSim)


#------------

## Productivity trends (observed)

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

rawDat <- recDatTrim1 %>%
  select(stk, yr, rollProd, prodZ, logProd) %>% 
  mutate(stk = as.factor(stk))

## Trends in productivity
stkN <- length(unique(rawDat$stk))
plotYrs <- c(min(rawDat$yr), max(rawDat$yr))

png(file = here("figs", "presentationFigs", "prodTrends.png"),
    height = 5, width = 6, units = "in", res = 450)
ggplot(rawDat, aes(x = yr, y = rollProd, col = stk)) + 
  labs(x = "", y = "Productivity (log(Recruits/Spawner))") + 
  ylim(min(rawDat$rollProd),  max(rawDat$rollProd)) +
  geom_line(size = 0.55) +
  xlim(plotYrs) +
  scale_color_manual(values = rep("grey", length.out = stkN), guide = FALSE) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", size = 1.35)  +
  theme_sleekX(axisSize = 13) +
  theme(axis.text.y = element_text(hjust = 0.5))
dev.off()


### Make three panel plot with productivity, comp cv and synch

# Trim SR data
longStks <- recDatTrim1 %>% 
  group_by(stk) %>% 
  summarize(startDate = min(yr)) %>% 
  filter(!(startDate > 1951 | stk == "11")) %>% 
  select(stk) %>% 
  as.vector()
recDatTrim <- recDatTrim1 %>% 
  dplyr::filter(stk %in% longStks$stk)
yrs <- unique(recDatTrim$yr)

rawDat <- recDatTrim1 %>%
  select(stk, yr, rollProd) %>% 
  filter(stk %in% longStks$stk) %>% 
  mutate(stk = as.factor(stk))

# Import TMB model outputs
modOut <- readRDS(here("outputs", "generatedData", "tmbSynchEst.rds"))
modOut$broodYr <- rep(yrs[12:61], times = 3)
retroCVc <- modOut %>% 
  filter(term == "log_cv_s")  
retroPhi <- modOut %>% 
  filter(term == "logit_phi")

# Make fig
stkN <- length(unique(recDatTrim1$stk))

rawProdPlot <- ggplot(rawDat, aes(x = yr, y = rollProd, col = stk)) + 
  labs(x = "", y = "log(Recruits/Spawner)") + 
  ylim(min(rawDat$rollProd),  max(rawDat$rollProd)) +
  geom_line(size = 0.45) +
  xlim(plotYrs) +
  scale_color_manual(values = rep("grey", length.out = stkN), guide = FALSE) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", size = 1.25)  +
  theme_sleekX(axisSize = 10, position = "bottom") +
  theme(axis.text.y = element_text(hjust = 0.5))
compCVPlot <- ggplot(retroCVc, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.15) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 10, position = "bottom") +
  theme(axis.text.y = element_text(hjust = 0.5)) +
  labs(x = "", y = "Component Variability") 
synchPlot <- ggplot(retroPhi, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.25) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 10) +
  theme(axis.text.y = element_text(hjust = 0.5)) +
  labs(x = "", y = "Synchrony")

png(here("figs/presentationFigs/threePanelRetro.png"), height = 2.25, width = 7,
    units = "in", res = 450)
# png(here("figs/presentationFigs/threePanelRetro.png"), width = 2.25, height = 7,
#     units = "in", res = 450)
ggarrange(rawProdPlot, compCVPlot, synchPlot, ncol = 3)
dev.off()

#------------


### Simulated population dynamics plots

### Load in data from synchrony analyses to visualize impacts of lower 
# productivity
cuPar <- read.csv(here("data/sox/fraserCUpars.csv"), stringsAsFactors = F)

synchDirNames <- c("medSig_sockeye", "medSigLowA_sockeye")
stkNames <- genOutputList(synchDirNames[1], 
                          agg = FALSE)[["medSynch_TAM"]][["stkName"]]
nCUs <- length(stkNames)
#matrix of array names to be passed
arrayNames <- sapply(synchDirNames[1], function(x) { 
  list.files(paste(here("outputs/simData"), x, sep="/"), 
             pattern = "\\Arrays.RData$")
})
prodNames <- c("ref", "low")

# Function to pull and clean array data
outList2 <- lapply(seq_along(synchDirNames), function(i) {
  datList <- readRDS(paste(here("outputs/simData"), synchDirNames[i], 
                           arrayNames[3], 
                           sep = "/"))
  set.seed(222)
  trialID <- sample.int(dim(datList[["S"]])[3], size = 1)
  datList[["S"]][ , , trialID] %>% 
    reshape2::melt() %>% 
    dplyr::rename("yr" = "Var1", "cu" =  "Var2", "S" = "value") %>% 
    mutate(prod = as.factor(prodNames[i]))
})
outDat <- do.call(rbind, outList2) 

#For plotting purposes constrain to smallest populations 
subCU <- cuPar %>% 
  filter(medianRec < 0.05, 
         !stkName == "Harrison") %>% 
  select(stkName) %>%
  mutate(abbStkName = abbreviate(stkName, minlength = 4)) %>% 
  dplyr::sample_n(., size = 5)

plotDat <- outDat %>% 
  mutate(stkName = as.factor(plyr::mapvalues(outDat$cu, 
                                             from = unique(outDat$cu),
                                             to = stkNames))) %>%
  filter(yr > 60,
         stkName %in% subCU$abbStkName) %>% 
  mutate(stkName = factor(stkName),
         prod = recode(factor(prod), "low" = "Low Prod.",
                       "ref" = "Reference Prod."),
         stkName = as.factor(plyr::mapvalues(stkName, 
                                             from = unique(stkName),
                                             to = subCU$stkName)))
  
# colPal <- viridis(nrow(subCU), begin = 0, end = 1)
colPal <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
names(colPal) <- unique(plotDat$stkName)

png(file = paste(here(),"/figs/presentationFigs/spwn_lines.png", sep = ""),
    height = 3, width = 6, units = "in", res = 300)
ggplot(plotDat, aes(x = yr, y = S, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, name = "Stock") +
  labs(x = "Year", y = "Spawner Abundance (millions)") +
  theme_sleekX() + 
  facet_wrap(~prod)
dev.off()

## Load in data from synchrony analyses to visualize impacts of greater 
#variance/covariance
synchDirNames <- c("lowSig_sockeye", "highSig_sockeye")
stkNames <- genOutputList(synchDirNames[1], 
                          agg = FALSE)[["medSynch_TAM"]][["stkName"]]
#matrix of array names to be passed
arrayNames <- sapply(synchDirNames[1], function(x) { 
  list.files(paste(here("outputs/simData"), x, sep="/"), 
             pattern = "\\Arrays.RData$")
})
synchNames <- c("high", "low", "med")
sigNames <- c("low", "high", "low", "high")

# Function to pull and clean array data
outList2 <- lapply(seq_along(synchDirNames), function(i) {
  outList1 <-  lapply(seq_along(arrayNames), function(h) {
    datList <- readRDS(paste(here("outputs/simData"), synchDirNames[i], 
                             arrayNames[h], 
                             sep = "/"))
    recDev <- datList[["recDev"]][ , , trialID] %>% 
      reshape2::melt() %>% 
      dplyr::rename("yr" = "Var1", "cu" =  "Var2", "recDev" = "value") %>% 
      mutate(sigma = as.factor(sigNames[i]), 
             synch = as.factor(synchNames[h]))
    spwn <- datList[["S"]][ , , trialID] %>% 
      reshape2::melt() %>% 
      dplyr::rename("yr" = "Var1", "cu" =  "Var2", "S" = "value") %>% 
      mutate(prod = as.factor(prodNames[i]))
    rec <- datList[["recBY"]][ , , trialID] %>% 
      reshape2::melt() %>% 
      dplyr::rename("yr" = "Var1", "cu" =  "Var2", "R" = "value") %>% 
      mutate(prod = as.factor(prodNames[i]))
    recDev %>% 
      mutate(S = spwn$S, 
             R = rec$R)
  })
  do.call(rbind, outList1)
})
outDat <- do.call(rbind, outList2) 
outDat <- outDat %>% 
  mutate(stkName = as.factor(plyr::mapvalues(outDat$cu, 
                                             from = unique(outDat$cu),
                                             to = stkNames)))
outDat$synch <- forcats::fct_relevel(outDat$synch, "high", after = Inf) 

## Plot example SR curve
srDat <- outDat %>% 
  filter(yr > 60,
         stkName == "UppP", 
         synch == "med") %>% 
  mutate(sigma = recode(factor(sigma), "low" = "Low Var.",
                        "high" = "High Var."))
pars <- cuPar %>% 
  filter(stkName == "Upper Pitt") %>% 
  select(alpha, beta0)


png(file = paste(here(),"/figs/presentationFigs/srCurveEmpty.png", sep = ""),
    height = 3, width = 4, units = "in", res = 300)
ggplot(srDat, aes(x = S, y = R)) +
  geom_point(colour = "white") +
  labs(x = "Spawners", y = "Recruits") +
  theme_sleekX() +
  stat_function(fun = function(x) x * exp(pars$alpha - pars$beta0 * x)) 
dev.off()

png(file = paste(here(),"/figs/presentationFigs/srCurve.png", sep = ""),
    height = 3, width = 6, units = "in", res = 300)
ggplot(srDat, aes(x = S, y = R)) +
  geom_point() +
  labs(x = "Spawners", y = "Recruits") +
  theme_sleekX() +
  stat_function(fun = function(x) x * exp(pars$alpha - pars$beta0 * x)) +
  facet_wrap(~sigma)
dev.off()

## Plot variance treatments first
varDat <- outDat %>% 
  filter(yr > 60,
         stkName %in% subCU$abbStkName, 
         synch == "med") %>% 
  mutate(stkName = factor(stkName),
         sigma = recode(factor(sigma), "low" = "Low Var.",
                        "high" = "High Var."),
         stkName = as.factor(plyr::mapvalues(stkName,
                                             from = unique(stkName),
                                             to = subCU$stkName)))

png(file = paste(here(),"/figs/presentationFigs/prod_var.png", sep = ""),
    height = 3, width = 6, units = "in", res = 300)
ggplot(varDat, aes(x = yr, y = recDev, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, name = "Stock") +
  labs(x = "Year", y = "Recruitment Deviation") +
  theme_sleekX() + 
  facet_wrap(~sigma)
dev.off()

## Plot synchrony treatments second
synchDat <- outDat %>% 
  filter(yr > 60,
         stkName %in% subCU$abbStkName, 
         sigma == "low",
         !synch == "med") %>% 
  mutate(stkName = factor(stkName),
         synch = recode(factor(synch), "low" = "Low Synch.",
                        "high" = "High Synch."),
         stkName = as.factor(plyr::mapvalues(stkName,
                                             from = unique(stkName),
                                             to = subCU$stkName)))

png(file = paste(here(),"/figs/presentationFigs/prod_synch.png", sep = ""),
    height = 3, width = 6, units = "in", res = 300)
ggplot(synchDat, aes(x = yr, y = recDev, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, name = "Stock") +
  labs(x = "Year", y = "Recruitment Deviation") +
  theme_sleekX() + 
  facet_wrap(~synch)
dev.off()


## Empty time series
dum <-  outDat %>% 
  filter(yr > 60,
         cu == "4",
         sigma == "low",
         synch == "med")
png(file = paste(here(),"/figs/presentationFigs/blankTS.png", sep = ""),
    height = 1.5, width = 3, units = "in", res = 300)
ggplot(dum, aes(x = yr, y = recDev)) +
  geom_line(size = 1.5) +
  labs(x = "Year", y = "Recruitment Deviation") +
  theme_sleekX()
dev.off()


#------------


### Synchrony simulation analysis
# Plots of a) each PM at reference and b) juxtaposed w/ low

## Load data
plotDat <- read.csv(here("outputs/generatedData", "summaryTable_noSkew.csv"), 
                    stringsAsFactors = F) %>% 
  select(-X) %>% 
  filter(!om == "Low Prod. Heavy Tails",
         var %in% c("medRecRY", "ppnMixedOpen", "ppnYrsHighCatch", 
                    "stabilityCatch")) %>% 
  mutate(synch = as.factor(synch),
         om = as.factor(om),
         sigma = as.factor(sigma)) %>% 
  mutate(synch = fct_relevel(synch, "high", after = 2),
         sigma = fct_relevel(sigma, "high", after = 2),
         om = fct_relevel(om, "Low Prod.", after = 1))

colPal <- viridis(length(levels(plotDat$synch)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$synch)

# Make individual plots
pmPlot <- function(temp, catchYLab, facetCol = 3, pos = NULL, dotSize = 3.25, 
                   lineSize = 0.8, legSize = 0.7, axSize = 10, 
                   facetSize = 0.95) {
  q <- ggplot(temp, aes(x = sigma, y = medn, ymin = lowQ, ymax = highQ,
                   fill = synch)) +
    labs(x = "Variance", y = catchYLab, 
         fill = "Sim.\nParameter\nValue") +
    geom_pointrange(shape = 21, fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.65)) +
    scale_x_discrete(labels = c("low" = "Low",
                                "med" = "Reference",
                                "high" = "High")) +
    scale_fill_manual(name = "Synchrony", values = colPal,
                      labels = c("low" = "Low",
                                 "med" = "Moderate",
                                 "high" = "High")) +
    theme_sleekX(legendSize = legSize, axisSize = axSize, facetSize = facetSize)
  
  if (length(unique(temp$om)) > 1) {
    q <- q +
      facet_wrap(~om, scales = "fixed", ncol = facetCol, nrow = 1)
  }
  if (is.null(pos) == FALSE) {
    q <- q + theme_sleekX(position = pos, axisSize = axSize, 
                          facetSize = facetSize)
  }
  return(q)
}

retDat <- plotDat %>% 
  filter(var == "medRecRY")
png(file = paste(here(),"/figs/presentationFigs/retPMRef.png", sep = ""),
    height = 2.5, width = 5, units = "in", res = 300)
pmPlot(retDat %>% filter(om == "Reference Prod."),
       catchYLab = "Agg. Return\nAbundance", dotSize = 4.5)
dev.off()

png(file = paste(here(),"/figs/presentationFigs/retPM.png", sep = ""),
    height = 2.5, width = 6, units = "in", res = 300)
retPlot <- pmPlot(retDat, catchYLab = "Agg. Return\nAbundance", dotSize = 4.5,
                  facetSize = 1.1)
print(retPlot)
dev.off()

escDat <- plotDat %>% 
  filter(var == "ppnMixedOpen")
png(file = paste(here(),"/figs/presentationFigs/escPMRef.png", sep = ""),
    height = 2.5, width = 5, units = "in", res = 300)
pmPlot(escDat %>% filter(om == "Reference Prod."),
       catchYLab = "Ppn. of Esc.\nGoals Met",
       dotSize = 4.5)
dev.off()

png(file = paste(here(),"/figs/presentationFigs/escPM.png", sep = ""),
    height = 2.5, width = 6, units = "in", res = 300)
escPlot <- pmPlot(escDat, catchYLab = "Ppn. of Esc.\nGoals Met",
                  dotSize = 4.5, facetSize = 1.1)
print(escPlot)
dev.off()

stbDat <- plotDat %>% 
  filter(var == "stabilityCatch")
png(file = paste(here(),"/figs/presentationFigs/stbPMRef.png", sep = ""),
    height = 2.5, width = 5, units = "in", res = 300)
pmPlot(stbDat %>% filter(om == "Reference Prod."),
       catchYLab = "Catch Stability",
       dotSize = 4.5)
dev.off()

png(file = paste(here(),"/figs/presentationFigs/stbPM.png", sep = ""),
    height = 2.5, width = 6, units = "in", res = 300)
stbPlot <- pmPlot(stbDat, catchYLab = "Catch Stability", dotSize = 4.5, 
                  facetSize = 1.1)
print(stbPlot)
dev.off()

retPlot <- pmPlot(retDat, catchYLab = "Agg. Return\nAbundance", dotSize = 4.5,
                  facetSize = 1.2, axSize = 13, pos = "top")
escPlot <- pmPlot(escDat, catchYLab = "Ppn. of Years\nEsc. Goal Met",
                  dotSize = 4.5, facetSize = 1.2, axSize = 13, pos = "mid")
stbPlot <- pmPlot(stbDat, catchYLab = "Catch Stability", dotSize = 4.5, 
                  facetSize = 1.2, axSize = 13, pos = "bottom")

png(file = paste(here(),"/figs/presentationFigs/allPMs.png", sep = ""),
    height = 7, width = 6, units = "in", res = 300)
ggarrange(retPlot, escPlot, stbPlot,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1.2))
dev.off()
