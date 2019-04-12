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
  # filter(stk %in% longStks$stk) %>% 
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
  
colPal <- viridis(nrow(subCU), begin = 0, end = 1)
names(colPal) <- unique(plotDat$stkName)

png(file = paste(here(),"/figs/presentationFigs/spwn_lines.png", sep = ""),
    height = 3, width = 6, units = "in", res = 300)
ggplot(plotDat, aes(x = yr, y = S, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, name = "Stock") +
  labs(x = "Year", y = "Spawner Abundance") +
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
    datList[["recDev"]][ , , trialID] %>% 
      reshape2::melt() %>% 
      dplyr::rename("yr" = "Var1", "cu" =  "Var2", "recDev" = "value") %>% 
      mutate(sigma = as.factor(sigNames[i]), 
             synch = as.factor(synchNames[h]))
  })
  do.call(rbind, outList1)
})
outDat <- do.call(rbind, outList2) 
outDat$synch <- forcats::fct_relevel(outDat$synch, "high", after = Inf)

## Plot variance treatments first
varDat <- outDat %>% 
  mutate(stkName = as.factor(plyr::mapvalues(outDat$cu, 
                                             from = unique(outDat$cu),
                                             to = stkNames))) %>%
  filter(yr > 60,
         stkName %in% subCU$abbStkName, 
         synch == "med") %>% 
  mutate(stkName = factor(stkName),
         sigma = recode(factor(sigma), "low" = "Low Var.",
                        "high" = "High Var."),
         stkName = as.factor(plyr::mapvalues(stkName,
                                             from = unique(stkName),
                                             to = subCU$stkName)))

png(file = paste(here(),"/figs/presentationFigs/spwn_lines.png", sep = ""),
    height = 3, width = 6, units = "in", res = 300)
ggplot(varDat, aes(x = yr, y = recDev, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, name = "Stock") +
  labs(x = "Year", y = "Recruitment Deviation") +
  theme_sleekX() + 
  facet_wrap(~sigma)
dev.off()


colPal <- viridis(length(subCU), begin = 0, end = 1)
names(colPal) <- unique(standDat$stkName)



#For plotting purposes constrain to smallest populations 
plotDat <- outDat %>% 
  mutate(stkName = as.factor(plyr::mapvalues(outDat$cu, 
                                             from = unique(outDat$cu),
                                             to = stkNames))) %>%
  filter(yr > 60,
         stkName %in% subCU$abbStkName) %>% 
  mutate(stkName = factor(stkName),
         sigma = recode(factor(sigma), "low" = "Low Var.",
                       "high" = "High Var."),
         synch = recode(factor(synch), "high" = "High Synch.", 
                        "low" = "Low Synch.", "med" = "Med. Synch."),
         stkName = as.factor(plyr::mapvalues(stkName, 
                                             from = unique(stkName),
                                             to = subCU$stkName)))





# png(file = paste(here(),"/figs/april2019Meeting/recDev_lines.png", sep = ""),
#     height = 4, width = 6, units = "in", res = 300)
ggplot(standDat, aes(x = yr, y = recDev, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal) +
  labs(x = "Year", y = "Recruitment Deviations") +
  theme_sleekX() + 
  facet_wrap(~scen)
# dev.off()