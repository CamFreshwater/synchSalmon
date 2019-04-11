#*******************************************************************************
# pbsSeminarFigures.R
# Date revised: April 10, 2019
# Explainer: Tweaks previously generated plots for use in April seminar at
# Pacific Biological Station.
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

calcRollMean <- function(datMat, w = 4) {
  zoo::rollapplyr(datMat, width = w, FUN = mean, by = 1, fill = NA)
}

dd <- sapply(tt, function(x) calcRollMean(x, w = 4))

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
trialID <- sample.int(unique(outDat$trial), size = 1)

outList2 <- lapply(seq_along(synchDirNames), function(i) {
  datList <- readRDS(paste(here("outputs/simData"), synchDirNames[i], 
                           arrayNames[3], 
                           sep = "/"))
  datList[["recBY"]][ , , trialID] %>% 
    # calcRollMean() %>% 
    reshape2::melt() %>% 
    dplyr::rename("yr" = "Var1", "cu" =  "Var2", "rec" = "value") %>% 
    mutate(prod = prodNames[i]) 
})
outDat <- do.call(rbind, outList2) 

subCU <- cuPar %>% 
  select(stk, stkName, medianRec) %>% 
  filter(medianRec < 0.05)

  sample.int(subCU$stk, size = 5)
plotDat <- outDat %>% 
  mutate(stkName = as.factor(plyr::mapvalues(outDat$cu, 
                                             from = unique(outDat$cu),
                                             to = stkNames))) %>%
  filter(yr > 60,
         cu %in% subCU) %>% 
  mutate(stkName = factor(stkName))
colPal <- viridis(length(subCU), begin = 0, end = 1)
names(colPal) <- unique(plotDat$stkName)


# png(file = paste(here(),"/figs/april2019Meeting/recDev_lines.png", sep = ""),
#     height = 4, width = 6, units = "in", res = 300)
ggplot(plotDat, aes(x = yr, y = smoothRec, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal) +
  labs(x = "Year", y = "Recruitment") +
  theme_sleekX() + 
  facet_wrap(~prod)
# dev.off()