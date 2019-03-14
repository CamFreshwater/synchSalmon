#*******************************************************************************
# productivityTimeSeries.R
# Date revised: March 13, 2019
# Explainer: Gnerates time series of recruitment deviations figure. Originally
# intended for management workshop but discarded.
#*******************************************************************************

require(here)

## Load in data from synchrony analyses to visualize impacts of greater 
#variance/covariance
synchDirNames <- c("lowSig_sockeye", "highSig_sockeye")
stkNames <- genOutputList(synchDirNames[1], 
                          agg = FALSE)[["medSynch_TAM"]][["stkName"]]
nCUs <- length(stkNames)
#matrix of array names to be passed
arrayNames <- sapply(synchDirNames[1], function(x) { 
  list.files(paste(here("outputs/simData"), x, sep="/"), 
             pattern = "\\Arrays.RData$")
})
synchNames <- c("high", "low", "med")
sigNames <- c("low", "high")

# Function to pull and clean array data
outList2 <- lapply(seq_along(synchDirNames), function(i) {
  outList1 <-  lapply(seq_along(arrayNames), function(h) {
    datList <- readRDS(paste(here("outputs/simData"), synchDirNames[i], 
                             arrayNames[h], 
                             sep = "/"))
    datList[["recDev"]] %>% 
      reshape2::melt() %>% 
      dplyr::rename("yr" = "Var1", "cu" =  "Var2", "trial" = "Var3", 
                    "recDev" = "value") %>% 
      mutate(sigma = sigNames[i], synch = synchNames[h])
  })
  do.call(rbind, outList1)
})
outDat <- do.call(rbind, outList2) %>% 
  mutate(scen = paste(sigma, "Sig", synch, "Synch", sep = ""),
         stkName = as.factor(plyr::mapvalues(outDat$cu, 
                                             from = unique(outDat$cu),
                                             to = stkNames)))
# standardize abundance within a CU and trim
subCU <- sample.int(nCUs, size = 6)
trialID <- sample.int(unique(outDat$trial), size = 1)
standDat <- outDat %>% 
  filter(yr > 60,
         cu %in% subCU,
         trial == trialID,
         !scen %in% c("lowSigmedSynch", "highSigmedSynch")) %>% 
  # group_by(stkName, trial, scen) %>% 
  # mutate(interTrialMean = mean(recBY, na.rm = TRUE),
  #        interTrialSD = sd(recBY, na.rm = TRUE)) %>% 
  # mutate(standRecBY = (recBY - interTrialMean) / interTrialSD) %>% 
  # ungroup() %>% 
  mutate(stkName = factor(stkName), 
         scen = recode(factor(scen), "lowSiglowSynch" = "Low Var. - Low Synch.",
                       "highSiglowSynch" = "High Var. - Low Synch.",
                       "highSighighSynch" = "High Var. - High Synch.",
                       "lowSighighSynch" = "Low Var. - High Synch.")) %>% 
  mutate(scen = factor(scen, levels(scen)[c(4, 2, 3, 1)]))

colPal <- viridis(length(subCU), begin = 0, end = 1)
names(colPal) <- unique(standDat$stkName)


# png(file = paste(here(),"/figs/april2019Meeting/recDev_lines.png", sep = ""),
#     height = 4, width = 6, units = "in", res = 300)
ggplot(standDat, aes(x = yr, y = recDev, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal) +
  labs(x = "Year", y = "Recruitment Deviations") +
  theme_sleekX() + 
  facet_wrap(~scen)
# dev.off()