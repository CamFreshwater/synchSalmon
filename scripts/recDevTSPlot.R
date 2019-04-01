#*******************************************************************************
# presTimeSeries.R
# Date revised: April 1, 2019
# Explainer: Gnerates time series of recruitment deviations and of catch
# thresholds being met. 
#*******************************************************************************

require(here); require(samSim); require(tidyverse)

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
sigNames <- c("low", "high", "low", "high")

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

#### Equivalent analysis to above but focusing on years where aggregate return
### abundance above threshold
synchDirNames <- c("lowSig_sockeye", "highSig_sockeye", "lowSigLowA_sockeye",
                   "highSigLowA_sockeye")

#matrix of array names to be passed
aggTSNames <- sapply(synchDirNames[1], function(x) { 
  list.files(paste(here("outputs/simData"), x, sep="/"), 
             pattern = "\\TimeSeries.RData$")
})
synchNames <- c("high", "low", "med")
sigNames <- c("low", "high", "low", "high")
prodNames <- c("ref", "ref", "low", "low")

outList2 <- lapply(seq_along(synchDirNames), function(i) {
  outList1 <-  lapply(seq_along(aggTSNames), function(h) {
    datList <- readRDS(paste(here("outputs/simData"), synchDirNames[i], 
                             aggTSNames[h], 
                             sep = "/"))
    datList[["Ag Catch"]] %>% 
      reshape2::melt() %>% 
      dplyr::rename("yr" = "Var1", "trial" =  "Var2", "pm" = "value") %>% 
      mutate(prod = prodNames[i], sigma = sigNames[i], synch = synchNames[h])
  })
  do.call(rbind, outList1)
})
outAggDat <- do.call(rbind, outList2) %>% 
  mutate(scen = paste(prod, "Prod", sigma, "Sig", synch, "Synch", sep = ""))
# standardize abundance within a CU and trim
# trialID <- sample(unique(as.numeric(outAggDat$trial)), size = 1)
plotDat <- outAggDat %>% 
  filter(yr > 60,
         trial == "586",
         prod == "low",
         !scen %in% c("lowProdlowSigmedSynch", "lowProdhighSigmedSynch")
         ) %>% 
  mutate(scen = recode(factor(scen), "lowProdlowSiglowSynch" = "Low Var. - Low Synch.",
                       "lowProdhighSiglowSynch" = "High Var. - Low Synch.",
                       "lowProdhighSighighSynch" = "High Var. - High Synch.",
                       "lowProdlowSighighSynch" = "Low Var. - High Synch.")) %>%
  mutate(scen = factor(scen, levels(scen)[c(4, 2, 3, 1)]),
         col = as.factor(case_when(
           pm <= 1 ~ "red",
           pm > 1 ~ "green"
         )))

colPal <- c("red", "green")
names(colPal) <- levels(plotDat$col)
ggplot(plotDat, aes(x = yr, y = pm, colour = pm)) +
  geom_line() +
  # scale_color_manual(values = colPal, guide = FALSE) +
  scale_colour_gradientn(colours = c("red", "yellow", "green"), 
                         values = scales::rescale(c(0.01, 1, 2)),
                         limits = c(0.01, max(plotDat$pm))) +
  geom_hline(yintercept = 1) +
  theme_sleekX() +
  facet_wrap(~scen)
