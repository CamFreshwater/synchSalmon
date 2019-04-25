#*******************************************************************************
# runModel_mainSim.R
# Date revised: ONGOING
# Inputs: recoverySim.R
# Outputs: pdf plots
# Explainer: Runs closed loop simulation model with two different productivity 
# OMs, three component variability and three synchrony OMs; constant MP (TAM
# rule)
# Analysis for: Weakened portfolio effects constrain management effectiveness 
# for population aggregates
#*******************************************************************************

# Check if required packages are installed and run
listOfPackages <- c("plyr", "here", "parallel", "doParallel", "foreach", 
                    "reshape2", "tidyverse", "gsl", "tictoc", "stringr", 
                    "synchrony", "zoo", "Rcpp", "RcppArmadillo", "sn", 
                    "sensitivity", "mvtnorm", "forcats", "ggpubr", "viridis", 
                    "samSim")

here <- here::here

newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if(length(newPackages)) install.packages(newPackages)
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data/sox/fraserOMInputs_varyCorr.csv"), 
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/sox/fraserCUpars.csv"), stringsAsFactors = F)
srDat <- read.csv(here("data/sox/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/sox/fraserCatchDatTrim.csv"), 
                     stringsAsFactors = F)
ricPars <- read.csv(here("data/sox/pooledRickerMCMCPars.csv"), 
                    stringsAsFactors = F)
larkPars <- read.csv(here("data/sox/pooledLarkinMCMCPars.csv"), 
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/sox/tamRefPts.csv"), stringsAsFactors = F)


### SET UP MODEL RUN -----------------------------------------------------

## Define simulations to be run
nTrials <- 1500

simParTrim <- subset(simPar, scenario %in% c("lowSig", "medSig", "highSig",
                                             "lowSigLowA", "medSigLowA",
                                             "highSigLowA",
                                             "lowSigLowStudT", "medSigLowStudT",
                                             "highSigLowStudT"
                                             ))

scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species),
                                                sep = "_"))

# recoverySim(simParTrim[24, ], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             dirName = "test", nTrials = 5, makeSubDirs = FALSE, random = FALSE)
# for(i in seq_along(dirNames)) {
# d <- subset(simParTrim, scenario == scenNames[i])
# simsToRun <- split(d, seq(nrow(d)))
# lapply(simsToRun, function(x) {
#   recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
#               ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirNames[i],
#               nTrials=125, makeSubDirs=TRUE)
# })
# }
 
for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 1) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MASS),
                     library(here),
                     library(sensitivity),
                     library(mvtnorm),
                     library(scales), #shaded colors for figs
                     library(viridis), #color blind gradient palette
                     library(gsl),
                     library(dplyr),
                     library(Rcpp),
                     library(RcppArmadillo),
                     library(sn),
                     library(samSim)))
  #export custom function and objects
  clusterExport(cl, c("simsToRun", "recoverySim", "cuPar", "dirName", "nTrials",
                      "catchDat", "srDat", "ricPars", "dirName", "larkPars",
                      "tamFRP"), envir = environment())
  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    recoverySim(x, cuPar, catchDat = catchDat, srDat = srDat,
                variableCU = FALSE, ricPars, larkPars = larkPars,
                tamFRP = tamFRP, cuCustomCorrMat = NULL, dirName = dirName,
                nTrials = nTrials, makeSubDirs = FALSE, random = FALSE)
  })
  stopCluster(cl) #end cluster
  toc()
}


#_________________________________________________________________________
vars <- c("medRecRY", "ppnCUUpper", "ppnMixedOpen",
          "medCatch", "ppnYrsHighCatch", "stabilityCatch")
omNames <- rep(c("ref", "lowA", "lowStudT"), each = 3)
sigNames <- rep(c("low", "med", "high"), length.out = length(omNames))

plotDat1 = NULL
for(h in seq_along(dirNames)) {
  agList <- genOutputList(dirNames[h], agg = TRUE)
  keyVar <- sapply(agList, function(x) unique(x$keyVar))
  plotOrder <- sapply(agList, function(x) unique(x$plotOrder))
  singleScen = NULL
  for (i in seq_along(vars)) {
    dum <- data.frame(sigma = rep(sigNames[h], length.out = length(agList)),
                      om = rep(omNames[h], length.out = length(agList)),
                      var = rep(vars[i], length.out = length(agList)),
                      synch = as.factor(keyVar),
                      medn = sapply(agList, function(x) median(x[,vars[i]])),
                      lowQ = sapply(agList, function(x) qLow(x[,vars[i]])),
                      highQ = sapply(agList, function(x) qHigh(x[,vars[i]])),
                      row.names = NULL
    )
    singleScen <- rbind(singleScen, dum)
  }
  rownames(singleScen) <- c()
  #merge multiple scenarios into one dataframe
  plotDat1 <- rbind(plotDat1, singleScen) 
}
plotDat1 <- plotDat1 %>%
  mutate(synch = recode(synch, "0.05" = "low", "0.5" = "med", "0.75" = "high", 
                     .default = levels(synch)))



### Generate standardized CU-specific return abundances; load in CU-specific 
# time series data and scale each trial's results relative to low CVc, low 
# synch scenario then average across CUs for each trial
arrayNames <- sapply(dirNames, function(x) {
 list.files(paste(here("outputs/simData"), x, sep="/"),
            pattern = "\\Arrays.RData$")
})

#iterate over productivity scenarios
listByOM <- lapply(seq_along(dirNames), function (h) { 
 #iterate over synch scenarios
 listBySynch <- lapply(seq_along(arrayNames[, h]), function (x) {
   datList <- readRDS(paste(here("outputs/simData"), dirNames[h],
                            arrayNames[x, h], sep = "/"))
   datList$recRY %>% 
     reshape2::melt() %>% 
     dplyr::rename("yr" = "Var1", "cu" =  "Var2", "trial" = "Var3", 
                   "recRY" = "value") %>% 
     mutate(sigma = as.factor(sigNames[h]), synch = as.factor(datList$nameOM), 
            om = as.factor(omNames[h]),
            scenID = as.factor(paste(omNames[h], sigNames[h], datList$nameOM, 
                                     sep = "_"))) %>% 
     filter(yr > datList$nPrime)
 }) #iterate across different OMs within a scenario
 dat <- do.call(rbind, listBySynch) #merge synch list elements into DF
})
fullDat <- do.call(rbind, listByOM)

# standardize data relative to low synch/low sigma in reference OM
#First calculate within CU medians for low synch/low sigma/reference dataset
lowAggV <- fullDat %>% 
  filter(scenID == "ref_low_lowSynch") %>% 
  dplyr::group_by(cu, trial, sigma, synch, om, scenID) %>% 
  dplyr::summarize(medR = median(recRY)) %>%
  dplyr::group_by(cu, sigma, synch, om, scenID) %>% 
  mutate(lowScenMedR = median(medR)) %>% 
  as.data.frame
#Trimmed version that can be merged and used to calculate relative differences
trimLowV <- lowAggV %>% 
  select(cu, trial, scenID, lowScenMedR)

scens <- unique(fullDat$scenID)
#standardize by subtracting median (lowScenMedR)
stdInnerList <- lapply(seq_along(scens), function (i) {
  if (scens[i] == "ref_low_lowSynch") {
    out <- lowAggV %>% 
      mutate(stdMedR = medR - lowScenMedR) %>% 
      select(-scenID, -lowScenMedR)
  } else {
    out <- fullDat %>% 
      filter(scenID == scens[i]) %>% 
      group_by(cu, trial, sigma, synch, om) %>% 
      #calculate median CU-specific recruitment within a trial
      dplyr::summarize(medR = median(recRY)) %>% 
      #join so rel. diff can be calc
      inner_join(., trimLowV, by = c("cu", "trial")) %>%
      group_by(cu) %>% 
      #calc rel. diff
      mutate(stdMedR = (medR - lowScenMedR)) %>% 
      #remove values pulled from lowAggV
      select(-lowScenMedR, -scenID) %>%
      as.data.frame
  }
  return(out)
})

#merge and calculate means across CUs per trial
temp <- do.call(rbind, stdInnerList) %>% 
  group_by(trial, sigma, synch, om) %>% 
  dplyr::summarize(meanStdRecRY = mean(stdMedR))

#finally calculate medians and quantiles to generate dataset equivalent to 
#plotdat
finalOut <- temp %>% 
  group_by(sigma, synch, om) %>% 
  mutate(var = "stdRecRY",
         mean = mean(meanStdRecRY),
         medn = median(meanStdRecRY),
         lowQ = qLow(meanStdRecRY),
         highQ = qHigh(meanStdRecRY)) %>%
  select(sigma, om, var, synch, medn, lowQ, highQ) 

stdFullDat <- finalOut %>% 
  as.data.frame() %>% 
  mutate(synch = recode(synch, "lowSynch" = "low", "medSynch" = "med",
                        "highSynch" = "high", .default = levels(synch)))
saveRDS(stdFullDat, here("outputs/generatedData/stdMedRecRY.rds"))

plotDat <- rbind(plotDat1, stdFullDat) %>%
  mutate(om = recode(om, "ref" = "Reference Prod.", "lowA" = "Low Prod.",
                     "lowStudT" = "Low Prod. Heavy Tails", 
                     .default = levels(om))) 

#Save summary data to pass to Rmd
write.csv(plotDat, here("outputs/generatedData", "summaryTable_3ProdOMs.csv"),
          row.names = FALSE)


# If available, skip above and load in summary data directly
plotDat <- read.csv(here("outputs/generatedData", "summaryTable_3ProdOMs.csv"),
                    stringsAsFactors = T)

## Subset plotDat to remove low productivity, heavy tailed scenario (presented
# in supplement only)
plotDat <- plotDat %>% 
  filter(!om == "Low Prod. Heavy Tails") %>% 
  mutate(synch = fct_relevel(synch, "high", after = Inf),
         sigma = fct_relevel(sigma, "high", after = Inf),
         om = fct_relevel(om, "Low Prod.", after = Inf))


### Make some figures 
# Default plotting variable settings
dotSize = 3.25; lineSize = 0.8; legSize = 0.8; axSize = 10; facetSize = 0.95

# Default color palette
colPal <- viridis(length(levels(plotDat$synch)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$synch)

consVars <- c("medRecRY", "stdRecRY", "ppnMixedOpen", "ppnCUUpper")
consYLabs <- c("Agg. Return\nAbundance", "CU-Specific Std.\nReturn Abundance", 
               "Prop. MUs\nAbove Esc. Goal", "Prop. CUs\nAbove Benchmark")
#make dataframe of labels to annotate facets
consLabs <- data.frame(om = rep(factor(unique(plotDat$om),
                                       levels = unique(plotDat$om)),
                                each = 4),
                       lab = c("a)", "c)", "e)", "g)", "b)", "d)", "f)", "h)"),
                       var = rep(factor(consVars, levels = unique(consVars)),
                                 times = 2)
)

consPlots <- lapply(seq_along(consVars), function(i) {
  temp <- plotDat %>%
    filter(var == consVars[i])
  q <- ggplot(temp, aes(x = sigma, y = medn, ymin = lowQ, ymax = highQ,
                        fill = synch)) +
    labs(x = "Component Variance", y = consYLabs[i], 
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
    geom_text(data = consLabs %>% filter(var == consVars[i]),
              mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab,
                            hjust = 1, vjust = 0),
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~om, scales = "fixed", ncol = 4, nrow = 1)
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 2 | i == 3) {
    q <- q + theme_sleekX(position = "mid", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 4) {
    q <- q + theme_sleekX(position = "bottom", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  return(q)
})

catchVars <- c("medCatch", "ppnYrsHighCatch", "stabilityCatch")
catchYLabs <- c("Catch\nAbundance", "Prop. Years Above\nCatch Threshold",
                "Catch\nStability")
catchLabs <- data.frame(om = rep(factor(unique(plotDat$om),
                                        levels = unique(plotDat$om)),
                                 each = 3),
                        lab = c("a)", "c)", "e)", "b)", "d)", "f)"),
                        var = rep(factor(catchVars, levels = unique(catchVars)),
                                  times = 2))

catchPlots <- lapply(seq_along(catchVars), function(i) {
  temp <- plotDat %>%
    filter(var == catchVars[i])
  q <- ggplot(temp, aes(x = sigma, y = medn, ymin = lowQ, ymax = highQ,
                        fill = synch)) +
    labs(x = "Component Variance", y = catchYLabs[i], 
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
    geom_text(data = catchLabs %>% filter(var == catchVars[i]),
              mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab,
                            hjust = 1, vjust = 0),
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~om, scales = "fixed", ncol = 4, nrow = 1)
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 2) {
    q <- q + theme_sleekX(position = "mid", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 3) {
    q <- q + theme_sleekX(position = "bottom", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  return(q)
})

png(file = paste(here(),"/figs/Fig3_consGroupedPlots.png", sep = ""),
    height = 5.5, width = 5.25, units = "in", res = 600)
ggarrange(consPlots[[1]], consPlots[[2]], consPlots[[3]], consPlots[[4]],
          ncol = 1, nrow = 4, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,0.9,0.9,1.15))
dev.off()
png(file = paste(here(),"/figs/Fig4_catchGroupedPlots.png", sep = ""),
    height = 4.5, width = 5.25, units = "in", res = 600)
ggarrange(catchPlots[[1]], catchPlots[[2]], catchPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,0.9,1.2))
dev.off()


#_________________________________________________________________________
# Generate box plots (originally time series) of component CV, synch and ag CV;
# Note that to make comparable to retrospective analysis only show 10 
# refDirNames <- dirNames[1:3]
# 
# plotList = vector("list", length = length(refDirNames))
# arrayNames <- sapply(refDirNames, function(x) {
#  list.files(paste(here("outputs/simData"), x, sep="/"),
#             pattern = "\\Arrays.RData$")
# })
# 
# tic("runParallel")
# Ncores <- detectCores()
# cl <- makeCluster(Ncores - 1) #save two cores
# registerDoParallel(cl)
# clusterEvalQ(cl, c(library(here), library(synchrony), library(zoo),
#                   library(parallel), library(doParallel), library(foreach),
#                   library(samSim)))
# newAgTSList <- lapply(seq_along(refDirNames), function (h) {
#  #export custom function and objects
#  clusterExport(cl, c("refDirNames", "arrayNames", "calcSynchMetrics", "wtdCV",
#                      "genOutputList", "h"), envir = environment())
#  listSynchLists <- parLapply(cl, 1:length(arrayNames[, h]), function(x) {
#    datList <- readRDS(paste(here("outputs/simData"), refDirNames[h],
#                             arrayNames[x, h], sep = "/"))
#    synchList <- calcSynchMetrics(datList, corr = FALSE,
#                                  weight = TRUE, windowSize = 12)
#    synchList <- c(datList$nameOM, synchList)
#    names(synchList)[1] <- "opMod"
#    return(synchList)
#  }) #iterate across different OMs within a scenario
#  # pull list of time series metrics estimated internally, then merge with
#  # synch list generated above based on common op model
#  agTSList <- genOutputList(refDirNames[h], agg = TRUE, aggTS = TRUE)
#  for(j in seq_along(agTSList)) {
#    for(k in seq_along(listSynchLists)) {
#      om <- agTSList[[j]]$opMod
#      if (listSynchLists[[k]]$opMod == om) {
#        agTSList[[j]] <- c(agTSList[[j]], listSynchLists[[k]][-1])
#      }
#    }
#  }
#  plotList[[h]] <- agTSList
# }) #iterate across different scenarios
# names(newAgTSList) <- refDirNames
# stopCluster(cl)
# toc()

## Save
# saveRDS(newAgTSList, here("outputs/generatedData/synchTS/synchTSList.rda"))
newAgTSList <- readRDS(here("outputs/generatedData/synchTS/synchTSList.rda"))

### Manipulate lists to create plottable data structure 
## In this case that is one median value per year per unique combo of OMs
refOmNames <- rep("ref", each = 3)
refSigNames <- c("lowSigma", "medSigma", "highSigma")
fullList <- sapply(seq_along(refDirNames), function(h) {
  d <- newAgTSList[[h]]
  nYears <- d[[1]]$`nYears`
  simLength <- d[[1]]$`nYears` - d[[1]]$`nPrime`
  firstYear <- d[[1]]$`firstYr`
  start <- d[[1]]$`nPrime` + firstYear
  prodNames <- refOmNames[h]
  #subset list so it contains only sigma, synch and vars of interest,
  # calculate medians, and combine
  trimList <- lapply(d, function(x) {
    dat1 <- data.frame(sigmaOM = rep(refSigNames[h], length.out = nYears),
                      synchOM = rep(x[["opMod"]], length.out = nYears),
                      prodOM = rep(prodNames, length.out = nYears),
                      year = seq(from = firstYear,
                                 to = (firstYear + nYears - 1))
                      ) %>%
      mutate(sigmaOM = as.character(sigmaOM),
             synchOM = as.character(synchOM),
             medSynchRecBY = apply(x[["synchRecBY"]], 1, median),
             medCompCVRecBY = apply(x[["compCVRecBY"]], 1, median)
      )
    dat1[dat1$year < start, c("sigmaOM", "synchOM")] <- "obs"
    dat2 <- dat1 %>%
      mutate(sigmaOM = factor(factor(sigmaOM),
                              levels = c("obs", "lowSigma", "medSigma",
                                         "highSigma")),
             synchOM = factor(factor(synchOM),
                              levels = c("obs", "lowSynch", "medSynch",
                                         "highSynch"))
      )
    return(dat2)
  })
})
plotDat <- do.call(rbind, fullList)

# saveRDS(plotDat, file = here("outputs/generatedData/synchTS/fullSynch_refOMs.rda"))
plotDat <- readRDS(file = here("outputs/generatedData/synchTS/fullSynch_refOMs.rda"))
start <- plotDat %>%
  filter(!sigmaOM == "obs") %>%
  summarise(min(year))
start <- start[[1]]

dum <- plotDat %>%
  filter(synchOM == "medSynch" | synchOM == "obs",
         prodOM == "ref",
         !is.na(medCompCVRecBY))
colPal <- c("black", viridis(length(unique(dum$sigmaOM)) - 1, begin = 0, 
                             end = 1))
names(colPal) <- levels(dum$sigmaOM)

## Faceted boxplot
q2 <- ggplot(dum, aes(x = sigmaOM, y = medCompCVRecBY, fill = sigmaOM)) +
  labs(x = "Component Variability Scenario", y = "Median CVc of Returns", 
       title = NULL) +
  geom_boxplot(size = 0.8, alpha = 0.5) +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(name = "Operating Model", values = colPal) +
  scale_x_discrete(breaks = waiver(), labels = c("Observed", 
                                                 expression(paste("0.75", 
                                                                  sigma["R"])),
                                                 expression(paste(sigma["R"])),
                                                 expression(paste("1.25", 
                                                                  sigma["R"]))
  )) +
  theme_sleekX(position = "standard", axisSize = axSize * 0.9)

dum2 <- plotDat %>%
  dplyr::filter(sigmaOM == "medSigma" | sigmaOM == "obs",
                prodOM == "ref",
                !is.na(medSynchRecBY))
colPal2 <- c("black", viridis(length(unique(dum2$synchOM)) - 1, begin = 0, 
                              end = 1))
names(colPal2) <- levels(dum2$synchOM)
p2 <- ggplot(dum2, aes(x = synchOM, y = medSynchRecBY, fill = synchOM)) +
  labs(x = "Synchrony Scenario",
       y = expression(paste("Median ", phi, " of Returns")),
       title = NULL) +
  geom_boxplot(size = 0.8, alpha = 0.5) +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = colPal2) +
  scale_x_discrete(breaks = waiver(), labels = c("Observed", 
                                                 expression(paste(rho["R"], 
                                                                  " = 0.05")),
                                                 expression(paste(rho["R"], 
                                                                  " = 0.50")),
                                                 expression(paste(rho["R"], 
                                                                  " = 0.75"))
  )) +
  theme_sleekX(position = "standard", axisSize = axSize * 0.9)

png(file = paste(here(),"/figs/Fig2_SynchCompBoxPlots.png", sep = ""),
    height = 4, width = 3, units = "in", res = 600)
ggarrange(q2, p2, nrow = 2, ncol = 1, heights = c(1, 1.2))
dev.off()


#_________________________________________________________________________
# Generate CU-specific spawner abundance violin plots for Bowron and Chilko
selectedCUs <- c("Bwrn" , "Chlk")
nCUs <- length(selectedCUs)
colPal <- c("#e41a1c", "#4daf4a")
#make DF to contain CU-specific benchmark estimates from sim run
bmDat <- data.frame(cu = selectedCUs, 
                    highBM = NA,
                    lowBM  = NA)

### Sigma data
plotList <- lapply(seq_along(dirNames), function(i) { 
  cuList <- genOutputList(dirNames[i], selectedCUs = selectedCUs,
                          agg = FALSE)[["medSynch_TAM"]]
  nTrials <- nrow(cuList[["medSpawners"]])
  spwnDat <- data.frame(om = rep(omNames[i], length.out = nTrials * nCUs),
                        sigma = rep(sigNames[i], length.out = nTrials * nCUs)
  )
  spwn <- cuList[["medSpawners"]] %>%
    as.data.frame() %>%
    gather(key = cu, value = spawners)
  spwn$cu <- plyr::revalue(as.factor(spwn$cu), c("V1" = selectedCUs[1],
                                                 "V2" = selectedCUs[2]))
  spwn$lowBM <- rep(cuList[["meanSGen"]], each = nTrials)
  spwn$highBM <- rep(0.8*cuList[["meanSMSY"]], each = nTrials)
  plotDat <- cbind(spwnDat, spwn)
  plotDat %>%
    mutate(cu = recode(cu, "Bwrn" = "Bowron", "Chlk" = "Chilko",
                       .default = levels(cu)),
           om = recode(om, "ref" = "Reference", "lowA" = "Low",
                       "lowStudT" = "Low +\nHeavy Tails", 
                       .default = levels(om)))
})

sigPlotDat <- do.call(rbind, plotList) 
row.names(sigPlotDat) <- NULL

### Synch data
fullPlotList <- list()
for (i in c(2,5,8)) { #make dataframe
  cuList <- genOutputList(dirNames[i], selectedCUs = selectedCUs, agg = FALSE)
  plotList <- lapply(cuList, function(h) {
    plotDat <- NULL
    nTrials <- nrow(h[["medSpawners"]])
    spwnDat <- data.frame(om = rep(omNames[i], length.out = nTrials * nCUs),
                          synch = rep(unique(h[["opMod"]]), 
                                      length.out = nTrials * nCUs)
    )
    spwn <- h[["medSpawners"]] %>%
      as.data.frame() %>%
      gather(key = cu, value = spawners)
    spwn$cu <- plyr::revalue(as.factor(spwn$cu), c("V1" = selectedCUs[1],
                                                   "V2" = selectedCUs[2]))
    spwn$lowBM <- rep(h[["meanSGen"]], each = nTrials)
    spwn$highBM <- rep(0.8*h[["meanSMSY"]], each = nTrials)
    plotDat <- rbind(plotDat, cbind(spwnDat, spwn))
    return(plotDat)
  })
  plotDat2 <- do.call(rbind, plotList)
  fullPlotList[[paste0(omNames[i])]] <- plotDat2
}

synchPlotDat <- do.call(rbind, fullPlotList) %>%
  mutate(cu = recode(cu, "Bwrn" = "Bowron", "Chlk" = "Chilko", 
                     .default = levels(cu)), 
         om = recode(om, "ref" = "Reference", "lowA" = "Low",
                     "lowStudT" = "Low +\nHeavy Tails", .default = levels(om)),
         synch = factor(factor(synch),
                        levels = c("lowSynch", "medSynch", "highSynch")))
row.names(synchPlotDat) <- NULL

# Identify y axis limits using so facets have same bounds (synch
# instead of sigma because it has a larger range)
trimSigDat <- sigPlotDat %>% 
  select(cu, spawners)
trimSynchDat <- synchPlotDat %>% 
  select(cu, spawners)
yLimDat <- rbind(trimSigDat, trimSynchDat) %>% 
  group_by(cu) %>% 
  summarize(minY = min(spawners),
            maxY = max(spawners))
sigPlotDat <- sigPlotDat %>% 
  left_join(yLimDat)
synchPlotDat <- synchPlotDat %>% 
  left_join(yLimDat)

q <- ggplot(sigPlotDat, aes(x = om, y = spawners, fill = cu, alpha = sigma)) +
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.75)) +
  geom_hline(plotDat, mapping = aes(yintercept = highBM), linetype = 2) +
  scale_alpha_manual(values = c(1, 0.575, 0.15),
                     guide = FALSE) +
  scale_fill_manual(values = colPal, guide = FALSE) +
  theme_sleekX(axisSize = axSize - 1) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Component Variability") +
  facet_wrap(~cu, scales = "free_y", nrow = 2) +
  geom_blank(aes(y = minY)) +
  geom_blank(aes(y = maxY))

p <- ggplot(synchPlotDat, aes(x = om, y = spawners, fill = cu, alpha = synch)) +
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.75)) +
  geom_hline(plotDat3, mapping = aes(yintercept = highBM), linetype = 2) +
  scale_alpha_manual(name = "Operating Model", values = c(1, 0.575, 0.15),
                     labels = c("lowSynch" = "Low",
                                "medSynch" = "Ref/Mod.",
                                "highSynch" = "High")) +
  scale_fill_manual(name = "Conservation Unit", values = colPal) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  theme_sleekX(axisSize = axSize - 1, 
               legendSize = legSize * 1.1)+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("Synchrony") +
  facet_wrap(~cu, scales = "free_y", nrow = 2) +
  geom_blank(aes(y = minY)) +
  geom_blank(aes(y = maxY))

png(file = paste(here(),"/figs/Fig5_spawnerViolin.png", sep = ""),
    height = 4.5, width = 6, units = "in", res = 600)
fig <- ggarrange(q, p, nrow = 1, ncol = 2, widths = c(1, 1.55))
annotate_figure(fig, 
                bottom = text_grob("Productivity Scenario", size = axSize, 
                                   color = "grey30"),
                left = text_grob("Spawner Abundance (millions)", size = axSize,
                                 rot = 90, color = "grey30"))
dev.off()
