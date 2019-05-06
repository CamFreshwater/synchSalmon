#*******************************************************************************
# runModel_multiHCR.R
# Date revised: May 6, 2019
# Inputs: recoverySim.R
# Outputs: pdf plots
# Explainer: Runs closed loop sim for synch manuscript with reference prod
# and alternative multistock HCR that is less conservative
# CONCLUSION: Impacts are relatively modest (similar to decline in productivity)
#*******************************************************************************

# Check if required packages are installed and run
listOfPackages <- c("plyr", "here", "parallel", "doParallel", "foreach", 
                    "reshape2", "tidyverse", "gsl", "tictoc", "stringr", 
                    "synchrony", "zoo", "Rcpp", "RcppArmadillo", "sn", 
                    "sensitivity", "mvtnorm", "forcats", "ggpubr", "viridis", 
                    "samSim")
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data/sox/fraserOMInputs_varyCorrMultiHCR.csv"), 
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/sox/fraserCUParsFRP.csv"), stringsAsFactors = F) 
# %>% 
#   filter(!stkName == "E. Shuswap") #remove for now because recursive est not avail
srDat <- read.csv(here("data/sox/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/sox/fraserCatchDatTrim.csv"), 
                     stringsAsFactors = F)
# ricPars <- read.csv(here("data/sox/hiddenData/trimRecursiveRickerMCMCPars.csv"), 
#                     stringsAsFactors = F)
# larkPars <- read.csv(here("data/sox/hiddenData/trimRecursiveLarkinMCMCPars.csv"), 
#                      stringsAsFactors = F)
ricPars <- read.csv(here("data/sox/pooledRickerMCMCPars.csv"),
                    stringsAsFactors = F)
larkPars <- read.csv(here("data/sox/pooledLarkinMCMCPars.csv"),
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/sox/tamRefPts.csv"), stringsAsFactors = F)


### SET UP MODEL RUN -----------------------------------------------------

## Define simulations to be run
nTrials <- 400

simParTrim <- simPar %>% 
  filter(prodRegime == "low") %>% 
  mutate(scenario = paste("lowProd", scenario, nameMP, sep = "_"))

scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, sep = "_"))

# recoverySim(simParTrim[12, ], cuPar, catchDat = catchDat, srDat = srDat,
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
mpNames <- rep(c("TAM", "genPA"), each = 3)
sigNames <- rep(c("low", "med", "high"), length.out = length(mpNames))

plotDat1 = NULL
for(h in seq_along(dirNames)) {
  agList <- genOutputList(dirNames[h], agg = TRUE)
  keyVar <- sapply(agList, function(x) unique(x$keyVar))
  plotOrder <- sapply(agList, function(x) unique(x$plotOrder))
  singleScen = NULL
  for (i in seq_along(vars)) {
    dum <- data.frame(sigma = rep(sigNames[h], length.out = length(agList)),
                      mp = rep(mpNames[h], length.out = length(agList)),
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
             mp = as.factor(mpNames[h]),
             scenID = as.factor(paste(mpNames[h], sigNames[h], datList$nameOM, 
                                      sep = "_"))) %>% 
      filter(yr > datList$nPrime)
  }) #iterate across different OMs within a scenario
  dat <- do.call(rbind, listBySynch) #merge synch list elements into DF
})
fullDat <- do.call(rbind, listByOM)

# standardize data relative to low synch/low sigma in reference OM
#First calculate within CU medians for low synch/low sigma/reference dataset
lowAggV <- fullDat %>% 
  filter(scenID == "TAM_low_lowSynch") %>% 
  dplyr::group_by(cu, trial, sigma, synch, mp, scenID) %>% 
  dplyr::summarize(medR = median(recRY)) %>%
  dplyr::group_by(cu, sigma, synch, mp, scenID) %>% 
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
      group_by(cu, trial, sigma, synch, mp) %>% 
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
  group_by(trial, sigma, synch, mp) %>% 
  dplyr::summarize(meanStdRecRY = mean(stdMedR))

#finally calculate medians and quantiles to generate dataset equivalent to 
#plotdat
finalOut <- temp %>% 
  group_by(sigma, synch, mp) %>% 
  mutate(var = "stdRecRY",
         mean = mean(meanStdRecRY),
         medn = median(meanStdRecRY),
         lowQ = qLow(meanStdRecRY),
         highQ = qHigh(meanStdRecRY)) %>%
  select(sigma, mp, var, synch, medn, lowQ, highQ) 

stdFullDat <- finalOut %>% 
  as.data.frame() %>% 
  mutate(synch = recode(synch, "lowSynch" = "low", "medSynch" = "med",
                        "highSynch" = "high", .default = levels(synch)))
# saveRDS(stdFullDat, here("outputs/generatedData/stdMedRecRY.rds"))

plotDat <- rbind(plotDat1, stdFullDat) 


# Default plotting variable settings
dotSize = 3.25; lineSize = 0.8; legSize = 0.8; axSize = 10; facetSize = 0.95

# Default color palette
colPal <- viridis(length(levels(plotDat$synch)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$synch)

consVars <- c("medRecRY", "stdRecRY", "ppnMixedOpen", "ppnCUUpper")
consYLabs <- c("Agg. Return\nAbundance", "CU-Specific Std.\nReturn Abundance", 
               "Prop. MUs\nAbove Esc. Goal", "Prop. CUs\nAbove Benchmark")
#make dataframe of labels to annotate facets
consLabs <- data.frame(mp = rep(factor(unique(plotDat$mp),
                                       levels = unique(plotDat$mp)),
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
    facet_wrap(~mp, scales = "fixed", ncol = 4, nrow = 1)
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
catchLabs <- data.frame(mp = rep(factor(unique(plotDat$mp),
                                        levels = unique(plotDat$mp)),
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
    facet_wrap(~mp, scales = "fixed", ncol = 4, nrow = 1)
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

png(file = paste(here(),"/figs/multiHCR/lowProd_consGroupedPlots.png", sep = ""),
    height = 5.5, width = 5.25, units = "in", res = 600)
ggarrange(consPlots[[1]], consPlots[[2]], consPlots[[3]], consPlots[[4]],
          ncol = 1, nrow = 4, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,0.9,0.9,1.15))
dev.off()
png(file = paste(here(),"/figs/multiHCR/lowProd_catchGroupedPlots.png", sep = ""),
    height = 4.5, width = 5.25, units = "in", res = 600)
ggarrange(catchPlots[[1]], catchPlots[[2]], catchPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,0.9,1.2))
dev.off()
