#*******************************************************************************
# runModelRobustOM_synch.R
# Date revised: ONGOING
# Inputs: recoverySim.R
# Outputs: pdf plots
# Explainer: Runs closed loop simulation model with different OMs and constant 
# MP
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
ricPars <- read.csv(here("data/sox/pooledRickerMCMCPars.csv"), stringsAsFactors = F)
larkPars <- read.csv(here("data/sox/pooledLarkinMCMCPars.csv"), stringsAsFactors = F)
tamFRP <- read.csv(here("data/sox/tamRefPts.csv"), stringsAsFactors = F)


### SET UP MODEL RUN -----------------------------------------------------

## Define simulations to be run
nTrials <- 1200

## General robustness runs
simParTrim <- subset(simPar, scenario %in% c("lowSig", "medSig", "highSig"
                                             # , 
                                             # "lowSigLowA", "medSigLowA", 
                                             # "highSigLowA", "lowSigStudT",
                                             # "medSigStudT", "highSigStudT",
                                             # "lowSigLowStudT", "medSigLowStudT",
                                             # "highSigLowStudT"
                                             ))

scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species),
                                                sep = "_"))

# recoverySim(simParTrim[31, ], cuPar, catchDat = catchDat, srDat = srDat,
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
## Ugly custom code to generate proportional dot plots split across three 
# different productivity scenarios
vars <- c("medRecRY", "ppnYrsCUsUpper", "ppnCUExtant",
          "medCatch", "ppnYrsHighCatch", "stabilityCatch")
omNames <- rep(c("ref", "lowA", "studT", "lowStudT"), each = 3)
sigNames <- rep(c("low", "med", "high"), length.out = length(omNames))

plotDat = NULL
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
                      cat = as.factor(plotOrder),
                      avg = sapply(agList, function(x) median(x[,vars[i]])),
                      lowQ = sapply(agList, function(x) qLow(x[,vars[i]])),
                      highQ = sapply(agList, function(x) qHigh(x[,vars[i]])),
                      row.names = NULL
    )
    singleScen <- rbind(singleScen, dum)
  }
  rownames(singleScen) <- c()
  #merge multiple scenarios into one dataframe
  plotDat <- rbind(plotDat, singleScen) 
}
plotDat <- plotDat %>%
  mutate(cat = recode(cat, "1" = "low", "2" = "med", "3" = "high", 
                      .default = levels(cat)),
         om = recode(om, "ref" = "Reference Prod.", "lowA" = "Low Prod.",
                     "studT" = "Ref. Heavy Tails", 
                     "lowStudT" = "Low Heavy Tails", .default = levels(om))
  )

#Save summary data to pass to Rmd
write.csv(plotDat, here("outputs/generatedData", "summaryTable_noSkew.csv"))


colPal <- viridis(length(levels(plotDat$cat)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$cat)
dotSize = 3; lineSize = 0.8; legSize = 0.7; axSize = 10; facetSize = 0.95

consVars <- c("medRecRY", "ppnYrsCUsUpper", "ppnCUExtant")
consYLabs <- c("Return\nAbundance", "Prop. CUs\nAbove Benchmark", 
               "Prop. CUs\nExtant")
## Correct based on column number
# consLabs <- data.frame(om = rep(factor(unique(plotDat$om), 
#                                        levels = unique(plotDat$om)),
#                                 each = 3),
#                        lab = c("a)", "d)", "g)", "b)", "e)", "h)", "c)", "f)", 
#                                "i)"),
#                        var = rep(factor(consVars, levels = unique(vars)), 
#                                  times = 3)
# )#make dataframe of labels to annotate facets
consPlots <- lapply(seq_along(consVars), function(i) {
  temp <- plotDat %>%
    filter(var == consVars[i])
  q <- ggplot(temp, aes(x = sigma, y = avg, ymin = lowQ, ymax = highQ,
                        color = cat)) +
    labs(x = "Component Variance", y = consYLabs[i], 
         color = "Sim.\nParameter\nValue") +
    geom_pointrange(fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.65)) +
    scale_x_discrete(labels = c("lowSigma" = "Low",
                                "medSigma" = "Reference",
                                "highSigma" = "High")) +
    scale_colour_manual(name = "Synchrony", values = colPal,
                        labels = c("low" = "Low",
                                   "med" = "Moderate",
                                   "high" = "High")) +
    # geom_text(data = consLabs %>% filter(var == consVars[i]),
    #           mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab,
    #                         hjust = 0.5, vjust = 0),
    #           show.legend = FALSE, inherit.aes = FALSE) +
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

catchVars <- c("medCatch", "ppnYrsHighCatch", "stabilityCatch")
catchYLabs <- c("Catch\nAbundance", "Prop. Years Above\nCatch Threshold",
                "Catch\nStability")
# catchLabs <- data.frame(om = rep(factor(unique(plotDat$om), 
#                                         levels = unique(plotDat$om)),
#                                  each = 3),
#                         lab = c("a)", "d)", "g)", "b)", "e)", "h)", "c)", "f)", 
#                                 "i)"),
#                         var = rep(factor(catchVars, levels = unique(vars)), 
#                                   times = 3)
# )#make dataframe of labels to annotate facets
catchPlots <- lapply(seq_along(catchVars), function(i) {
  temp <- plotDat %>%
    filter(var == catchVars[i])
  q <- ggplot(temp, aes(x = sigma, y = avg, ymin = lowQ, ymax = highQ,
                        color = cat)) +
    labs(x = "Component Variance", y = catchYLabs[i], 
         color = "Sim.\nParameter\nValue") +
    geom_pointrange(fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.65)) +
    scale_x_discrete(labels = c("lowSigma" = "Low",
                                "medSigma" = "Reference",
                                "highSigma" = "High")) +
    scale_colour_manual(name = "Synchrony", values = colPal,
                        labels = c("low" = "Low",
                                   "med" = "Moderate",
                                   "high" = "High")) +
    # geom_text(data = catchLabs %>% filter(var == catchVars[i]),
    #           mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab,
    #                         hjust = 0.5, vjust = 0),
    #           show.legend = FALSE, inherit.aes = FALSE) +
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

png(file = paste(here(),"/figs/Fig3_consGroupedPlots_noSkew.png", sep = ""),
    height = 4.5, width = 6, units = "in", res = 600)
ggarrange(consPlots[[1]], consPlots[[2]], consPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1.2))
dev.off()
png(file = paste(here(),"/figs/Fig4_catchGroupedPlots_noSkew.png", sep = ""),
    height = 4.5, width = 6, units = "in", res = 600)
ggarrange(catchPlots[[1]], catchPlots[[2]], catchPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1.2))
dev.off()


#_________________________________________________________________________
# Generate box plots (originally time series) of component CV, synch and ag CV;
# Note that to make comparable to retrospective analysis only show 10 
# plotList = vector("list", length = length(dirNames))
# arrayNames <- sapply(dirNames, function(x) {
#   list.files(paste(here("outputs/simData"), x, sep="/"),
#              pattern = "\\Arrays.RData$")
# })
# 
# tic("runParallel")
# Ncores <- detectCores()
# cl <- makeCluster(Ncores - 2) #save two cores
# registerDoParallel(cl)
# clusterEvalQ(cl, c(library(here), library(synchrony), library(zoo),
#                    library(parallel), library(doParallel), library(foreach),
#                    library(samSim)))
# newAgTSList <- lapply(seq_along(dirNames), function (h) {
#   #export custom function and objects
#   clusterExport(cl, c("dirNames", "arrayNames", "calcSynchMetrics", "wtdCV",
#                       "genOutputList", "h"), envir = environment())
#   listSynchLists <- parLapply(cl, 1:length(arrayNames[, h]), function(x) {
#     datList <- readRDS(paste(here("outputs/simData"), dirNames[h],
#                              arrayNames[x, h], sep = "/"))
#     synchList <- calcSynchMetrics(datList, corr = FALSE,
#                                   weight = TRUE, windowSize = 12)
#     synchList <- c(datList$nameOM, synchList)
#     names(synchList)[1] <- "opMod"
#     return(synchList)
#   }) #iterate across different OMs within a scenario
#   # pull list of time series metrics estimated internally, then merge with
#   # synch list generated above based on common op model
#   agTSList <- genOutputList(dirNames[h], agg = TRUE, aggTS = TRUE)
#   for(j in seq_along(agTSList)) {
#     for(k in seq_along(listSynchLists)) {
#       om <- agTSList[[j]]$opMod
#       if (listSynchLists[[k]]$opMod == om) {
#         agTSList[[j]] <- c(agTSList[[j]], listSynchLists[[k]][-1])
#       }
#     }
#   }
#   plotList[[h]] <- agTSList
# }) #iterate across different scenarios
# names(newAgTSList) <- dirNames
# stopCluster(cl)
# toc()

## Save
# saveRDS(newAgTSList, here("outputs/generatedData/synchTS/synchTSList.rda"))
# newAgTSList <- readRDS(here("outputs/generatedData/synchTS/synchTSList.rda"))

### Manipulate lists to create plottable data structure 
## In this case that is one median value per year per unique combo of OMs
# omNames <- rep(c("ref", "skewN", "studT"), each = 3)
# sigNames <- rep(c("lowSigma", "medSigma", "highSigma"), length.out = 9)
# fullList <- sapply(seq_along(dirNames), function(h) {
#   d <- newAgTSList[[h]]
#   nYears <- d[[1]]$`nYears`
#   simLength <- d[[1]]$`nYears` - d[[1]]$`nPrime`
#   firstYear <- d[[1]]$`firstYr`
#   start <- d[[1]]$`nPrime` + firstYear
#   prodNames <- omNames[h]
#   #subset list so it contains only sigma, synch and vars of interest,
#   # calculate medians, and combine
#   trimList <- lapply(d, function(x) {
#     dat1 <- data.frame(sigmaOM = rep(sigNames[h], length.out = nYears),
#                       synchOM = rep(x[["opMod"]], length.out = nYears),
#                       prodOM = rep(prodNames, length.out = nYears),
#                       year = seq(from = firstYear,
#                                  to = (firstYear + nYears - 1))
#                       ) %>%
#       mutate(sigmaOM = as.character(sigmaOM),
#              synchOM = as.character(synchOM),
#              medSynchRecBY = apply(x[["synchRecBY"]], 1, median),
#              medCompCVRecBY = apply(x[["compCVRecBY"]], 1, median)
#       )
#     dat1[dat1$year < start, c("sigmaOM", "synchOM")] <- "obs"
#     dat2 <- dat1 %>%
#       mutate(sigmaOM = factor(factor(sigmaOM),
#                               levels = c("obs", "lowSigma", "medSigma",
#                                          "highSigma")),
#              synchOM = factor(factor(synchOM),
#                               levels = c("obs", "lowSynch", "medSynch",
#                                          "highSynch"))
#       )
#     return(dat2)
#   })
# })
# plotDat <- do.call(rbind, fullList)

# saveRDS(plotDat, file = here("outputs/generatedData/synchTS/fullSynch_3OMs.rda"))
plotDat <- readRDS(file = here("outputs/generatedData/synchTS/fullSynch_3OMs.rda"))
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
                                                                  sigma)),
                                                 expression(paste(sigma)),
                                                 expression(paste("1.25", 
                                                                  sigma))
  )) +
  theme_sleekX(position = "standard", axisSize = axSize)

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
                                                 expression(paste(rho, 
                                                                  " = 0.05")),
                                                 expression(paste(rho, 
                                                                  " = 0.50")),
                                                 expression(paste(rho, 
                                                                  " = 0.75"))
  )) +
  theme_sleekX(position = "standard", axisSize = axSize)


png(file = paste(here(),"/figs/Fig2_SynchCompBoxPlots.png", sep = ""),
    height = 5, width = 3, units = "in", res = 600)
ggarrange(q2, p2, nrow = 2, ncol = 1, heights = c(1, 1.2))
dev.off()


#_________________________________________________________________________
# Generate CU-specific spawner abundance violin plots for Bowron and Chilko
selectedCUs <- c("Bwrn" , "Chlk")
nCUs <- length(selectedCUs)
colPal <- c("#7fc97f", "#beaed4", "#fdc086")
#make DF to contain CU-specific benchmark estimates from sim run
bmDat <- data.frame(cu = selectedCUs, 
                    highBM = NA,
                    lowBM  = NA)
plotDat <- NULL
for (i in seq_along(dirNames)) { #make dataframe
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
  plotDat <- rbind(plotDat, cbind(spwnDat, spwn))
  plotDat <- plotDat %>%
    mutate(cu = recode(cu, "Bwrn" = "Bowron", "Chlk" = "Chilko",
                       .default = levels(cu)),
           om = recode(om, "ref" = "Reference Prod.", "lowA" = "Low Prod.",
                       "studT" = "Ref. Heavy Tails", 
                       "lowStudT" = "Low Heavy Tails", .default = levels(om)))
}

q <- ggplot(plotDat, aes(x = om, y = spawners, fill = cu, alpha = sigma)) +
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.75)) +
  geom_hline(plotDat, mapping = aes(yintercept = highBM), linetype = 2) +
  scale_alpha_manual(name = "Component CV", values = c(1, 0.65, 0.3),
                     labels = c("low" = "Low",
                                "med" = "Ref.",
                                "high" = "High")) +
  scale_fill_manual(values = colPal, guide = FALSE) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  labs(y = "Median Spawner Abundance (millions)",
       x = "Productivity Scenario") +
  theme_sleekX(facetSize = 1.1, axisSize = axSize - 1, 
               legendSize = legSize * 0.75) +
  facet_wrap(~cu, scales = "free_y", nrow = 2)


png(file = here("figs/Fig6_spawnerViolinSigma.png"),
    height = 3.5, width = 3, units = "in", res = 600)
print(q)
dev.off()


fullPlotList <- list()
for (i in c(2,5,8,11)) { #make dataframe
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

plotDat3 <- do.call(rbind, fullPlotList) %>%
  mutate(cu = recode(cu, "Bwrn" = "Bowron", "Chlk" = "Chilko", 
                     .default = levels(cu)), 
         om = recode(om, "ref" = "Reference Prod.", "lowA" = "Low Prod.",
                     "studT" = "Ref. Heavy Tails",
                     "lowStudT" = "Low Heavy Tails", .default = levels(om)))
row.names(plotDat3) <- NULL

p <- ggplot(plotDat3, aes(x = om, y = spawners, fill = cu, alpha = synch)) +
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.75)) +
  geom_hline(plotDat3, mapping = aes(yintercept = highBM), linetype = 2) +
  scale_alpha_manual(name = "Synchrony", values = c(1, 0.65, 0.3),
                     labels = c("lowSynch" = "Low",
                                "medSynch" = "Mod.",
                                "highSynch" = "High")) +
  scale_fill_manual(values = colPal, guide = FALSE) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  labs(y = "Median Spawner Abundance (millions)",
       x = "Productivity Scenario") +
  theme_sleekX(facetSize = 1.1, axisSize = axSize - 1, 
               legendSize = legSize * 0.75) +
  facet_wrap(~cu, scales = "free_y", nrow = 2)


png(file = paste(here(),"/figs/SFig_spawnerViolinSynch.png", sep = ""),
    height = 3.5, width = 3, units = "in", res = 300)
print(p)
dev.off()


fullPlotList <- list()
for (i in c(2,5,8,11)) { #make dataframe
  cuList <- genOutputList(dirNames[i], selectedCUs = selectedCUs, agg = FALSE)
  plotList <- lapply(cuList, function(h) {
    plotDat <- NULL
    nTrials <- nrow(h[["ppnYrsUpper"]])
    spwnDat <- data.frame(om = rep(omNames[i], length.out = nTrials * nCUs),
                          synch = rep(unique(h[["opMod"]]), 
                                      length.out = nTrials * nCUs)
    )
    spwn <- h[["ppnYrsUpper"]] %>%
      as.data.frame() %>%
      gather(key = cu, value = spawners)
    spwn$cu <- plyr::revalue(as.factor(spwn$cu), c("V1" = selectedCUs[1],
                                                   "V2" = selectedCUs[2]))
    plotDat <- rbind(plotDat, cbind(spwnDat, spwn))
    return(plotDat)
  })
  plotDat2 <- do.call(rbind, plotList)
  fullPlotList[[paste0(omNames[i])]] <- plotDat2
}

##### END PRIMARY ANALYSIS #####


#______________________________________________________________________________
#Stand alone analysis to see whether there are increases or decreases in median spawner
# abundance associated with greater CVc
sigNames <- c("low", "med", "high")
plotAgDat = NULL
plotCUDat = NULL
for (i in 1:3) {
  cuList <- genOutputList(dirNames[i], selectedCUs = NULL,
                          agg = FALSE)[["medSynch_TAM"]]
  nTrials <- nrow(cuList[["medRecRY"]])
  nCUs <- ncol(cuList[["medRecRY"]])
  spwnDat <- data.frame(sigma = rep(sigNames[i], length.out = nTrials * nCUs),
                        trial = rep(seq(from = 1, to = nTrials, by = 1),
                                    times = nCUs)
  )
  spwn <- cuList[["medRecRY"]] %>%
    as.data.frame() %>%
    gather(key = cu, value = recruits)
  plotCUDat <- rbind(plotCUDat, cbind(spwnDat, spwn))
  
  agDat <- genOutputList(dirNames[i], agg = TRUE)[["medSynch_TAM_aggDat.csv"]] %>%
    select(trial, medSpawners, medRecRY) %>%
    mutate(sigma = sigNames[i])
  plotAgDat <- rbind(plotAgDat, agDat)
}
plotCUDat %>%
  group_by(trial, sigma) %>%
  summarise(agRec = sum(recruits)) %>%
  group_by(sigma) %>%
  summarise(mean(agRec))
plotAgDat %>%
  group_by(sigma) %>%
  summarise(mean(medRecRY))

