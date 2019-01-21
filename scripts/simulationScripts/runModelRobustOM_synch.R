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
nTrials <- 25

## General robustness runs
simParTrim <- subset(simPar,
                     scenario == "lowSig" | scenario == "medSig" |
                       scenario == "highSig" | scenario == "lowSigSkew" |
                       scenario == "medSigSkew" | scenario == "highSigSkew" |
                     scenario == "lowSigSkewT" | scenario == "medSigSkewT" | 
                       scenario == "highSigSkewT" 
                     # | scenario == "lowSigLowA" |
                     # scenario == "medSigLowA" | scenario == "highSigLowA"
                     )

scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species),
                                                sep = "_"))

# recoverySim(simParTrim[1, ], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP, 
#             dirName = "test", nTrials = 5, multipleMPs = FALSE)
# for(i in seq_along(dirNames)) {
# d <- subset(simParTrim, scenario == scenNames[i])
# simsToRun <- split(d, seq(nrow(d)))
# lapply(simsToRun, function(x) {
#   recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
#               ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirNames[i],
#               nTrials=125, multipleMPs=TRUE)
# })
# }

for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 3) #save two cores
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
                nTrials = nTrials, multipleMPs = FALSE)
  })
  stopCluster(cl) #end cluster
  toc()
}


#_________________________________________________________________________
## Ugly custom code to generate proportional dot plots split across three 
# different productivity scenarios
vars <- c("medRecRY", "ppnCUUpper", "ppnCUExtant",
          "medCatch", "ppnYrsHighCatch", "stabilityCatch")
omNames <- rep(c("ref", "skewN", "skewT"), each = 3)
sigNames <- rep(c("low", "med", "high"), length.out = 9)

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
         om = recode(om, "ref" = "Reference Prod.", "skewN" = "Low Prod.",
                       "skewT" = "Very Low Prod.",
                       .default = levels(om))
  )

colPal <- viridis(length(levels(plotDat$sigma)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$sigma)
dotSize = 3.5; lineSize = 0.8; legendSize = 14

consVars <- c("medRecRY", "ppnCUUpper", "ppnCUExtant")
consYLabs <- c("Return\nAbundance", "Prop. CUs Above Benchmark", 
               "Prop. CUs\nExtant")
consLabs <- data.frame(om = rep(factor(unique(plotDat$om), 
                                       levels = unique(plotDat$om)),
                                each = 3),
                       lab = c("a)", "d)", "g)", "b)", "e)", "h)", "c)", "f)", 
                               "i)"),
                       var = rep(factor(consVars, levels = unique(vars)), 
                                 times = 3)
)#make dataframe of labels to annotate facets
consPlots <- lapply(seq_along(consVars), function(i) {
  temp <- plotDat %>%
    filter(var == consVars[i])
  q <- ggplot(temp, aes(x = sigma, y = avg, ymin = lowQ, ymax = highQ,
                        color = cat)) +
    labs(x = "Component Variance", y = consYLabs[i], 
         color = "Sim.\nParameter\nValue") +
    geom_pointrange(fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.5)) +
    scale_x_discrete(labels = c("low" = "Low",
                                "med" = "Reference",
                                "high" = "High")) +
    scale_colour_manual(name = "Synchrony", values = colPal,
                        labels = c("low" = "Low",
                                   "med" = "Moderate",
                                   "high" = "High")) +
    # scale_shape_manual(name = "Sigma", breaks = c("low", "med", "high"),
    #                    values = c(16, 17, 18), guide = FALSE) +
    geom_text(data = consLabs %>% filter(var == consVars[i]),
              mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab, 
                            hjust = 1.5, vjust = 0),
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~om, scales = "fixed")
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = 1.15, axisSize = 13)
  }
  if (i == 2) {
    q <- q + theme_sleekX(position = "mid", legendSize = 1.15, axisSize = 13)
  }
  if (i == 3) {
    q <- q + theme_sleekX(position = "bottom", legendSize = 1.15, axisSize = 13)
  }
  return(q)
})

catchVars <- c("medCatch", "ppnYrsHighCatch", "stabilityCatch")
catchYLabs <- c("Catch\nAbundance", "Prop. Years Above\nCatch Threshold",
                "Catch\nStability")
catchLabs <- data.frame(om = rep(factor(unique(plotDat$om), 
                                        levels = unique(plotDat$om)),
                                 each = 3),
                        lab = c("a)", "d)", "g)", "b)", "e)", "h)", "c)", "f)", 
                                "i)"),
                        var = rep(factor(catchVars, levels = unique(vars)), 
                                  times = 3)
)#make dataframe of labels to annotate facets
catchPlots <- lapply(seq_along(catchVars), function(i) {
  temp <- plotDat %>%
    filter(var == catchVars[i])
  q <- ggplot(temp, aes(x = sigma, y = avg, ymin = lowQ, ymax = highQ,
                        color = cat)) +
    labs(x = "Component Variance", y = catchYLabs[i], 
         color = "Sim.\nParameter\nValue") +
    geom_pointrange(fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.5)) +
    scale_x_discrete(labels = c("low" = "Low",
                                "med" = "Reference",
                                "high" = "High")) +
    scale_colour_manual(name = "Synchrony", values = colPal,
                        labels = c("low" = "Low",
                                   "med" = "Moderate",
                                   "high" = "High")) +
    # scale_shape_manual(name = "Sigma", breaks = c("low", "med", "high"),
    #                    values = c(16, 17, 18), guide = FALSE) +
    geom_text(data = catchLabs %>% filter(var == catchVars[i]),
              mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab,
                            hjust = 1.5, vjust = 0),
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~om, scales = "fixed")
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = 1.15, axisSize = 13)
  }
  if (i == 2) {
    q <- q + theme_sleekX(position = "mid", legendSize = 1.15, axisSize = 13)
  }
  if (i == 3) {
    q <- q + theme_sleekX(position = "bottom", legendSize = 1.15, axisSize = 13)
  }
  return(q)
})

png(file = paste(here(),"/figs/Fig4_consGroupedPlots_3OMs.png", sep = ""),
    height = 5.5, width = 8.5, units = "in", res = 300)
ggarrange(consPlots[[1]], consPlots[[2]], consPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1.2))
dev.off()
png(file = paste(here(),"/figs/Fig5_catchGroupedPlots_3OMs.png", sep = ""),
    height = 5.5, width = 8.5, units = "in", res = 300)
ggarrange(catchPlots[[1]], catchPlots[[2]],
          catchPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1.2))
dev.off()


#_________________________________________________________________________
# Generate box plots (originally time series) of component CV, synch and ag CV;
# Note that to make comparable to retrospective analysis only show 10 
plotList = vector("list", length = length(dirNames))
#matrix of array names to be passed
arrayNames <- sapply(dirNames, function(x) {
  list.files(paste(here("outputs/simData"), x, sep="/"), 
             pattern = "\\Arrays.RData$")
})
#stocks that are included in retrospective analysis 
#(note these are matrix positions not stock ID numbers)
# longStks <- c(1, 3, 4, 5, 6, 7, 8, 9, 16, 17)

tic("runParallel")
Ncores <- detectCores()
cl <- makeCluster(Ncores - 2) #save two cores
registerDoParallel(cl)
clusterEvalQ(cl, c(library(here), library(synchrony), library(zoo),
                   library(parallel), library(doParallel), library(foreach),
                   library(samSim)))
newAgTSList <- lapply(seq_along(dirNames), function (h) {
  #export custom function and objects
  clusterExport(cl, c("dirNames", "arrayNames", "calcSynchMetrics", "wtdCV",
                      "genOutputList", "longStks", "h"), envir = environment())
  listSynchLists <- parLapply(cl, 1:length(arrayNames[, h]), function(x) {
    datList <- readRDS(paste(here("outputs/simData"), dirNames[h],
                             arrayNames[x, h], sep = "/"))
    # Subset datList to only include longStks
    trimDatList <- datList
    # for (i in 4:length(datList)) {
    #   trimDatList[[i]] <- datList[[i]][ , longStks, ]
    # }
    synchList <- calcSynchMetrics(trimDatList, log = FALSE, corr = FALSE,
                                  weight = TRUE, windowSize = 12)
    synchList <- c(trimDatList$nameOM, synchList)
    names(synchList)[1] <- "opMod"
    return(synchList)
  }) #iterate across different OMs within a scenario
  # pull list of time series metrics estimated internally, then merge with
  # synch list generated above based on common op model
  agTSList <- genOutputList(dirNames[h], agg = TRUE, aggTS = TRUE)
  for(j in seq_along(agTSList)) {
    for(k in seq_along(listSynchLists)) {
      om <- agTSList[[j]]$opMod
      if (listSynchLists[[k]]$opMod == om) {
        agTSList[[j]] <- c(agTSList[[j]], listSynchLists[[k]][-1])
      }
    }
  }
  plotList[[h]] <- agTSList
}) #iterate across different scenarios
names(newAgTSList) <- dirNames
stopCluster(cl)
toc()

## Save 
#saveRDS(newAgTSList, here("outputs/generatedData/synchTS/synchTSList.rda"))
newAgTSList <- readRDS(here("outputs/generatedData/synchTS/synchTSList.rda"))

## Import tmb data to replace observed data with
modOut <- readRDS(here("outputs", "generatedData", "tmbSynchEst.rds"))
modCVc <- modOut %>% 
  filter(term == "log_cv_s") %>% 
  select(est) %>% 
  unlist()
modSynch <- modOut %>% 
  filter(term == "logit_phi") %>% 
  select(est) %>% 
  unlist()
#add NAs to front end to match length of sim time series
modCVc <- c(rep(NA, length.out = 10), modCVc) 
modSynch <- c(rep(NA, length.out = 10), modSynch)

### Manipulate lists to create plottable data structure 
## In this case that is one median value per year per unique combo of OMs
omNames <- rep(c("ref", "skewN", "skewT"), each = 3)
sigNames <- rep(c("lowSigma", "medSigma", "highSigma"), length.out = 9)
fullList <- sapply(seq_along(dirNames), function(h) {
  d <- newAgTSList[[h]]
  nYears <- d[[1]]$`nYears`
  simLength <- d[[1]]$`nYears` - d[[1]]$`nPrime`
  firstYear <- d[[1]]$`firstYr`
  start <- d[[1]]$`nPrime` + firstYear
  prodNames <- omNames[h]
  #subset list so it contains only sigma, synch and vars of interest,
  # calculate medians, and combine
  trimList <- lapply(d, function(x) {
    dat1 <- data.frame(sigmaOM = rep(sigNames[h], length.out = nYears),
                      synchOM = rep(x[["opMod"]], length.out = nYears),
                      prodOM = rep(prodNames, length.out = nYears),
                      year = seq(from = firstYear, to = (firstYear + nYears - 1))
                      ) %>%
      mutate(sigmaOM = as.character(sigmaOM),
             synchOM = as.character(synchOM),
             medSynchRecBY = apply(x[["synchRecBY"]], 1, median),
             medCompCVRecBY = apply(x[["compCVRecBY"]], 1, median),
             medCorrRecBY = apply(x[["corrRecBY"]], 1, median)
      )
    dat1[dat1$year < start, c("sigmaOM", "synchOM")] <- "obs"
    dat2 <- dat1 %>%
      mutate(sigmaOM = factor(factor(sigmaOM), levels = c("obs", "lowSigma",
                                                          "medSigma", "highSigma")),
             synchOM = factor(factor(synchOM), levels = c("obs", "lowSynch",
                                                          "medSynch", "highSynch"))
      ) 
    # %>%
    #   filter(!is.na(medSynchRecBY)) #remove yrs where obs synch couldn't be calc
    # replace observed data with estimates from TMB fitting
    dat2[dat2$sigmaOM == "obs" & dat2$synchOM == "obs", ]$medSynchRecBY <- modSynch
    dat2[dat2$sigmaOM == "obs" & dat2$synchOM == "obs", ]$medCompCVRecBY <- modCVc
    return(dat2)
  })
})
plotDat <- do.call(rbind, fullList)


# saveRDS(plotDat, file = here("outputs/generatedData/synchTS/fullSynchTSDF_3OMs.rda"))
plotDat <- readRDS(file = here("outputs/generatedData/synchTS/fullSynchTSDF_3OMs.rda"))
start <- plotDat %>%
  filter(!sigmaOM == "obs") %>%
  summarise(min(year))
start <- start[[1]]

dum <- plotDat %>%
  filter(synchOM == "medSynch" | synchOM == "obs",
         prodOM == "ref") %>% 
  mutate(sigmaOM = recode(sigmaOM, "obs" = "Observed", 
                         "lowSigma" = "Low", 
                         "medSigma" = "Moderate",
                         "highSigma" = "High")
  )
colPal <- c("black", viridis(length(unique(dum$sigmaOM)) - 1, begin = 0, end = 1))
names(colPal) <- levels(dum$sigmaOM)
### Instead of plotting time series, make faceted box plot
q2 <- ggplot(dum, aes(x = sigmaOM, y = medCompCVRecBY, fill = sigmaOM)) +
  labs(x = "Component Variability Scenario", y = "Median CVc of Returns", 
       title = NULL) +
  geom_boxplot(size = 1, alpha = 0.5) +
  geom_vline(xintercept = start, color = "black", linetype = 3, size = 1) +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(name = "Operating Model", values = colPal) +
  theme_sleekX(position = "standard", legendSize = 0.9, axisSize = 12)

dum2 <- plotDat %>%
  dplyr::filter(sigmaOM == "medSigma" | sigmaOM == "obs",
                prodOM == "ref") %>% 
  mutate(synchOM = recode(synchOM, "obs" = "Observed", 
                          "lowSynch" = "Low", 
                          "medSynch" = "Moderate",
                          "highSynch" = "High")
  )
colPal2 <- c("black", viridis(length(unique(dum2$synchOM)) - 1, begin = 0, 
                              end = 1))
names(colPal2) <- levels(dum2$synchOM)
p2 <- ggplot(dum2, aes(x = synchOM, y = medSynchRecBY, fill = synchOM)) +
  labs(x = "Synchrony Scenario",
       y = expression(paste("Median, ", phi, " of Returns")),
       title = NULL) +
  geom_boxplot(size = 1, alpha = 0.5) +
  geom_vline(xintercept = start, color = "black", linetype = 3, size = 1) +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = colPal2) +
  theme_sleekX(position = "standard", legendSize = 0.9, axisSize = 12)


png(file = paste(here(),"/figs/Fig3_SynchCompBoxPlots.png", sep = ""),
    height = 6, width = 4, units = "in", res = 300)
ggarrange(q2, p2, nrow = 2, ncol = 1, heights = c(1, 1.2))
dev.off()

boxplot(modOut$est ~ modOut$term)

#_________________________________________________________________________
# Generate CU-specific spawner abundance violin plots for Bowron and Chilko
selectedCUs <- c("Bwrn" , "Chlk")
nCUs <- length(selectedCUs)
colPal <- c("#7fc97f", "#beaed4", "#fdc086")
omNames <- rep(c("ref", "skewN", "skewT"), each = 3)
sigNames <- rep(c("low", "med", "high"), length.out = 9)
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
           om = recode(om, "ref" = "Ref.", "skewN" = "Low",
                       "skewT" = "V. Low",
                       .default = levels(om)))
}

q <- ggplot(plotDat, aes(x = om, y = spawners, fill = cu, alpha = sigma)) +
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.75)) +
  geom_hline(plotDat, mapping = aes(yintercept = highBM), linetype = 2) +
  scale_alpha_manual(name = "Component CV", values = c(1, 0.65, 0.3),
                     labels = c("low" = "Low",
                                "med" = "Reference",
                                "high" = "High")) +
                     # labels = c("low" = expression(paste("0.75", sigma)),
                     #            "med" = expression(paste("1.0", sigma)),
                     #            "high" = expression(paste("1.25", sigma)))) +
  scale_fill_manual(values = colPal, guide = FALSE) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  labs(y = "Median Spawner\n Abundance (millions)",
       x = "Productivity Scenario") +
  theme_sleekX(facetSize = 1.2, axisSize = 12, legendSize = 0.85) +
  facet_wrap(~cu, scales = "free_y")


png(file = here("figs/Fig6_spawnerViolinSigma.png"),
    height = 3, width = 5.5, units = "in", res = 300)
print(q)
dev.off()


fullPlotList <- list()
for (i in c(2,5,8)) { #make dataframe
  cuList <- genOutputList(dirNames[i], selectedCUs = selectedCUs, agg = FALSE)
  plotList <- lapply(cuList, function(h) {
    plotDat <- NULL
    nTrials <- nrow(h[["medSpawners"]])
    spwnDat <- data.frame(om = rep(omNames[i], length.out = nTrials * nCUs),
                          synch = rep(unique(h[["opMod"]]), length.out = nTrials * nCUs)
    )
    spwn <- h[["medSpawners"]] %>%
      as.data.frame() %>%
      gather(key = cu, value = spawners)
    spwn$cu <- plyr::revalue(as.factor(spwn$cu), c("V1" = selectedCUs[1],
                                                   "V2" = selectedCUs[2]))
    spwn$lowBM <- rep(h[["meanSGen"]], each = nTrials)
    spwn$highBM <- rep(0.8*h[["meanSMSY"]], each = nTrials)
    plotDat <- rbind(plotDat, cbind(spwnDat, spwn))
    plotDat <- plotDat %>%
      mutate(cu = recode(cu, "Bwrn" = "Bowron", "Chlk" = "Chilko",
                         .default = levels(cu)),
             om = recode(om, "ref" = "Ref.", "skewN" = "Low",
                         "skewT" = "V. Low",
                         .default = levels(om))
             )
    return(plotDat)
  })
  plotDat2 <- do.call(rbind, plotList)
  fullPlotList[[paste0(omNames[i])]] <- plotDat2
}
plotDat3 <- do.call(rbind, fullPlotList) %>%
  mutate(synch =  factor(synch, levels = c("lowSynch", "medSynch", "highSynch")))
row.names(plotDat3) <- NULL

p <- ggplot(plotDat3, aes(x = om, y = spawners, fill = cu, alpha = synch)) +
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.75)) +
  geom_hline(plotDat3, mapping = aes(yintercept = highBM), linetype = 2) +
  scale_alpha_manual(name = "Synchrony", values = c(1, 0.65, 0.3),
                     labels = c("low" = "Low",
                                "med" = "Moderate",
                                "high" = "High")) +
                     # labels = c("lowSynch" = expression(paste(rho, " = 0.05")),
                     #            "medSynch" = expression(paste(rho, " = 0.50")),
                     #            "highSynch" = expression(paste(rho, " = 0.75")))
                     # ) +
  scale_fill_manual(values = colPal, guide = FALSE) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  labs(y = "Median Spawner\n Abundance (millions)",
       x = "Productivity Scenario") +
  theme_sleekX(facetSize = 1.2, axisSize = 12, legendSize = 0.85) +
  facet_wrap(~cu, scales = "free_y")


png(file = paste(here(),"/figs/Fig7_spawnerViolinSynch.png", sep = ""),
    height = 3, width = 5.5, units = "in", res = 300)
print(p)
dev.off()


##### END PRIMARY ANALYSIS #####


#______________________________________________________________________________
## Stand alone analysis looking at relative performance of cyclic vs. non cylic stocks
cuList <- genOutputList("medSig_sockeye", selectedCUs = NULL, agg = FALSE)[["lowSynch_TAM"]]
statusDat <- cuList$ppnYrsLower %>%
  set_colnames(cuList$stkName) %>%
  as.data.frame() %>%
  gather(key = cu, value = status) %>%
  mutate(cu = factor(cu))
larkinCUs <- c("E.St", "L.St", "Qsnl", "L.Sh", "Symr")
statusDat$model <- "ricker"
statusDat[statusDat$cu %in% larkinCUs, ]$model <- "larkin"
statusDat2 <- statusDat %>%
  mutate(model = as.factor(model)) %>%
  arrange(model)

ggplot(statusDat, aes(x = cu, y = status, fill = model)) +
  geom_boxplot() +
  theme_sleekX()

ggplot(statusDat, aes(x = model, y = status, fill = model)) +
  geom_boxplot() +
  theme_sleekX()


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


#_________________________________________________________________________
#### Junky code to make equivalent plots two scenarios

axisSize = 15; dotSize = 3.5; lineSize = 0.8

vars <- c("medRecRY", "ppnCULower", "ppnCUUpper", "ppnCUExtant",
          "medCatch", "ppnFisheriesOpen", "ppnYrsHighCatch", "stabilityCatch")
omNames <- c("ref", "ref", "ref", "skew", "skew", "skew")
# omNames <- c("TAM", "TAM", "TAM", "fixed", "fixed", "fixed")
sigNames <- c("low", "med", "high", "low", "med", "high")
plotDat = NULL
for (h in seq_along(dirNames)) {
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
  plotDat <- rbind(plotDat, singleScen) #merge multiple scenarios into one dataframe
}
plotDat <- plotDat %>%
  mutate(cat = recode(cat, "1" = "low", "2" = "med", "3" = "high", .default = levels(cat)),
         om = recode(om, "ref" = "Reference Productivity", "skew" = "Skewed Productivity",
                     .default = levels(om))
  )
colPal <- viridis(length(levels(plotDat$sigma)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$sigma)

# Plot
consVars <- c("medRecRY", "ppnCULower", "ppnCUUpper", "ppnCUExtant")
consYLabs <- c("Recruit\nAbundance", "Prop. CUs\nLower", "Prop. CUs\nUpper", "Prop. CUs\nExtant")
consPlots <- lapply(seq_along(consVars), function(i) {
  temp <- plotDat %>%
    filter(var == consVars[i])
  q <- ggplot(temp, aes(x = sigma, y = avg, ymin = lowQ, ymax = highQ,
                        color = cat,
                        shape = om
  )) +
    labs(x = "Component Variance", y = consYLabs[i], color = "Sim.\nParameter\nValue") +
    geom_pointrange(fatten = dotSize, size = lineSize, position = position_dodge(width = 0.5)) +
    scale_x_discrete(labels = c("low" = expression(paste("0.75", sigma)),
                                "med" = expression(paste("1.0", sigma)),
                                "high" = expression(paste("1.25", sigma)))) +
    scale_colour_manual(name = "Synchrony", values = colPal,
                        labels = c("low" = expression(paste(rho, " = 0.05")),
                                   "med" = expression(paste(rho, " = 0.50")),
                                   "high" = expression(paste(rho, " = 0.75")))) +
    scale_shape_manual(name = "Productivity", breaks = c("ref", "skew"), values = c(16, 17),
                       guide = FALSE) +
    facet_wrap(~om, scales = "fixed")
  if (i == 1) {
    q <- q + theme_sleekX(position = "top")
  }
  if (i == 2 | i == 3) {
    q <- q + theme_sleekX(position = "mid")
  }
  if (i == 4) {
    q <- q + theme_sleekX(position = "bottom")
  }
  return(q)
})

catchVars <- c("medCatch", "stabilityCatch", "ppnFisheriesOpen", "ppnYrsHighCatch")
catchYLabs <- c("Catch\nAbundance", "Catch Stability", "Prop.\nFisheries Open", "Prop. Years\nHigher Catch")
catchPlots <- lapply(seq_along(catchVars), function(i) {
  temp <- plotDat %>%
    filter(var == catchVars[i])
  q <- ggplot(temp, aes(x = sigma, y = avg, ymin = lowQ, ymax = highQ,
                        color = cat,
                        shape = om
  )) +
    labs(x = "Component Variance", y = catchYLabs[i], color = "Sim.\nParameter\nValue") +
    geom_pointrange(fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.5)) +
    scale_x_discrete(labels = c("low" = expression(paste("0.75", sigma)),
                                "med" = expression(paste("1.0", sigma)),
                                "high" = expression(paste("1.25", sigma)))) +
    scale_colour_manual(name = "Synchrony", values = colPal,
                        labels = c("low" = expression(paste(rho, " = 0.05")),
                                   "med" = expression(paste(rho, " = 0.50")),
                                   "high" = expression(paste(rho, " = 0.75")))) +
    scale_shape_manual(name = "Productivity", breaks = c("ref", "skew"), values = c(16, 17),
                       guide = FALSE) +
    facet_wrap(~om, scales = "fixed")
  if (i == 1) {
    q <- q + theme_sleekX(position = "top")

    # q <- q + theme_sleekX()
  }
  if (i == 2 | i == 3) {
    q <- q + theme_sleekX(position = "mid")
  }
  if (i == 4) {
    q <- q + theme_sleekX(position = "bottom")
  }
  return(q)
})

# pdf(file = paste(here(),"/outputs/synchrony/consGroupedPlots.pdf", sep = ""),
#     height = 6.5, width = 8.5)
png(file = paste(here(),"/outputs/synchrony/consGroupedPlots.png", sep = ""),
    height = 6.5, width = 8.5, units = "in", res = 150)
ggarrange(consPlots[[1]], consPlots[[2]],
          consPlots[[3]], consPlots[[4]],
          ncol = 1, nrow = 4, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1,1.2))
dev.off()
# pdf(file = paste(here(),"/outputs/synchrony/catchGroupedPlots.pdf", sep = ""),
#     height = 6.5, width = 8.5)
png(file = paste(here(),"/outputs/summaryFigs/synchTrials/catchGroupedPlots.png", sep = ""),
    height = 6.5, width = 8.5, units = "in", res = 150)
ggarrange(catchPlots[[1]], catchPlots[[2]],
          catchPlots[[3]], catchPlots[[4]],
          ncol = 1, nrow = 5, common.legend = TRUE, legend = "right",
          align = "v",
          heights = c(1.1,1,1,1.2))
dev.off()

# _____________________________________________________________________________
#
# # Equivalent to above but with synchrony OMs
# plotSynchDat <- NULL
# for (i in c(2, 5)) { #make dataframe
#   cuList <- genOutputList(dirNames[i], selectedCUs = selectedCUs, agg = FALSE)
#   nTrials <- nrow(cuList[[1]][["medSpawners"]])
#   for (h in seq_along(cuList)) {
#     dat <- cuList[[h]]
#     spwnDat <- data.frame(om = rep(omNames[i], length.out = nTrials * nCUs),
#                           synch = rep(dat[["opMod"]], length.out = nTrials * nCUs))
#     spwn <- dat[["medSpawners"]] %>%
#       as.data.frame() %>%
#       gather(key = cu, value = spawners)
#     spwn$cu <- plyr::revalue(as.factor(spwn$cu), c("V1" = selectedCUs[1], "V2" = selectedCUs[2]))
#     plotSynchDat <- rbind(plotSynchDat, cbind(spwnDat, spwn))
#   }
# }
# plotSynchDat$synch = factor(plotSynchDat$synch, levels(plotSynchDat$synch)[c(2, 3, 1)])
#
# spwnSynchHists <- lapply(seq_along(selectedCUs), function(i) {
#   d <- plotSynchDat %>%
#     filter(cu == selectedCUs[i])
#   upp <- subset(bmDat, cu == selectedCUs[i] & bm == "high")[["s"]]
#   low <- subset(bmDat, cu == selectedCUs[i] & bm == "low")[["s"]]
#   q <- ggplot(d, aes(x = spawners, alpha = synch)) +
#     geom_histogram(fill = colPal[i], position = "identity", bins = 40) +
#     facet_wrap(~om) +
#     geom_vline(xintercept = upp, color = "black", linetype = "dashed") +
#     scale_alpha_manual(name = "Synchrony", values = c(0.6, 0.4, 0.2),
#                        labels = c("obs" = "Observed",
#                                   "lowSynch" = expression(paste(rho, " = 0.05")),
#                                   "medSynch" = expression(paste(rho, " = 0.50")),
#                                   "highSynch" = expression(paste(rho, " = 0.75")))) +
#     scale_fill_manual(name = "Conservation\n  Unit", values = colPal,
#                       labels = c("Chilko", "Cultus")) +
#     guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
#     labs(x = "Median Spawner Abundance (millions)", y = "") +
#     geom_text(data = labDat[i ,],
#               mapping = aes(x = Inf, y = Inf, label = lab, hjust = 1.5, vjust = 2),
#               show.legend = FALSE)
#   if (i == 1) {
#     q <- q + theme_sleekX(position = "topWithX")
#   }
#   if (i == 2) {
#     q <- q  + theme_sleekX(position = "bottom")
#   }
#   return(q)
# })
#
# png(file = paste(here(),"/outputs/summaryFigs/synchTrials/spawnerHistsSynch.png", sep = ""),
#     height = 6.5, width = 8.5, units = "in", res = 150)
# p <- ggarrange(spwnSynchHists[[1]], spwnSynchHists[[2]], nrow = 2, ncol = 1, common.legend = TRUE,
#                legend = "right", label.y = "Number of Trials")
# annotate_figure(p,
#                 left = text_grob("Number of Trials", color = "grey30", rot = 90, size = 15))
# dev.off()
