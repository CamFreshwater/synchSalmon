#*******************************************************************************
# runModelRobustOM_synchSensitivity.R
# Date revised: ONGOING
# Inputs: recoverySim.R
# Outputs: pdf plots
# Explainer: Equivalent to runModelRobustOM_synch.R, but focused on sensitivity 
# analyses so smaller range of OMs examined (i.e. only median levels sigma and 
# correlations); plots are stripped down as a result
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

simPar <- read.csv(here("data/sox/fraserOMInputs_varyCorrSensitivity.csv"), 
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

## Check SD among CU-specific uncertainty parameters to ensure that proposed
# sensitivity values are reasonable
sdER <- sd(cuPar$meanDBE)
par(mfrow = c(2, 2))
hist(cuPar$sdDBE * 1.25)
hist(cuPar$sdDBE * 0.75)
hist(cuPar$sdDBE + sdER)
hist(cuPar$sdDBE - sdER) #lines up well

sdTau <- sd(cuPar$tauCycAge)
hist(cuPar$tauCycAge * 1.25)
hist(cuPar$tauCycAge * 0.5)
hist(cuPar$tauCycAge + sdTau)
hist(cuPar$tauCycAge - sdTau) #lines up well

### SET UP MODEL RUN -----------------------------------------------------

## Define simulations to be run
nTrials <- 1500

## General robustness runs
simParTrim <- simPar
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))
 
# recoverySim(simParTrim[1, ], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             cuCustomCorrMat = cuCustomCorrMat, dirName = "test", nTrials = 5,
#             multipleMPs = FALSE)

for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 2) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MASS),
                     library(here),
                     library(sensitivity),
                     library(mvtnorm),
                     library(scales), #shaded colors for figs
                     library(here), #use this package so wd are common across different computers
                     library(synchrony),
                     library(zoo), #synch and zoo used to calculate rolling estimates of synchrony
                     library(viridis), #color blind gradient palette
                     library(ggplot2),
                     library(gsl), #to calculate exact estimate of MSY following Scheuerell 2016 PeerJ
                     library(dplyr),
                     library(Rcpp),
                     library(RcppArmadillo),
                     library(sn)))
  clusterExport(cl, c("simsToRun", "recoverySim", "cuPar", "dirName", "nTrials",
                      "catchDat", "srDat",
                      "ricPars", "dirName", "larkPars", "tamFRP"), 
                envir = environment()) #export custom function and objects
  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    recoverySim(x, cuPar, catchDat = catchDat, srDat = srDat, variableCU = FALSE, 
                ricPars, larkPars = larkPars, tamFRP = tamFRP,  
                dirName = dirName, nTrials = nTrials, makeSubDirs = FALSE, 
                random = FALSE)
  })
  stopCluster(cl) #end cluster
  toc()
}


#_________________________________________________________________________
## Modified version of multi-OM grouped box plots
vars <- c("medRecRY", "ppnCUUpper", "ppnMixedOpen",
          "medCatch", "ppnYrsHighCatch", "stabilityCatch")

plotDat = NULL
for (h in seq_along(dirNames)) {
  agList <- genOutputList(dirNames[h], agg = TRUE)
  singleScen = NULL
  for (i in seq_along(vars)) {
    dum <- data.frame(scen = as.factor(rep(scenNames[h], 
                                           length.out =  length(agList))),
                      var = rep(vars[i], length.out = length(agList)),
                      om = as.factor(sapply(agList, function(x) unique(x$opMod))),
                      avg = sapply(agList, function(x) median(x[,vars[i]])),
                      lowQ = sapply(agList, function(x) qLow(x[,vars[i]])),
                      highQ = sapply(agList, function(x) qHigh(x[,vars[i]])),
                      row.names = NULL
    )
    singleScen <- rbind(singleScen, dum)
  }
  #merge multiple scenarios into one dataframe
  plotDat <- rbind(plotDat, singleScen)
}
plotDat <- plotDat  %>% 
  mutate(om =  factor(om, levels = c("reference", "low", "high")),
         scen = factor(scen, levels = c("rho", "ageTau", "enRouteSig", "ouSig", 
                                        "refSensitivity")))


colPal <- c("black", "blue", "orange")
names(colPal) <- levels(plotDat$om)

dotSize = 3; lineSize = 0.8; legSize = 0.8; axSize = 10
consVars <- c("medRecRY", "ppnCUUpper", "ppnMixedOpen") 
consYLabs <- c("Return\nAbundance", "Prop. MUs\nAbove Esc. Goal",
               "Prop. CUs\nAbove Benchmark")
consPlots <- lapply(seq_along(consVars), function(i) {
  temp <- plotDat %>% 
    filter(var == consVars[i],
           !scen == "refSensitivity")
  refDat <- plotDat %>% 
    filter(var == consVars[i],
           scen == "refSensitivity")
  q <- ggplot(temp, aes(x = scen, y = avg, ymin = lowQ, ymax = highQ, 
                        color = om)) +
    labs(x = "Scenario", y = consYLabs[i], color = "Operating\nModel") +
    geom_pointrange(fatten = dotSize, size = lineSize, position = 
                      position_dodge(width = 0.65)) +
    scale_x_discrete(labels = c("refSensitivity" = "Reference",
                                "ageTau" = "Maturation\nAge", "ouSig" = 
                                  "Outcome\nUncertainty", "rho" = 
                                  "Temporal\nAutocorrelation", 
                                "enRouteSig" = "En Route\nMortality")) +
    scale_colour_manual(values = colPal) +
    geom_hline(refDat, mapping = aes(yintercept = avg), linetype = 1) +
    geom_hline(refDat, mapping = aes(yintercept = lowQ), linetype = 2, 
               size = 0.5) +
    geom_hline(refDat, mapping = aes(yintercept = highQ), linetype = 2, 
               size = 0.5) 
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = legSize,
                          axisSize = axSize) 
  } 
  if (i == 2) {
    q <- q + theme_sleekX(position = "mid", legendSize = legSize,
                          axisSize = axSize)
  }
  if (i == 3) {
    q <- q + theme_sleekX(position = "bottom", legendSize = legSize,
                          axisSize = axSize)
  }
  return(q)
})

catchVars <- c("medCatch", "ppnYrsHighCatch", "stabilityCatch")
catchYLabs <- c("Catch\nAbundance", "Prop. Years\nHigher Catch",
                "Catch Stability")
catchPlots <- lapply(seq_along(catchVars), function(i) {
  temp <- plotDat %>% 
    filter(var == catchVars[i],
           !scen == "refSensitivity")
  refDat <- plotDat %>% 
    filter(var == catchVars[i],
           scen == "refSensitivity")
  q <- ggplot(temp, aes(x = scen, y = avg, ymin = lowQ, ymax = highQ, 
                        color = om)) +
    labs(x = "Scenario", y = catchYLabs[i], color = "Operating\nModel") +
    geom_pointrange(fatten = dotSize, size = lineSize, position = 
                      position_dodge(width = 0.65)) +
    scale_x_discrete(labels = c("refSensitivity" = "Reference",
                                "ageTau" = "Maturation\nAge", "ouSig" = 
                                  "Outcome\nUncertainty", "rho" = 
                                  "Temporal\nAutocorrelation", "enRouteSig" = 
                                  "En Route\nMortality")) +
    scale_colour_manual(values = colPal) +
    geom_hline(refDat, mapping = aes(yintercept = avg), linetype = 1) +
    geom_hline(refDat, mapping = aes(yintercept = lowQ), linetype = 2, 
               size = 0.5) +
    geom_hline(refDat, mapping = aes(yintercept = highQ), linetype = 2, 
               size = 0.5)
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = legSize, 
                          axisSize = axSize)
  }
  if (i == 2) {
    q <- q + theme_sleekX(position = "mid", legendSize = legSize, 
                          axisSize = axSize)
  }
  if (i == 3) {
    q <- q + theme_sleekX(position = "bottom", legendSize = legSize, 
                          axisSize = axSize)
  }
  return(q)
})

png(file = paste(here("/figs/SFig_ConsPMs.png"),
                 sep = ""), 
    height = 4.5, width = 6, units = "in", res = 600)
ggarrange(consPlots[[1]], consPlots[[2]], consPlots[[3]], 
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right", 
          align = "v", heights = c(1,1,1.25))
dev.off()
png(file = paste(here("/figs/SFig_CatchPMs.png"),
                 sep = ""), 
    height = 4.5, width = 6, units = "in", res = 600)
ggarrange(catchPlots[[1]], catchPlots[[2]], catchPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right", 
          align = "v", heights = c(1,1,1.25))
dev.off()

