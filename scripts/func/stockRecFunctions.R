## Functions to estimate stock recruit pararameters as inputs in closed-loop 
# simulation model
## Trimmed version of stockRecFunctions.R originally from salmon-sim repo that
# that contains only single stock models for synchrony manuscript.

library("dplyr"); library("parallel"); library("doParallel"); library("foreach")


## Single stock Bayesian Ricker or Larkin models
# dat=srDatTrim
# Niter=12500; Nthin=50; burnin=2500
# fname="TEST"; savedCores=1; model = "larkin"
# cap = 1; betaPrior = "logN"; tauPrior = 0.001
RBayesSingle <- function(dat = SRdat, Niter = 75000, Nthin = 50, burnin = 25000,
                         fname, model = "ricker", betaPrior = "logN", 
                         cap = 3, tauPrior = 0.001) {
  
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/figs/", sep="")), 
         dir.create(paste(here(),"/outputs/srAnalysis/figs/", sep="")), FALSE)
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/data/", sep="")), 
         dir.create(paste(here(),"/outputs/srAnalysis/data/", sep="")), FALSE)
  
  # prepare data to feed into model
  # number of CUs
  stkNames <- unique(dat$stkName)
  CUs <- unique(dat$CU)
  nCUs <- length(CUs)
  
  inits <- function(){
    list(tau_R = 1, alpha=rep(1, nsites))
  }
  
  #indices of non NA data points for R/S and for just S
  inds <- which(is.na(dat$lnRS)==F & is.na(dat$Escape)==F)
  inds2 <- which(is.na(dat$Escape)==F)
  # Recruits and spawners data
  dataGood <- dat[inds,] 
  lnRSdat <- dataGood$lnRS
  Sdat <- dataGood$Escape
  escGood <- dat[inds2,] 
  # Years and CUs corresponding to data points, as indices
  # need to change so that missing years are accomodated, differing numbers of years for CUs 
  sampledCUs <- NULL
  for(i in 1:length(inds)){
    sampledCUs[i] <- which(CUs==dataGood$CU[i])
  }
  # need upper bound for Smax prior for each site, number of years each site
  upp <- NULL
  Years <- NULL
  nyears <- NULL
  MuSmax <- NULL
  TauSmax <- NULL
  for(i in 1:length(CUs)){
    #use data frame from ricker fit to use later to initialize model
    cuDat <- dataGood[which(dataGood$CU==CUs[i]), ]
    EscDat <- escGood[which(escGood$CU==CUs[i]), ]
    #calculate upper bound on Smax prior
    upp[i] <- max(EscDat$Escape) #
    nyears[i] <- dim(cuDat)[1]
    Years <- c(Years, 1:nyears[i])
  }    
  
  # convert obs spawners and recruits to matrices
  Smat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  lnRSmat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  for(k in 1:length(inds)){
    Smat[sampledCUs[k], Years[k]] <- Sdat[k]
    lnRSmat[sampledCUs[k], Years[k]] <- lnRSdat[k]
  }
  
  # create data input for JAGS
  Dlist <- list("Smat" = Smat, "lnRSmat" = lnRSmat, "nyears" = nyears, 
                "upp" = upp, "nCUs" = nCUs, "cap" = cap, "tauPrior" =  tauPrior)
  if (model == "larkin") {
    if (betaPrior != "unif") {
      simoutJAGS <- jags(Dlist,
                         model.file = here("scripts", "jagsModels", 
                                           "singleStockLarkin_logN.txt"),
                         parameters.to.save=c("alpha", "beta", "beta1", "beta2", 
                                              "beta3", "tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    } else {
      simoutJAGS <- jags(Dlist,
                         model.file = here("scripts", "jagsModels", 
                                           "singleStockLarkin_unif.txt"),
                         parameters.to.save=c("alpha", "beta", "beta1", "beta2", 
                                              "beta3","tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    }
  } #end if model = Larkin
  if (model == "ricker") {
    if (betaPrior != "unif") {
      simoutJAGS <- jags(Dlist,
                         model.file = here("scripts", "jagsModels", 
                                         "singleStockRicker_logN.txt"),
                         parameters.to.save=c("alpha", "beta", "tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    } else {
      simoutJAGS <- jags(Dlist,
                         model.file = here("scripts", "jagsModels", 
                                           "singleStockRicker_unif.txt"),
                         parameters.to.save=c("alpha", "beta", "tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    }
  } #end if model = Rciker
  # change to bugs output to easily extract values
  simoutBUGS <- simoutJAGS$BUGSoutput
  
  # extract mcmc list to plot usuing coda functions
  ToPlot <- as.mcmc(simoutJAGS)
  # plot traces and posterior densities
  pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"TraceDens.pdf", sep="") )
  plot(ToPlot)
  dev.off()
  # plot gelman plots
  pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"Gelman.pdf", sep="" ))
  tryCatch( {
    gelman.plot(ToPlot)
  }, error=function(x=1){plot.new() })
  dev.off()
  
  # Re-arrange and store each iteration's parameter estimates
  if(nrow(simoutBUGS$sims.matrix) > 1000){
    simSamples <- sample_n(as.data.frame(simoutBUGS$sims.matrix), size=1000) #subset to 1000 samples equivalent to Fraser outputs generated by FRSSI modelers
  } else{
    simSamples <- simoutBUGS$sims.matrix
  }
  #Change header for simSamples since default is different when nCU=1
  if (nCUs == 1) {
    colnames(simSamples) <- paste(colnames(simSamples), "[1]", sep="") 
    simSamples <- simSamples %>% 
      dplyr::rename(deviance = "deviance[1]")
  }
  simPars <- NULL
  for(i in seq_along(CUs)){
    CU <- rep(CUs[i], each=nrow(simSamples))
    alpha <- simSamples[, paste('alpha', "[", i , "]", sep="")]
    beta0 <- simSamples[,paste('beta', "[", i , "]", sep="")]
    sigma <- simSamples[,paste('sig_R', "[", i , "]", sep="")]
    deviance <- simSamples[, "deviance"]
    allParsTemp <- cbind(CU, alpha, beta0, sigma, deviance)
    if (model == "larkin") {
      beta1 <- simSamples[,paste('beta1', "[", i , "]", sep="")]
      beta2 <- simSamples[,paste('beta2', "[", i , "]", sep="")]
      beta3 <- simSamples[,paste('beta3', "[", i , "]", sep="")]
      allParsTemp <- cbind(allParsTemp, beta1, beta2, beta3)
    }
    allParsTemp <- allParsTemp %>% 
      as.data.frame() %>% 
      mutate(stk = stkNames[i])
    simPars <- rbind(simPars, allParsTemp)
  }
  
  #save output for later
  saveRDS(simoutJAGS, file=paste(here(),"/outputs/srAnalysis/data/",
                                 fname,"SimOutput",".rds", sep=""))
  write.csv(simPars, file=paste(here(),"/outputs/srAnalysis/data/",
                                fname,"MCMCPars.csv", sep=""), row.names=FALSE)
}