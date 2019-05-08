#*******************************************************************************
# frStockRecModels_synchDir.r
# Date last revised: 4 May, 2019
# Copy of script (frStockRecModel.R) moved from salmon-sim directory. Trimmed to
# remove preliminary analyses. 
#*******************************************************************************

require(here); require(R2jags); require(ggplot2); require(rstan); require(dplyr)

source(here("scripts/func/stockRecFunctions.R"))

# Clean data to match format for stockRecFunctions.R
# Use total
srDat <- read.csv(here("data/sox/fraserRecDat.csv"), 
                  stringsAsFactors=F) %>% 
  select(stk, stkName, yr, totalSpwn, rec) %>% 
  rename(CU = stk, Year = yr, Escape = totalSpwn, Recruit = rec)

### Fit all CUs with full dataset; exception is Cultus which has strong hatchery
## impacts since 2000
srDatTrim <- srDat %>% 
  filter(!(stkName == "Cultus" & Year > 1999)) %>% 
  mutate(lnRS = log(Recruit/Escape))

RBayesSingle(dat = srDatTrim, Niter=250000, Nthin=50, burnin=50000, 
             fname="fraserPooledRicker", model = "ricker")
RBayesSingle(dat = srDatTrim, Niter=250000, Nthin=50, burnin=50000, 
             fname="fraserPooledLarkin", model = "larkin")


modelSummRick <- readRDS(here("outputs/srAnalysis/data", 
                              paste("fraserPooledRicker", "SimOutput.RDS", 
                                    sep = "")))
modelSummLark <- readRDS(here("outputs/srAnalysis/data", 
                              paste("fraserPooledLarkin", "SimOutput.RDS", 
                                    sep = "")))

pooledRicker <- read.csv(here("outputs/srAnalysis/data", 
                              paste("fraserPooledRicker", "MCMCPars.csv", 
                                    sep = ""))) %>% 
  select(CU, alpha, beta0, sigma, deviance) %>% 
  mutate(deviance = -1 * deviance) %>% 
  rename(stk = CU)

pooledLarkin <- read.csv(here("outputs/srAnalysis/data", 
                          paste("fraserPooledLarkin", "MCMCPars.csv", 
                                sep = ""))) %>% 
  select(CU, alpha, beta0, beta1, beta2, beta3, sigma, deviance) %>% 
  mutate(deviance = -1 * deviance) %>% 
  rename(stk = CU)


## Look at medians to check if they seem reasonable
pooledRicker %>% 
  group_by(stk) %>% 
  summarize(medA = median(alpha),
            medB = median(beta0))
pooledLarkin %>% 
  group_by(stk) %>% 
  summarize(medA = median(alpha),
            medB = median(beta0))

write.csv(pooledRicker,
          here("data", "sox", "pooledRickerMCMCPars.csv"),
          row.names = FALSE)
write.csv(pooledLarkin,
          here("data", "sox", "pooledLarkinMCMCPars.csv"),
          row.names = FALSE)