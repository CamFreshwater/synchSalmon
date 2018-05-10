#*************************************************************************************
# prepCovariates.R
# Date revised: April 20, 2018
# Inputs: environmental and catch data 
# Outputs: single clean csv to be passed to synch models
#*************************************************************************************


# setwd("C:/github/synchSalmon/")
setwd("/Users/cam/github/synchSalmon") #Cam's Mac wd

require(here); require(dplyr); require(tidyr)

source(here("scripts/synchFunctions.R"))


## Import covariates
sstPca <- read.table(here("/data/env/sstPCA.txt")) #principal components of SST variation in NE Pacific (170E to 240E, 40N-65N)
sstRaw <- read.table(here("/data/env/sstRaw.txt")) # pacific ocean SST
pdo <- read.csv(here("/data/env/pdo.csv"), stringsAsFactors=F) 
npgo <- read.csv(here("/data/env/npgo.csv"), stringsAsFactors=F)
meanAlpi <- read.csv(here("/data/env/alpi.csv"), stringsAsFactors=F) #stops at 2015
catchDat <- read.csv(here("/data/env/cleanCatch.csv"), stringsAsFactors=F)


#_________________________________________________________________________
months <- c("3","4","5","6")

colnames(sstPca) <- c("id", "year", "month", "pc1", "pc2", "pc3", "pc4", "pc5")
primePca <- sstPca[sstPca$month %in% months, c("year", "month", "pc2")]
pc2Prime <- primePca %>%
			group_by(year) %>%
			summarise(pc2Prime = mean(pc2))
pc2 <- sstPca %>%
		group_by(year) %>%
		summarise(pc2 = mean(pc2))

colnames(sstRaw) <- c("long", "lat", "year", "month", "temp")
sstPrime <- sstRaw[sstRaw$month %in% months, c("year","month","temp")] 
sstPrime <- sstPrime %>%
			group_by(year) %>%
			summarise(sstPrime = mean(temp))
sst <- sstRaw %>%
		group_by(year) %>%
		summarise(sst = mean(temp))

pdo <- data.frame(year = pdo$YEAR,
				  pdoPrime = apply(pdo[, c("MAR","APR","MAY","JUN")], 1, mean),
				  pdo = apply(pdo[, 2:13], 1, mean)
				  )

npgoPrime <- npgo[npgo$month %in% months, c("year", "month", "index")] 
npgoPrime <- npgoPrime %>%
			group_by(year) %>%
			summarise(npgoPrime = mean(index))
npgo <- npgo %>%
		group_by(year) %>%
		summarise(npgo = mean(index))

catchTrim <- subset(catchDat, region == "Alaska")
catch <- subset(catchTrim, fishery == "Commercial") %>%
			group_by(year) %>%
			select(year, species, totalCatch)
catch <- reshape2::dcast(catch, year ~ species)
names(catch)[2:3] <- c("pinkCatch", "soxCatch")

covarDat <- Reduce(function(x, y) merge(x, y, by = c("year")), list(pc2, pc2Prime, sst, sstPrime, pdo, npgo, npgoPrime, catch))

write.csv(covarDat, here("/data/env/cleanCovar.csv"), row.names = FALSE)