#*************************************************************************************
# synchFunctions.R
# Date revised: May 9, 2018; ONGOING
# Outputs: functions to be used in salmon synchrony analysis
#*************************************************************************************


#_________________________________________________________________________
# This function calculates weighted CV following eq. 4 in Thibault and Connolly 
# 2013 Ecol Let; recMat necessary to weight values appropriately
wtdCV <- function(datMat, recMat = NULL){
  if (is.null(recMat)) { #if recMat is NULL assume datMat is a matrix of abundance
    recMat <- datMat
  }
  if (ncol(datMat) != ncol(recMat)){
    stop("Input matrices have unequal number of components")
  }
  aggAbund <- sum(apply(recMat, 2, mean)) #temporal mean of aggregate abundance
  wtdAbund <- apply(recMat, 2, function (x) mean(x) / aggAbund) #wtd mean of aggregate abundance
  wtdCV <- sum(wtdAbund * apply(datMat, 2, function(x) sqrt(var(x))/ mean(x)))
  return(wtdCV)
}


#_________________________________________________________________________
# This function calculates aggregate CV following eq. 3 in Thibault and Connolly 
# 2013 Ecol Let
cvAgg <- function(datMat, recMat = NULL){
  if (is.null(recMat)) {
    recMat <- datMat
  }
  sqrt(community.sync(datMat)$obs) * wtdCV(datMat, recMat)
}


#_________________________________________________________________________
# This function calculates a mean weighted by abundance
wtdMean <- function(datMat, recMat = NULL){
  if (is.null(recMat)) { #if recMat is NULL assume datMat is a matrix of abundance
    recMat <- datMat
  }
  if (ncol(datMat) != ncol(recMat)){
    stop("Input matrices have unequal number of components")
  }
  aggAbund <- sum(apply(recMat, 2, mean)) #temporal mean of aggregate abundance
  wtdAbund <- apply(recMat, 2, function (x) mean(x) / aggAbund) #wtd mean of aggregate abundance
  wtdMean <- sum(wtdAbund * apply(datMat, 2, function(x) mean(x)))
  return(wtdMean)
}


#_________________________________________________________________________
# These functions calculate CVcurrent (eq. 1) or CVnull (eq. 2) used to calculate
# diversity deficit (DD; eq. 3) as in Yamane et al. 2018 J App Ecol
calcCV <- function(recMat, current = TRUE) {
  sumVar <- sum(apply(recMat, 2, var))
  sumCov <- calcCov(recMat)
  sumMean <- sum(apply(recMat, 2, mean))
  if (current == TRUE) {
    cv <- sqrt(sumVar + sumCov) / sumMean
  } 
  if (current == FALSE) {
    cv <- sqrt(sumVar) / sumMean
  }
  return(cv)
}

calcCov <- function(recMat) {
  n <- ncol(recMat)
  covMat <- matrix(nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (h in (i + 1):n) {
      covMat[i, h] <- cov(recMat[ , i], recMat[ , h])
    }
  }
  sumCov <- sum(covMat, na.rm = TRUE)
  return(sumCov)
}

