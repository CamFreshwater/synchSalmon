#*************************************************************************************
# synchFunctions.R
# Date revised: May 9, 2018; ONGOING
# Outputs: functions to be used in salmon synchrony analysis
#*************************************************************************************


#_________________________________________________________________________
# This function calculates weighted CV following eq. 4 in Thibault and Connolly 
# 2013 Ecol Let
wtdCV <- function(recMat){ 
  aggAbund <- sum(apply(recMat, 2, mean)) #temporal mean of aggregate abundance
  sum(apply(recMat, 2, function(x) (mean(x) / aggAbund) * (sqrt(var(x)) / mean(x))))
}


#_________________________________________________________________________
# This function calculates aggregate CV following eq. 3 in Thibault and Connolly 
# 2013 Ecol Let
cvAgg <- function(recMat){ 
  sqrt(community.sync(recMat)$obs) * wtdCV(recMat)
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

