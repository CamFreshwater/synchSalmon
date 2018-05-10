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
