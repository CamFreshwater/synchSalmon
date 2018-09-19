#*************************************************************************************
# synchFunctions.R
# Date revised: May 9, 2018; ONGOING
# Outputs: functions to be used in salmon synchrony analysis
#*************************************************************************************


#_________________________________________________________________________
# This function calculates weighted CV following eq. 4 in Thibault and Connolly 
# 2013 Ecol Let; recMat necessary to weight values appropriately
wtdCV <- function(datMat, wtMat = NULL, weight = TRUE) {
  if (is.null(wtMat)) { #if recMat is NULL assume datMat is a matrix of abundance
    wtMat <- datMat
  }
  if (ncol(datMat) != ncol(wtMat)){
    stop("Input matrices have unequal number of components")
  }
  if (any(is.na(datMat))) {
    warning("NAs present. This will affect estimates of weighted CV.")
  }
  #temporal mean of aggregate abundance
  aggAbund <- sum(apply(wtMat, 2, function (x) mean(x, na.rm = TRUE)))
  #wtd mean of aggregate abundance
  wtdAbund <- apply(wtMat, 2, function (x) mean(x, na.rm = TRUE) / aggAbund)
  wtdCV <- sum(wtdAbund * apply(datMat, 2, 
                                function(x) sqrt(var(x, na.rm = TRUE)) / 
                                  mean(x, na.rm = TRUE)), na.rm = TRUE)
  unWtdCV <- mean(apply(datMat, 2, 
                        function(x) sqrt(var(x, na.rm = TRUE)) / 
                          mean(x, na.rm = TRUE)), na.rm = TRUE)
  if (weight == TRUE) {
    return(wtdCV) 
  } else {
      return(unWtdCV)
    }
}

#_________________________________________________________________________
# This function calculates aggregate CV following eq. 3 in Thibault and Connolly 
# 2013 Ecol Let
cvAgg <- function(datMat, wtMat = NULL, weight = TRUE){
  if (is.null(wtMat)) {
    wtMat <- datMat
  }
  sqrt(community.sync(datMat)$obs) * wtdCV(datMat, wtMat, weight = weight)
}


#_________________________________________________________________________
# This function calculates a mean weighted by abundance
wtdMean <- function(datMat, wtMat = NULL, weight = TRUE){
  if (is.null(wtMat)) { #if wtMat is NULL assume datMat is a matrix of abundance
    wtMat <- datMat
  }
  if (ncol(datMat) != ncol(wtMat)){
    stop("Input matrices have unequal number of components")
  }
  aggAbund <- sum(apply(wtMat, 2, mean)) #temporal mean of aggregate abundance
  wtdAbund <- apply(wtMat, 2, function (x) mean(x) / aggAbund) #wtd mean of aggregate abundance
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


#___________________________________________________________________________________________________________
# This function is a copy of S. Anderson's now that ggsidekick is not supported; adds argument
# for top, bottom, middle for multipanel functionality
theme_sleekX <- function(base_size = 11, base_family = "", position = "standard") {
  half_line <- base_size/2
  q <- theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(1.1)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = 1),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
  if (position == "bottom") {
    q <- q + theme(strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   strip.text.y = element_blank(),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_text(size = 0.85 * axisSize),
                   axis.title = element_text(size = axisSize)
    )
  }
  if (position == "top") {
    q <- q + theme(strip.text = element_text(size = axisSize),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_blank(),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  if (position == "topWithX") {
    q <- q + theme(strip.text = element_text(size = axisSize),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_text(size = 0.85 * axisSize),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  if (position == "mid") {
    q <- q + theme(strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   strip.text.y = element_blank(),
                   axis.text.y = element_text(size = 0.85 * axisSize),
                   axis.text.x = element_blank(),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  return(q)
}

