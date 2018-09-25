### Dummy data set testing effect of weighting variable on outcome
require(here); require(synchrony); require(zoo); require(ggplot2); require(dplyr)
require(tidyr); require(viridis); require(ggpubr); require(MASS)

source(here("scripts/synchFunctions.R"))

### Generate productivity based on Fraser 
mu <- c(1.3522896, 1.683566, 1.844185, 1.6001229, 1.514986)
sig <- rnorm(n = 5, 0.7, 0.2)
sigMat <- matrix(as.numeric(sig), nrow = 1, ncol = 5) #calculate correlations among CUs
covMat <- t(sigMat) %*% sigMat #calculate shared variance
corMat <- covMat * 0.5 #correct based on correlation
diag(corMat) <- as.numeric(sig^2) #add variance
prodMat <- mvrnorm(n = 40, mu, corMat, tol = 1e-6)
agPMat <- apply(prodMat, 1, sum)

S <- c(0.5, 0.05, 0.1, 0.3, 1)
sig2 <- rnorm(n = 5, 0.5, 0.2)
sigMat <- matrix(as.numeric(sig2), nrow = 1, ncol = 5) #calculate correlations among CUs
covMat <- t(sigMat) %*% sigMat #calculate shared variance
corSMat <- covMat * 0.2 #correct based on correlation
diag(corSMat) <- as.numeric(sig2^2) #add variance
sMat <- mvrnorm(n = 40, S, corSMat, tol = 1e-6)
sMat[which(sMat <= 0)] <- 0.0001
agSMat <- apply(sMat, 1, sum)

trueAgCV = rollapplyr(agPMat, width = 12, function(x) sd(x) / mean(x))
wCV = rollapplyr(prodMat, width = 12,  function(x) wtdCV(x, weightMat = NULL),
                   fill = NA, by.column = FALSE)
synch = rollapplyr(prodMat, width = 12, function(x) community.sync(x)$obs, 
                   fill = NA, by.column = FALSE)
agCV = rollapplyr(prodMat, width = 12, function(x) cvAgg(x, weightMat = NULL),
                  fill = NA, by.column = FALSE)

wCV2 = rollapplyr(prodMat, width = 12,  function(x) wtdCV(x, weightMat = sMat),
                   fill = NA, by.column = FALSE)
agCV2 = rollapplyr(prodMat, width = 12, function(x) cvAgg(x, weightMat = sMat),
                  fill = NA, by.column = FALSE)
agCV <- agCV[!is.na(agCV)]
agCV2 <- agCV[!is.na(agCV2)]
plot(trueAgCV ~ agCV, pch = 21, bg = "black")
points(trueAgCV ~ agCV2, pch = 21, bg = "white")
