require(sn)

## Sean's exploration into centered parameterization
cp <- list(mean=c(0,0), var.cov=matrix(c(3,2,2,3)/3, 2, 2), gamma1=c(0.8, 0.4))
dp <- cp2dp(cp, "SN")
rnd <- rmsn(1e4, dp=dp)
mean(rnd)

rnd <- rmsn(1e4, xi = c(0, 0), Omega = matrix(c(3,2,2,3)/3, 2, 2), 
            alpha = c(0.8, 0.4))
mean(rnd)



## Original examples from Azzalini pres
data(wines)

olo.ph <- wines[wines$wine=="Barolo", "phenols"]
fitLM <- lm(olo.ph ~ 1)
fit <- selm(olo.ph ~ 1, family="SN")
plot(fit, which=2:3)
#
# try
summary(fit) # works with CP
summary(fit, param.type="DP")


fit2 <- selm(cbind(tartaric, malic, uronic) ~ colour + hue,
             family="ST", data=wines, subset=(wine != "Grignolino"))
plot(fit2, which=2:3)
summary(fit2) # works with CP
summary(fit2, param.type="DP")
# constraint on degrees of freedom:
fit3 <- selm(cbind(tartaric, malic, uronic) ~ colour + hue,
             family="ST", fixed.param=list(nu=8),
             data=wines, subset=(wine != "Grignolino"))


## Play around with variance covariance matrices based on residuals
require(tidyverse)

cuPars <- read.csv(here("data", "sox", "fraserCUPars.csv"))
residMat <- readRDS(here("outputs", "generatedData", "residMat.rds"))
colnames(residMat) <- cuPars$stkName
nStks <- ncol(residMat)

calcCorMat <- function(mat) {
  N <- ncol(mat)
  corMat <- matrix(NA, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      dum <- cbind(mat[ , i], mat[ , j])
      dum <- dum[complete.cases(dum), ]
      corMat[j, i] <- cor(dum)[2,1]
    }
  }
  return(corMat)
}

## Generate variance-covariance based on residuals
# estimate CU-specific variance
sigmas <- apply(residMat, 2, function (x) sd(x, na.rm = T))
corMat <- calcCorMat(residMat)
sigMat <- matrix(as.numeric(sigmas), nrow = 1, ncol = nStks)
covMat <- (t(sigMat) %*% sigMat) * corMat
diag(covMat) <- sigmas^2

## Compare means for full and later periods
mean(corMat[lower.tri(corMat)])

trim <- residMat[(nrow(residMat) - 25):nrow(residMat), ]
sigmasS <- apply(trim, 2, function (x) sd(x, na.rm = T))
corMatS <- calcCorMat(trim)
mean(corMatS[lower.tri(corMatS)])


skewFit <- selm(residMat ~ 1, family="SN")
skewFitTrim <- selm(trim ~ 1, family="SN")


## Look at histograms 
residDat <- residMat %>% 
  as.data.frame() %>%
  mutate(yr = row.names(residMat),
         dataset = "full") %>% 
  gather(key = stk, value = resid, -yr, -dataset)
recResidDat <- residDat %>% 
  filter(yr > 1990) %>% 
  mutate(dataset = "recent")

plotDat <- rbind(residDat, recResidDat) %>% 
  mutate(dataset = as.factor(dataset))

colPal <- viridis(length(levels(plotDat$dataset)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$dataset)

ggplot(plotDat, aes(x = resid, fill = dataset)) +
  geom_histogram() +
  facet_wrap(~stk)




f <- function(N = 15, r = 0.4, sigma = 0.2, skew = 1) {
  sig_mat <- matrix(as.numeric(sigma), nrow = 1, ncol = N)
  cov_mat <- (t(sig_mat) %*% sig_mat) * r
  diag(cov_mat) <- sigma^2
  x <- sn::rmsn(1e2, xi = rep(0, N), Omega = cov_mat,
                alpha = rep(log(skew), N))
  apply(x, 1, function(xx) sum(exp(xx)))
}

.r <- rep(seq(0, 0.9, 0.1), each = 200)
y <- purrr::map_df(.r, function(.x) {
  out <- f(r = .x)
  tibble::tibble(
    .median = median(out),
    .mean = mean(out)
  )
})
y$correlation <- as.factor(.r)

library(ggplot2)
library(dplyr)

group_by(y, correlation) %>%
  ggplot(aes(correlation, .median)) +
  geom_boxplot()

group_by(y, correlation) %>%
  ggplot(aes(correlation, .mean)) +
  geom_boxplot()
