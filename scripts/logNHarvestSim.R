## Quick simulation to look at effects of different parameterizations for log-
# normal uncertainty


exTAC <- c(0.0221943066, 0.0936264436, 0.0548203227, 0.0001836434, 
           0.0001376739, 0.7077153698, 0.2844746972, 0.0064115127, 0.0974588335, 
           0.0281936863, 0.0000000000, 0.0044093000, 0.0065105755, 0.0008444005,
           0.0008352194, 0.0007538029, 0.0004385476, 0.0007639658, 0.0051275717)
exLNErr <- exp(qnorm(runif(19, 0.0001, 0.9999), 0, 0.2))
exNErr <- qnorm(runif(19, 0.0001, 0.9999), 0, 0.2)

norm <- NULL
logNorm <- NULL
logNorm2 <- NULL
for(i in 1:9999) {
  C1 <- exTAC * exOutErr
  C2 <- exp(log(exTAC) * exOutErr)
  C3 <- exp(log(exTAC)) * exp(exOutErr)
  norm <- c(norm, C1)
  logNorm <- c(logNorm, C2)
  logNorm2 <- c(logNorm2, C3)
}

par(mfrow = c(3, 1))
hist(norm)
hist(logNorm)
hist(logNorm2)

set.seed(123)
hr1 <- 0.2*exp(qnorm(runif(1000, 0.0001, 0.9999), 0, 0.2))
set.seed(123)
hr2 <- -log(1-0.2)*exp(qnorm(runif(1000, 0.0001, 0.9999), 0, 0.2))
set.seed(123)
hr3 <-0.2 + qnorm(runif(1000, 0.0001, 0.9999), 0, 0.2)

par(mfrow = c(3, 1))
hist(hr1)
hist(hr2)
hist(hr3)

c1a = hr1 * 2000
set.seed(123)
c1b = (0.2 * 2000) * exp(qnorm(runif(1000, 0.0001, 0.9999), 0, 0.2))
c2a = hr2 * 2000
set.seed(123)
c2b = -log(1-(400/2000))*exp(qnorm(runif(1000, 0.0001, 0.9999), 0, 0.2))

par(mfrow = c(2, 2))
hist(c1a)
hist(c1b)
hist(c2a)
hist(c2b)

hrSeq = seq(from = 0.01, to = 0.8, by = 0.005)
hrSeq2 = 1 - exp(log(1 - hrSeq))