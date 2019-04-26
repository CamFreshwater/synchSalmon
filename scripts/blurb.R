require(here); require(tidyverse); require(ggplot2)
cuPar <- read.csv(here("data/sox/fraserCUpars.csv"), stringsAsFactors = F)

enRouteSig <- cuPar$sdDBE
enRouteMR <- cuPar$meanDBE
nCU <- length(unique(cuPar$stk))

# migMortErr <- exp(qnorm(runif(nCU, 0.0001, 0.9999), 0, enRouteSig))
qnorm(runif(nCU, 0.0001, 0.9999), 0, enRouteSig)
enRouteMR * migMortErr
