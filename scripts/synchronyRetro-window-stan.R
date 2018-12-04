library(dplyr)
library(ggplot2)
library(rstan)
library(here)
rstan_options(auto_write = TRUE)

X <- readRDS(here("outputs", "generatedData", "recMat.rds"))
# X <- readRDS(here("data", "generated", "recMatShort.rds"))

sm <- stan_model(here("scripts", "sync.stan"))

# var(rowSums(X))/(sum(apply(X, 2, sd))^2)

fit_synch <- function(x) {
  m <- sampling(sm,
    data = list(
      N = nrow(x), y_time = rowSums(x),
      G = ncol(x), y_ind = t(x)),
    iter = 500, chains = 1
  )
  m
}


m <- fit_synch(X)
phi_hat <- extract(m)$phi
range(phi_hat)

m <- list()
window <- 12
for (i in seq(window, nrow(X))) {
  m[[i - (window-1)]] <- fit_synch(X[(i - (window-1)):i, , drop = FALSE])
}

phi <- plyr::ldply(m, function(.x) extract(.x)$phi)
apply(phi, 1, median)
apply(phi, 1, quantile, probs = 0.05)
apply(phi, 1, quantile, probs = 0.95)

# Stopping for now!
