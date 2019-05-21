#******************************************************************************
# skewedModels.R
# Date revised: Oct 3 2018
# Inputs: stock-recruit data csv; stan models
# Outputs: estimates of stock recruit and supplementary parameters
# Explainer: Quantifies evidence of skewness and heavy-tailedness using models
# stan models. These estimates will be used to inform forward simulation runs
# in synchrony analysis using recoverySim. R Based on code in 
# synchRetroAnlaysis-stan.Rmd written by S. Anderson.
#******************************************************************************


library(dplyr)
library(ggplot2)
library(rstan)
library(here)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
# Specify models
sm_lm <- stan_model(here("scripts/lm.stan"))
sm_lmt <- stan_model(here("scripts/lm-t.stan"))
sm_lmtskew <- stan_model(here("scripts/lm-t-skew.stan"))

recDatTrim1 <- readRDS(here("data", "generated", "recDatTrim1.rds"))

# Function to fit model across stocks
fit_sr_stan <- function(dat,
                        model_type = c("lm", "lmt", "lmtskew"),
                        sr_type = c("ricker", "larkin"),
                        chains = 4, iter = 2000, cores = 1, adapt_delta = 0.95) {
  
  sr_type <- match.arg(sr_type)
  model_type <- match.arg(model_type)
  
  if (sr_type == "ricker")
    X_ij <- model.matrix(~ets, data = dat)
  else
    X_ij <- model.matrix(~ets + ets1 + ets2 + ets3, data = dat)
  
  model <- switch(model_type,
                  lm = sm_lm,
                  lmt = sm_lmt,
                  lmtskew = sm_lmtskew
  )
  
  sampling(model,
           data = list(X_ij = X_ij, y_i = dat$logProd, N = nrow(X_ij), J = ncol(X_ij)),
           iter = iter, chains = chains, cores = cores,
           control = list(adapt_delta = 0.99, max_treedepth = 20))
}

# Apply each model to each stock
out <- plyr::dlply(recDatTrim1, "stk", function(d) {
  fit_sr_stan(d, model_type = "lm", sr_type = unique(d$model), chains = 4,
  iter = 2000)
})
out_t <- plyr::dlply(recDatTrim1, "stk", function(d) {
  fit_sr_stan(d, model_type = "lmt", sr_type = unique(d$model), chains = 4,
              iter = 2000)
})
out_tskew <- plyr::dlply(recDatTrim1, "stk", function(d) {
  fit_sr_stan(d, model_type = "lmtskew", sr_type = unique(d$model), chains = 4,
              iter = 2000, adapt_delta = 0.99)
})

# Extract par estimates
post_n <- plyr::ldply(out, function(x) {
  e <- extract(x)
  sm_summ <- summary(x)$summary
  max_rhat <- max(sm_summ[, "Rhat"])
  min_neff <- min(sm_summ[, "n_eff"])
  data.frame(alpha = e$b_j[, 1], beta = e$b_j[, 2], max_rhat, min_neff)
})
saveRDS(post_n, file = here("data", "generated", "normModelOut.rds"))

post_t <- plyr::ldply(out_t, function(x) {
  e <- extract(x)
  sm_summ <- summary(x)$summary
  max_rhat <- max(sm_summ[, "Rhat"])
  min_neff <- min(sm_summ[, "n_eff"])
  data.frame(alpha = e$b_j[, 1], beta = e$b_j[, 2], nu = e$nu, max_rhat, min_neff)
})
saveRDS(post_t, file = here("data", "generated", "studentTModelOut.rds"))

e <- extract(out_t[[19]])

post_tskew <- plyr::ldply(out_tskew, function(x) {
  e <- extract(x)
  sm_summ <- summary(x)$summary
  max_rhat <- max(sm_summ[, "Rhat"])
  min_neff <- min(sm_summ[, "n_eff"])
  data.frame(alpha = e$b_j[, 1], beta = e$b_j[, 2], log_skew = e$log_skew, 
             nu = e$nu, max_rhat, min_neff)
})
saveRDS(post_tskew, file = here("data", "generated", "studentTSkewModelOut.rds"))


# Check diagnostics
post_t %>% 
  select(stk, max_rhat, min_neff) %>% 
  unique() %>% 
  knitr::kable(digits = c(0, 2, 1), row.names = FALSE)
post_tskew %>% 
  select(stk, max_rhat, min_neff) %>%
  unique() %>% 
  knitr::kable(digits = c(0, 2, 1), row.names = FALSE)

# Plot posterior estimates
ggplot(post_t, aes(as.factor(stk), nu)) + #nu from non-skewed model
  geom_violin() +
  geom_hline(yintercept = c(2, 10), lty = 2) +
  coord_flip(ylim = c(2, 65))
ggplot(post_t, aes(as.factor(stk), beta)) + #skewness 
  geom_violin() +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip()

ggplot(post_tskew, aes(as.factor(stk), alpha)) + #nu from skewed model
  geom_violin() +
  coord_flip()
ggplot(post_tskew, aes(as.factor(stk), beta)) + #skewness 
  geom_violin() +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip()
ggplot(post_tskew, aes(as.factor(stk), nu)) + #nu from skewed model
  geom_violin() +
  geom_hline(yintercept = c(2, 10), lty = 2) + 
  coord_flip(ylim = c(2, 65))
ggplot(post_tskew, aes(as.factor(stk), exp(log_skew))) + #skewness 
  geom_violin() +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip()
#averaged (equal weighting) posteriors
ggplot(post_tskew, aes(exp(log_skew))) + geom_density() +
  geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = mean(exp(post_tskew$log_skew)), col = "red") +
  geom_vline(xintercept = median(exp(post_tskew$log_skew)), col = "blue") +
  geom_vline(xintercept = quantile(exp(post_tskew$log_skew), probs = 0.1), col = "blue") +
  geom_vline(xintercept = quantile(exp(post_tskew$log_skew), probs = 0.9), col = "blue")
#Probability nu < 10 
filter(post_t, stk == 13) %>% summarise(p = mean(nu < 10)) %>% pull(p)
filter(post_t, stk == 3) %>% summarise(p = mean(nu < 10)) %>% pull(p)
filter(post_t, stk == 12) %>% summarise(p = mean(nu < 10)) %>% pull(p)

#quantiles of parameters
mean(exp(post_tskew$log_skew))
median(exp(post_tskew$log_skew))
quantile(exp(post_tskew$log_skew), probs = c(0.1, 0.25, 0.75, 0.9))


## Use log(0.67) as skewness parameter, approximately equivalent to 75th percentile of
# skewness distribution


## Compare alpha and beta estimates from skewed vs. skewed t
post_n <- saveRDS(post_n, file = here("data", "generated", "normModelOut.rds"))
post_t <- readRDS(file = here("data", "generated", "studentTModelOut.rds"))
post_tskew <- readRDS(file = here("data", "generated", "studentTSkewModelOut.rds"))


parsT <- post_t %>% 
  select(stk, alpha, beta) %>% 
  mutate(model = "t")
parsSkewT <- post_tskew %>% 
  select(stk, alpha, beta) %>% 
  mutate(model = "skewT")
modPars <- rbind(parsT, parsSkewT)
ggplot(modPars, aes(x = as.factor(model), y = beta, fill = model)) +
  geom_boxplot() +
  facet_wrap(~as.factor(stk))
