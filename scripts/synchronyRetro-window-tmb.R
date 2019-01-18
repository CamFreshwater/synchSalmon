library(dplyr)
library(ggplot2)
library(here)
library(TMB)

X <- readRDS(here("outputs", "generatedData", "recMat.rds"))
# X <- readRDS(here("outputs", "generatedData", "recMatShort.rds"))

compile("scripts/sync.cpp")
dyn.load(dynlib("scripts/sync"))

# Test our model:
obj <- TMB::MakeADFun(
  data = list(y_time = rowSums(X), y_ind = X),
  parameters = list(
    log_numerator_sigma = 0, log_denominator_sigma = rep(0, ncol(X)),
    numerator_mean = 0, denominator_mean = rep(0, ncol(X))
  ),
  DLL = "sync"
)
object <- stats::nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
rep <- TMB::sdreport(obj)
rep

synchrony <- function(x) {
  var(rowSums(x)) / (sum(apply(x, 2, sd))^2)
}
assertthat::are_equal(
  synchrony(X),
  obj$report()$phi,
  tol = 0.00001
)
obj$report()
summary(rep)

# A little function to make it easy to fit the windows of years:
fit_synch <- function(x) {
  obj <- TMB::MakeADFun(
    data = list(y_time = rowSums(x), y_ind = x),
    parameters = list(
      log_numerator_sigma = 0, log_denominator_sigma = rep(0, ncol(x)),
      numerator_mean = 0, denominator_mean = rep(0, ncol(x))
    ),
    DLL = "sync"
  )
  object <- stats::nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
  rep <- TMB::sdreport(obj)
  d <- as.data.frame(summary(rep))
  d$term <- row.names(d)
  d <- d[d$term %in% c("logit_phi", "log_cv_s", "log_cv_c"),
    c("Estimate", "Std. Error", "term")]
  row.names(d) <- NULL
  d
}

# test it:
fit_synch(X)

# Now fitted to the year windows:
out_list <- list()
window <- 12 # just as desired
for (i in seq(window, nrow(X)))
  out_list[[i - (window - 1)]] <- fit_synch(X[(i - (window - 1)):i, , drop = FALSE])
out <- purrr::map_df(out_list, as.data.frame)
out$year <- rep(seq(window, nrow(X)), each = 3)
out <- arrange(out, term, year) %>%
  rename(se_link = `Std. Error`, est_link = Estimate) %>%
  mutate(logit = grepl("logit", term))

# The model returns estimates in logit/log space.
# Add and subtract the standard errors in logit space and
# inverse logit transform:
out <- mutate(out,
  est = ifelse(logit, plogis(est_link), exp(est_link)),
  lwr = ifelse(logit, plogis(est_link + qnorm(0.025) * se_link),
    exp(est_link + qnorm(0.025) * se_link)),
  upr = ifelse(logit, plogis(est_link + qnorm(0.975) * se_link),
    exp(est_link + qnorm(0.975) * se_link))
)
saveRDS(out, here("outputs", "generatedData", "tmbSynchEst.rds"))
# saveRDS(out, here("outputs", "generatedData", "tmbSynchEst_suppShort.rds"))
# Note that those are 90% CIs. Adjust as desired.

out %>%
  mutate(term = gsub("logit_", "", term)) %>%
  mutate(term = gsub("log_", "", term)) %>%
  ggplot(aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  facet_wrap(~term, scales = "free_y")
