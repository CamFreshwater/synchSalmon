library(dplyr)
library(ggplot2)
library(rstan)
library(here)
rstan_options(auto_write = TRUE)

X <- readRDS(here("data", "generated", "prodMat.rds"))
X_raw <- X

y_group <- rowSums(X)
y_ind <- X

sm <- stan_model(here("scripts", "sync.stan"))

m <- sampling(sm,
  data = list(
    N = nrow(X), y_group = y_group,
    G = ncol(X), y_ind = t(y_ind)
  )
)

# rstan::traceplot(m)

sqrt_phi <- sqrt(extract(m)$phi)
cv_s <- extract(m)$cv_s
cv_c <- extract(m)$cv_c

add_label <- function(xfrac, yfrac, label, pos = 2, ...) {
  u <- par("usr")
  x <- u[2] - xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

G <- ncol(X)
pal <- rep(RColorBrewer::brewer.pal(8, "Dark2"), 999L)[seq_len(ncol(X))]
# pdf("thibaut2013-eg4.pdf", width = 5, height = 7)
par(mfrow = c(4, 1), mar = c(3.5, 3.5, 0, 0), oma = c(2, 1, 1, 1),
  mgp = c(1.6, 0.4, 0), tck = -0.02, las = 1)
matplot(X_raw, type = "l", col = pal, xlab = "Time", ylab = "Productivity",
  lty = 1, xaxs = "i")
plot(density(sqrt_phi), xlim = c(0, 1), xlab = "sqrt(phi)", main = "", xaxs = "i")
add_label(0.01, 0.1, "Synchrony component")
plot(density(cv_s), xlim = c(0, 2), xlab = "CV_s", main = "", xaxs = "i")
add_label(0.01, 0.1, "Component variability")
for (i in seq_len(G)) {
  polygon(density(extract(m)$cv_s_ind[, i]), border = paste0(pal[i], "60"))
}
plot(density(cv_c), xlim = c(0, 1), xlab = "CV_c", main = "", xaxs = "i")
add_label(0.01, 0.1, "Realized aggregate variability")
# hist(extract(m)$cv_s * sqrt(extract(m)$phi))
# dev.off()

m <- list()
X <- X / sd(as.matrix(X)) # scale
y_group <- rowSums(X)
y_ind <- X

window <- 10

for (i in seq(window, nrow(X))) {
  y_group <- rowSums(X[(i - (window-1)):i, ])
  y_ind <- X[(i - (window-1)):i, ]
  m_temp <- sampling(sm,
    data = list(
      N = length(y_group), y_group = y_group,
      G = ncol(X), y_ind = t(y_ind)
    ),
    iter = 1000, chains = 1
  )
  m[[i]] <- extract(m_temp)
}

out_cv_s <- plyr::ldply(m, function(y)
  quantile(y$cv_s, probs = c(0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9)))

plot_ts <- function(dat, column, ylim = c(0, 1), ylab = column, xlab = "Year") {
  out <- plyr::ldply(dat, function(y) quantile(y[[column]],
    probs = c(0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9)))
  xt <- rev(nrow(X) - seq_len(nrow(out)))
  plot(xt, sqrt(out[, "50%"]), type = "l", ylim = ylim, ylab = ylab,
    lwd = 1.5, yaxs = "i", xlab = xlab)
  polygon(c(xt, rev(xt)), c(sqrt(out[, "40%"]), rev(sqrt(out[, "60%"]))),
    border = NA, col = "#00000050")
  polygon(c(xt, rev(xt)), c(sqrt(out[, "25%"]), rev(sqrt(out[, "75%"]))),
    border = NA, col = "#00000040")
  polygon(c(xt, rev(xt)), c(sqrt(out[, "10%"]), rev(sqrt(out[, "90%"]))),
    border = NA, col = "#00000025")
}

# pdf("fraser-synchrony-trends.pdf", width = 3.5, height = 7)
par(mfrow = c(3, 1), mar = c(2, 3, 0, 0), oma = c(2, 1, 1, 1),
  mgp = c(2.0, 0.5, 0), tck = -0.02, las = 1)
plot_ts(m, "cv_s", ylim = c(0, 1.6), ylab = "Weighted stream-level CV", xlab = "")
plot_ts(m, "phi", ylab = "Synchrony", xlab = "")
plot_ts(m, "cv_c", ylim = c(0, 1.5), ylab = "Fraser CV")
# dev.off()
