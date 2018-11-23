#******************************************************************************
# distributionPlots.R
# Date revised: Oct 7 2018
# Outputs: figure showing i) distributions used in different productivity OMs for 
# synch analysis and ii) empty plot used to represent TAM rule
#******************************************************************************

require(here); require(ggplot2); require(viridis); require(sn); 
require(tidyverse)
source(here("scripts/synchFunctions.R"))

dat <- data.frame(norm = rst(n = 100000, xi = 0, alpha = 0, nu = Inf, 
                             omega = 0.75),
                  sNorm = rst(n = 100000, xi = 0, alpha = log(0.67), nu = Inf, 
                              omega = 0.75),
                  sT = rst(n = 100000, xi = 0, alpha = log(0.67), nu = 2, 
                           omega = 0.75)) %>% 
  gather(key = dist, value = values) %>% 
  mutate(dist = as.factor(dist),
         expDev = exp(values))


colPal <- c("black", "#fd8d3c", "#bd0026")
p <- ggplot(dat, aes(x = values)) +
  geom_density(aes(group = dist, color = dist, fill = dist), 
               position = "identity", alpha = 0.1) +
  scale_x_continuous(limits = c(-5, 5)) +
  labs(x = "Standard Deviations", y = "Probabiliy Density") +
  geom_vline(xintercept = 0, color = "grey30", linetype = 2) +
  scale_color_manual(name = "Distribution", values = colPal,
                     labels = c("norm" = "Normal",
                                "sNorm" = "Skewed Normal",
                                "sT" = "Skewed Student t")) +
  scale_fill_manual(name = "Distribution", values = colPal,
                    labels = c("norm" = "Normal",
                               "sNorm" = "Skewed Normal",
                               "sT" = "Skewed Student t")) +
  theme_sleekX(legendSize = 0.65)

png(here("figs/distFig.png"), height = 3, width = 5, res = 300, units = "in")
print(p)
dev.off()


### How often do extreme events occur?
nDist <- rst(n = 100000, xi = 0, alpha = 0, nu = Inf, omega = 1)
sTDist <- rst(n = 100000, xi = 0, alpha = log(0.67), nu = 3, omega = 1)

1 / (length(nDist[nDist < -3]) / length(nDist))
1 / (length(sTDist[sTDist < -3]) / length(sTDist))


### Make empty plot to illustrate TAM rules

# Generate data
ret = seq(0, 2, length.out = 999)
esc = rep(NA, length.out = 999)
hr = rep(NA, length.out = 999)
for (i in 1:length(ret)) {
  if (ret[i] < 0.4) {
    esc[i] <- ret[i]
    hr[i] <- 0.05
  }
  if (ret[i] > 0.4 & ret[i] < 1) {
    esc[i] <- 0.4
    hr[i] <- max(0.05, (ret[i] - 0.4) / ret[i])
  }
  if (ret[i] > 1) {
    hr[i] <- 0.6
    esc[i] <- ret[i] * 0.4
  }
}

png(here("figs/emptyTAM.png"), height = 6, width = 6, res = 600, units = "in")
par(mfrow = c(2, 1), mar = c(3,4,1.5,0.5), mgp = c(3, 0.75, 0))
plot(hr ~ ret, xlim=c(0, 2), ylim=c(0, 1), ylab = "", xlab = "", axes = F, 
     type ="l", lwd = 2)
axis(1, las = 1)
axis(2, las = 1)
mtext("Total Allowable Mortalty\n(TAM) Rate", 2, line = 2)
abline(v = 0.4, col = "red", lty = 2, lwd = 2)
abline(v = 1, col = "darkgreen", lty = 2, lwd = 2)
text(x = c(0.15, 0.7, 1.5), y = rep(0.9, length.out = 3), labels = 
       c("Zone 1", "Zone 2", "Zone 3"))
plot(esc ~ ret, xlim=c(0, 2), ylim=c(0, 1), ylab="", xlab="", axes = F, 
     type ="l", lwd = 2)
axis(1, las = 1)
axis(2, las = 1)
abline(v = 0.4, col = "red", lty = 2, lwd = 2)
abline(v = 1, col = "darkgreen", lty = 2, lwd = 2)
mtext("Escapement\n(millions)", 2, line = 2)
mtext("Return Abundance (millions)", 1, line = 2)
text(x = c(0.15, 0.7, 1.5), y = rep(0.9, length.out = 3), labels = 
       c("Zone 1", "Zone 2", "Zone 3"))
dev.off()
