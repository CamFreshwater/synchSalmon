#******************************************************************************
# distributionPlots.R
# Date revised: Oct 7 2018
# Outputs: figure showing i) distributions used in different productivity OMs 
# for synch analysis and ii) empty plot used to represent TAM rule
#******************************************************************************

require(here); require(ggplot2); require(viridis); require(sn); 
require(tidyverse); require(samSim); require(ggpubr)

dat <- data.frame(norm = rst(n = 5000000, xi = 0, alpha = 0, nu = Inf,
                             omega = 0.75),
                  sNorm = rst(n = 5000000, xi = 0, alpha = log(0.67), nu = Inf,
                              omega = 0.75),
                  sT = rst(n = 5000000, xi = 0, alpha = log(0.67), nu = 2,
                           omega = 0.75)) %>%
  gather(key = dist, value = values) %>%
  mutate(dist = as.factor(dist),
         expDev = exp(values))

plotDat <- data.frame(x = c(-5, 5))
colPal <- c("black", "#fd8d3c", "#bd0026")
p <- ggplot(plotDat, aes(x = x)) +
  # geom_density(aes(group = dist, color = dist, fill = dist), 
  #              position = "identity", alpha = 0.1) +
  # stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) +
  stat_function(fun = dst, n = 101, 
                args = list(xi = 0, alpha = 0, nu = Inf, omega = 0.75),
                aes(colour = colPal[1], fill = colPal[1])) +
  stat_function(fun = dst, n = 101, 
                args = list(xi = 0, alpha = log(0.67), nu = 1e6, omega = 0.75),
                aes(colour = colPal[2], fill = colPal[2])) +
  stat_function(fun = dst, n = 101, 
                args = list(xi = 0, alpha = log(0.67), nu = 2, omega = 0.75),
                aes(colour = colPal[3], fill = colPal[3])) +
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
  guides(fill = FALSE, color = FALSE) +
  theme_sleekX(legendSize = 0.65)

datTrim <- dat %>%
  filter(values < -2,
         values > -4)
q <- ggplot(datTrim, aes(x = values)) +
  geom_density(aes(group = dist, color = dist, fill = dist), 
               position = "identity", alpha = 0.1) +
  scale_x_continuous(limits = c(-4, -2)) +
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


### Generate simulated data for boxplots
# Default pars for Chilko
n = 1000
chlkA = 1.8315
chlkLowA = 1.58795
chlkB = 1.23182
meanChlkS = 0.3792771
chlkSig = 0.79075
recCap = 14.4 #cap on maximum recruitment

#Generate estimates of recruitment based on Chilko parameters and different
#deviation distributions
normRec <- rickerModel(meanChlkS, chlkA, chlkB, 
                       error = rst(n, xi = 0, omega = chlkSig, alpha = log(1), 
                                   nu = Inf), 
                       rho = 0, utminus = 0)[[1]]
normLowARec <- rickerModel(meanChlkS, chlkLowA, chlkB, 
                           error = rst(n, xi = 0, omega = chlkSig, 
                                       alpha = log(1), nu = Inf), 
                           rho = 0, utminus = 0)[[1]]
skewNormRec <- rickerModel(meanChlkS, chlkA, chlkB, 
                           error = rst(n, xi = 0, omega = chlkSig, 
                                       alpha = log(0.65), nu = Inf), 
                           rho = 0, utminus = 0)[[1]]
skewTRec <- rickerModel(meanChlkS, chlkA, chlkB, 
                        error = rst(n, xi = 0, omega = chlkSig, 
                                    alpha = log(0.65), nu = 3), 
                        rho = 0, utminus = 0)[[1]]

recDat <- cbind(normRec, normLowARec, skewNormRec, skewTRec) %>% 
  as.data.frame() %>% 
  gather(key = prodOM, value = recruits) %>% 
  mutate(prodOM = as.factor(prodOM)) %>% 
  mutate(prodOM = factor(prodOM, levels(prodOM)[c(2, 1, 3, 4)])) %>%
  mutate(prodOM = recode(prodOM, "normRec" = "Normal", 
                         "normLowARec" = "Low Alpha", 
                         "skewNormRec" = "Skewed Normal",
                         "skewTRec" = "Skewed Student t")) %>% 
  filter(recruits < recCap) #remove absurdly large positive recruitment events


r <- ggplot(recDat, aes(x = prodOM, y = recruits)) + 
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = "mean_sdl", geom = "pointrange", 
               color = "black") +
  labs(x = "Productivity OM", y = "Return Abundance") +
  theme_sleekX(legendSize = 0.65)


png(file = paste(here(),"/figs/Fig5_catchGroupedPlots_3OMs.png", sep = ""),
    height = 5.5, width = 8.5, units = "in", res = 300)
ggarrange(ggarrange(p, q, ncol = 2, labels = c("a)", "b)"), widths = c(1,1.25)),
          r, nrow = 2, labels = c("", "c)"))
  

### How often do extreme events occur?
nDist <- rst(n = 100000, xi = 0, alpha = 0, nu = Inf, omega = 1)
sNDist <- rst(n = 100000, xi = 0, alpha = log(0.67), nu = Inf, omega = 1)
sTDist <- rst(n = 100000, xi = 0, alpha = log(0.67), nu = 3, omega = 1)

1 / (length(nDist[nDist < -3]) / length(nDist))
1 / (length(sNDist[sNDist < -3]) / length(sNDist))
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
