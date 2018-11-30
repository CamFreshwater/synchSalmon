#******************************************************************************
# distributionPlots.R
# Date revised: Oct 7 2018
# Outputs: figure showing i) distributions used in different productivity OMs 
# for synch analysis and ii) empty plot used to represent TAM rule
#******************************************************************************

require(here); require(ggplot2); require(viridis); require(sn); 
require(tidyverse); require(samSim); require(ggpubr)

plotDat <- data.frame(x = seq(from = -5, to = 5, length.out = 9999))
plotDat <- plotDat %>% 
  mutate(norm = dst(plotDat$x, xi = 0, alpha = 0, nu = Inf,
                     omega = 0.75),
         skewN = dst(plotDat$x, xi = 0, alpha = log(0.67), nu = Inf,
                     omega = 0.75),
         skewT = dst(plotDat$x, xi = 0, alpha = log(0.67), nu = 2,
                        omega = 0.75),
         yMin = 0
         ) %>%
  gather(key = dist, value = probY, c(-x, -yMin)) %>%
  mutate(dist = as.factor(dist)) %>% 
  mutate(dist = factor(dist, levels(dist)[c(3, 2, 1)]),)


colPal <- c("#bd0026", "#fd8d3c", "black")
p <- ggplot() +
  geom_ribbon(data = plotDat,
              aes(x = x, ymin = yMin, ymax = probY, group = dist, color = dist,
                  fill = dist),
               position = "identity", alpha = 0.2) +
  labs(x = "Standard Deviations", y = "Probability Density") +
  geom_vline(xintercept = 0, color = "grey30", linetype = 2, size = 1.1) +
  scale_color_manual(name = "Distribution", values = colPal) +
  scale_fill_manual(name = "Distribution", values = colPal) +
  geom_rect(aes(xmin = -4, xmax = -2, ymin = 0, ymax = 0.2), lty = 2,
            fill = "transparent", color = "black") +
  guides(fill = FALSE, color = FALSE) +
  theme_sleekX(axisSize = 12)

q <- ggplot() +
  geom_ribbon(data = plotDat,
              aes(x = x, ymin = yMin, ymax = probY, group = dist, color = dist,
                  fill = dist),
              position = "identity", alpha = 0.2) +
  scale_x_continuous(limits = c(-4, -2)) +
  scale_y_continuous(limits = c(0, 0.2)) +
  labs(x = "Standard Deviations", y = "Probabiliy Density") +
  scale_color_manual(name = "Distribution", values = colPal,
                     labels = c("norm" = "Normal",
                                "skewN" = "Skew\nNormal",
                                "skewT" = "Skew\nStudent t")) +
  scale_fill_manual(name = "Distribution", values = colPal,
                    labels = c("norm" = "Normal",
                               "skewN" = "Skew\nNormal",
                               "skewT" = "Skew\nStudent t")) +
  guides(fill = FALSE, color = FALSE) +
  theme_sleekX(axisSize = 12)


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
  mutate(prodOM = factor(prodOM, levels(prodOM)[c(2, 1, 3, 4)]),
         dist = case_when(
           prodOM == "normRec" | prodOM == "normLowARec" ~ "Normal",
           prodOM == "skewNormRec" ~ "Skew Normal",
           prodOM == "skewTRec" ~ "Skew Student t"
         )) %>%
  mutate(prodOM = recode(prodOM, "normRec" = "Reference Alpha", 
                         "normLowARec" = "Low Alpha", 
                         "skewNormRec" = "Skew Normal",
                         "skewTRec" = "Skew Student t"),
         ) %>% 
  filter(recruits < recCap) #remove absurdly large positive recruitment events

r <- ggplot(recDat, aes(x = prodOM, y = recruits)) + 
  geom_violin(trim = FALSE, aes(color = dist, fill = dist), alpha = 0.2) +
  stat_summary(fun.data = "mean_sd", geom = "pointrange", 
               color = "black", size = 1.25) +
  labs(x = "Productivity OM", y = "Return Abundance") +
  scale_color_manual(name = "Distribution", values = colPal) +
  scale_fill_manual(name = "Distribution", values = colPal) +
  # guides(fill = FALSE, color = FALSE) +
  theme_sleekX(legendSize = 1.1, axisSize = 12)


png(file = paste(here(),"/figs/SFig1_distPlots.png", sep = ""),
    height = 6.5, width = 8.5, units = "in", res = 300)
ggarrange(
  ggarrange(p, q, ncol = 2, labels = c("a)", "b)"), widths = c(1,1),
            font.label = list(size = 13, face = "plain")),
          r, nrow = 2, labels = c("", "c)"), 
  font.label = list(size = 13, face = "plain")
  )
dev.off()
    

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
